#!/usr/bin/env python3
"""Export reconstructed T reference windows. Standard library only; no genotype data."""

import argparse
import collections
import csv
import gzip
import hashlib
import json
import pathlib
import struct
import zipfile
import zlib
from typing import Iterator

# --- Source archive -------------------------------------------------------

REFERENCE_ENTRY = 'data/hg38.hipstr_reference.cage.500bp.around3end.fa'
EXPECTED_SOURCE_CRC32 = 3712841224
EXPECTED_SOURCE_BYTES = 3179903896
READ_BLOCK_BYTES = 4 * 1024 * 1024

# Layout of the cached 150-chunk copy of the archive: chunk_000 begins with a
# 63-byte prefix followed by the ZIP local file header of the reference entry.
CHUNK_COUNT = 150
LOCAL_HEADER_OFFSET = 63
LOCAL_HEADER_SIGNATURE = b'PK\x03\x04'
LOCAL_HEADER_SIZE = 30
NAME_LENGTHS_OFFSET = LOCAL_HEADER_OFFSET + 26  # file name + extra field lengths

# --- Selection and filtering ----------------------------------------------

TARGET_CLASS = 'T'
WINDOW_START = 450  # zero-based, inclusive
WINDOW_END = 551  # zero-based, exclusive
WINDOW_LENGTH = WINDOW_END - WINDOW_START
VALID_BASES = set('ACGT')
MAX_RETAINED_ERROR = 0.05
GROUPS = ('retained', 'excluded')

SCORE_COLUMNS = [
    'id',
    'raw_cage_tag',
    'scaled_average_rank',
    'model_prediction',
    'absolute_error',
    'group',
    'AT_fraction',
]


def scaled_ranks(values):
    """Map each distinct value to its average rank, scaled to [0, 1].

    Tied values share the mean of the ranks they occupy. The smallest rank is
    scaled to 0 and the largest to 1.
    """
    counts = collections.Counter(values)
    ranks = {}
    preceding = 0
    for value in sorted(counts):
        n = counts[value]
        ranks[value] = preceding + (n - 1) / 2
        preceding += n

    lowest, highest = min(ranks.values()), max(ranks.values())
    if highest == lowest:
        raise ValueError('Cannot scale a constant target')
    return {value: (rank - lowest) / (highest - lowest) for value, rank in ranks.items()}


# --- Reading the reference FASTA ------------------------------------------


def read_archive_blocks(archive_path: pathlib.Path) -> Iterator[bytes]:
    """Yield the decompressed reference entry from the original ZIP archive."""
    with zipfile.ZipFile(archive_path) as archive:
        with archive.open(REFERENCE_ENTRY) as source:
            while block := source.read(READ_BLOCK_BYTES):
                yield block


def read_chunk_blocks(chunk_dir: pathlib.Path) -> Iterator[bytes]:
    """Yield the decompressed reference entry from the 150-chunk archive cache."""
    decoder = zlib.decompressobj(-zlib.MAX_WBITS)  # raw DEFLATE stream
    for index in range(CHUNK_COUNT):
        block = (chunk_dir / f'chunk_{index:03}.bin').read_bytes()
        if index == 0:
            block = strip_local_file_header(block)
        yield decoder.decompress(block)
    if not decoder.eof:
        raise ValueError('Incomplete compressed source')


def strip_local_file_header(block: bytes) -> bytes:
    """Drop the prefix and ZIP local file header so only DEFLATE data remains."""
    signature_end = LOCAL_HEADER_OFFSET + len(LOCAL_HEADER_SIGNATURE)
    if block[LOCAL_HEADER_OFFSET:signature_end] != LOCAL_HEADER_SIGNATURE:
        raise ValueError('Unexpected cached ZIP layout')
    name_length, extra_length = struct.unpack_from('<HH', block, NAME_LENGTHS_OFFSET)
    data_start = LOCAL_HEADER_OFFSET + LOCAL_HEADER_SIZE + name_length + extra_length
    return block[data_start:]


class IntegrityTracker:
    """Accumulate CRC32 and byte count over a stream of blocks."""

    def __init__(self):
        self.crc32 = 0
        self.size = 0

    def track(self, blocks: Iterator[bytes]) -> Iterator[bytes]:
        for block in blocks:
            self.crc32 = zlib.crc32(block, self.crc32)
            self.size += len(block)
            yield block

    def verify(self):
        if (self.crc32, self.size) != (EXPECTED_SOURCE_CRC32, EXPECTED_SOURCE_BYTES):
            raise ValueError(f'Source integrity mismatch: {self.crc32}, {self.size}')


def iter_lines(blocks: Iterator[bytes]) -> Iterator[str]:
    """Split a stream of byte blocks into stripped text lines."""
    remainder = b''
    for block in blocks:
        lines = (remainder + block).split(b'\n')
        remainder = lines.pop()
        for line in lines:
            yield line.decode().strip()
    if remainder:
        yield remainder.decode().strip()


def iter_fasta_records(lines: Iterator[str]) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs; the header keeps its leading '>'."""
    header = None
    sequence_parts = []
    for line in lines:
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(sequence_parts)
            header = line
            sequence_parts = []
        else:
            sequence_parts.append(line)
    if header is not None:
        yield header, ''.join(sequence_parts)


def record_id(header: str) -> str:
    """Return the ID part of a FASTA header such as '>ID|extra'."""
    return header[1:].split('|')[0]


def extract_window(sequence: str, ident: str) -> str:
    """Return the validated, upper-cased model input window."""
    window = sequence[WINDOW_START:WINDOW_END].upper()
    if len(window) != WINDOW_LENGTH or set(window) - VALID_BASES:
        raise ValueError(f'Invalid window: {ident}')
    return window


def at_fraction(window: str) -> float:
    return (window.count('A') + window.count('T')) / WINDOW_LENGTH


# --- Predictions ----------------------------------------------------------


def load_target_predictions(path: pathlib.Path):
    """Load T-class predictions keyed by ID, plus scaled ranks of their expression."""
    with path.open() as f:
        rows = [row for row in csv.DictReader(f) if row['oriented_class'] == TARGET_CLASS]
    ranks = scaled_ranks(float(row['raw_cage_tag']) for row in rows)
    predictions = {row['id']: row for row in rows}
    if len(predictions) != len(rows):
        raise ValueError('Duplicate prediction IDs')
    return predictions, ranks


# --- Export ---------------------------------------------------------------


def export_windows(records, predictions, ranks, out_dir: pathlib.Path):
    """Write grouped FASTA files and the score table; return (group counts, seen IDs)."""
    counts = collections.Counter()
    seen = set()
    fasta_files = {group: gzip.open(out_dir / f'T_{group}.fa.gz', 'wt') for group in GROUPS}
    scores_file = gzip.open(out_dir / 'T_scores.csv.gz', 'wt', newline='')
    try:
        scores = csv.writer(scores_file)
        scores.writerow(SCORE_COLUMNS)
        for header, sequence in records:
            ident = record_id(header)
            row = predictions.get(ident)
            if row is None:
                continue
            if ident in seen:
                raise ValueError(f'Duplicate source ID: {ident}')

            window = extract_window(sequence, ident)
            target = ranks[float(row['raw_cage_tag'])]
            error = abs(float(row['model_prediction']) - target)
            group = 'retained' if error <= MAX_RETAINED_ERROR else 'excluded'

            fasta_files[group].write(f'>{ident}\n{window}\n')
            scores.writerow([
                ident,
                row['raw_cage_tag'],
                target,
                row['model_prediction'],
                error,
                group,
                at_fraction(window),
            ])
            seen.add(ident)
            counts[group] += 1
    finally:
        for handle in fasta_files.values():
            handle.close()
        scores_file.close()
    return counts, seen


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as f:
        while block := f.read(READ_BLOCK_BYTES):
            digest.update(block)
    return digest.hexdigest()


def build_manifest(counts, total, source, out_dir: pathlib.Path, predictions_path: pathlib.Path):
    files = {
        path.name: {'bytes': path.stat().st_size, 'sha256': sha256_file(path)}
        for path in sorted(out_dir.iterdir())
    }
    return {
        'population': (
            'Reconstructed hg38 oriented T reference class; '
            'not original GTEx/motif eligible population'
        ),
        'counts': dict(counts),
        'total': total,
        'source_crc32': source.crc32,
        'source_bytes': source.size,
        'window_zero_based_half_open': [WINDOW_START, WINDOW_END],
        'filter': 'absolute error <= 0.05 against within-class scaled average ranks',
        'variance_filter_applied': False,
        'files': files,
        'predictions_sha256': sha256_file(predictions_path),
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions', type=pathlib.Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--archive', type=pathlib.Path)
    source.add_argument(
        '--chunks', type=pathlib.Path, help='Existing verified 150-chunk reconstruction cache'
    )
    parser.add_argument('--out', type=pathlib.Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    predictions, ranks = load_target_predictions(args.predictions)
    args.out.mkdir(parents=True, exist_ok=False)

    blocks = read_archive_blocks(args.archive) if args.archive else read_chunk_blocks(args.chunks)
    source = IntegrityTracker()
    records = iter_fasta_records(iter_lines(source.track(blocks)))
    counts, seen = export_windows(records, predictions, ranks, args.out)

    source.verify()
    if seen != set(predictions):
        raise ValueError('Missing source sequences')

    manifest = build_manifest(counts, len(seen), source, args.out, args.predictions)
    manifest_text = json.dumps(manifest, indent=2)
    (args.out / 'manifest.json').write_text(manifest_text + '\n')
    print(manifest_text)


if __name__ == '__main__':
    main()
