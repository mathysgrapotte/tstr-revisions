#!/usr/bin/env python3
"""Export reconstructed T reference windows. Standard library only; no genotype data."""
import argparse, collections, csv, gzip, hashlib, io, json, pathlib, struct, zlib, zipfile

def scaled_ranks(values):
    counts = collections.Counter(values)
    ranks, before = {}, 0
    for value in sorted(counts):
        n = counts[value]
        ranks[value] = before + (n - 1) / 2
        before += n
    lo, hi = min(ranks.values()), max(ranks.values())
    if hi == lo:
        raise ValueError('Cannot scale a constant target')
    return {v: (r - lo) / (hi - lo) for v, r in ranks.items()}

def source_blocks(args):
    if args.archive:
        with zipfile.ZipFile(args.archive) as archive:
            with archive.open('data/hg38.hipstr_reference.cage.500bp.around3end.fa') as source:
                while block := source.read(4 * 1024 * 1024):
                    yield block
    else:
        decoder = zlib.decompressobj(-15)
        for k in range(150):
            block = (args.chunks / f'chunk_{k:03}.bin').read_bytes()
            if k == 0:
                if block[63:67] != b'PK\x03\x04':
                    raise ValueError('Unexpected cached ZIP layout')
                n, e = struct.unpack_from('<HH', block, 89)
                block = block[93+n+e:]
            yield decoder.decompress(block)
        if not decoder.eof:
            raise ValueError('Incomplete compressed source')

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--predictions', type=pathlib.Path, required=True)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--archive', type=pathlib.Path)
    source.add_argument('--chunks', type=pathlib.Path, help='Existing verified 150-chunk reconstruction cache')
    p.add_argument('--out', type=pathlib.Path, required=True)
    args = p.parse_args()
    with args.predictions.open() as f:
        rows = [r for r in csv.DictReader(f) if r['oriented_class'] == 'T']
    ranks = scaled_ranks(float(r['raw_cage_tag']) for r in rows)
    lookup = {r['id']: r for r in rows}
    if len(lookup) != len(rows): raise ValueError('Duplicate prediction IDs')
    args.out.mkdir(parents=True, exist_ok=False)
    counts = collections.Counter(); seen = set(); crc = size = 0
    handles = {group: gzip.open(args.out / f'T_{group}.fa.gz', 'wt') for group in ('retained', 'excluded')}
    scores = gzip.open(args.out / 'T_scores.csv.gz', 'wt', newline='')
    writer = csv.writer(scores)
    writer.writerow(['id','raw_cage_tag','scaled_average_rank','model_prediction','absolute_error','group','AT_fraction'])
    def finish(header, sequence):
        if header is None: return
        ident = header[1:].split('|')[0]
        row = lookup.get(ident)
        if row is None: return
        if ident in seen: raise ValueError(f'Duplicate source ID: {ident}')
        window = ''.join(sequence)[450:551].upper()
        if len(window) != 101 or set(window) - set('ACGT'): raise ValueError(f'Invalid window: {ident}')
        target = ranks[float(row['raw_cage_tag'])]
        error = abs(float(row['model_prediction']) - target)
        group = 'retained' if error <= .05 else 'excluded'
        handles[group].write(f'>{ident}\n{window}\n')
        writer.writerow([ident,row['raw_cage_tag'],target,row['model_prediction'],error,group,(window.count('A')+window.count('T'))/101])
        seen.add(ident); counts[group] += 1
    buf = b''; header = None; sequence = []
    try:
        for block in source_blocks(args):
            crc = zlib.crc32(block, crc); size += len(block)
            lines = (buf + block).split(b'\n'); buf = lines.pop()
            for raw in lines:
                line = raw.decode().strip()
                if line.startswith('>'):
                    finish(header, sequence); header = line; sequence = []
                else: sequence.append(line)
        if buf: sequence.append(buf.decode().strip())
        finish(header, sequence)
    finally:
        for handle in handles.values(): handle.close()
        scores.close()
    if (crc,size) != (3712841224,3179903896): raise ValueError(f'Source integrity mismatch: {crc}, {size}')
    if seen != set(lookup): raise ValueError('Missing source sequences')
    manifest = {'population':'Reconstructed hg38 oriented T reference class; not original GTEx/motif eligible population', 'counts':dict(counts),'total':len(seen),'source_crc32':crc,'source_bytes':size,'window_zero_based_half_open':[450,551],'filter':'absolute error <= 0.05 against within-class scaled average ranks','variance_filter_applied':False,'files':{}}
    for path in sorted(args.out.iterdir()):
        manifest['files'][path.name] = {'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    manifest['predictions_sha256'] = hashlib.sha256(args.predictions.read_bytes()).hexdigest()
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__ == '__main__': main()
