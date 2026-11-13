# tSTR revisions: reconstructed T sequence comparison

Code for generating two FASTA files to compare the A/T content of reference sequences that pass or fail the reconstructed prediction-error filter. This repository contains code, documentation, and the two reconstructed reference-sequence FASTA files below. Prediction tables are not included.

## Sequence downloads and generated outputs

- `data/T_retained.fa.gz`: 175,581 sequences, absolute error ≤ 0.05.
- `data/T_excluded.fa.gz`: 571,161 sequences, absolute error > 0.05.
- `data/T_scores.csv.gz` (generated, not included): ID, expression, scaled rank, prediction, error, group and A/T fraction.
- `data/manifest.json` (generated, not included): counts, input integrity checks, output sizes and SHA-256 hashes.

These are 101-base model input windows, not full 1,001-base source sequences. IDs preserve the original motif/strand labels. Selection is by **oriented T model class**, so an original A-labelled sequence can belong here. The window is source positions 450–550 inclusive (Python `[450:551]`); no additional reverse-complement operation is applied.

## Download

Download the gzip-compressed FASTA files directly:

- [T_retained.fa.gz](https://github.com/mathysgrapotte/tstr-revisions/raw/refs/heads/main/data/T_retained.fa.gz): 5,048,865 bytes; 175,581 sequences; SHA-256 `336785456bc4f61fb25115f443ac70a3e3a8b8c999e3fc00d01c5cd646bf7dcb`.
- [T_excluded.fa.gz](https://github.com/mathysgrapotte/tstr-revisions/raw/refs/heads/main/data/T_excluded.fa.gz): 16,043,430 bytes; 571,161 sequences; SHA-256 `f6aabcd9656e14d2f5bf36c18a77f354e635cb3a07a60cc6c136b73b3c8e6d60`.

These files are **reconstructed hg38 reference model windows**, not the recovered original retained/excluded GTEx sequences. Download and decompress with `gzip -d`, or use software that reads `.fa.gz` directly.

## Important limits

This is the reconstructed **hg38 reference population**, not a recovered copy of the original filtered GTEx association population or the manuscript's 22,491-STR motif population. No individual genotypes or GTEx variance filter are included. It can answer a reference-sequence composition comparison for this reconstruction; it cannot by itself establish the original analysis's exact counts or validate the claimed absence of A/T bias. Do not call these the original retained/excluded sequences.

Observed expression values are converted to average ranks within the T class, with ties sharing an average rank. Those ranks are scaled so the smallest is 0 and the largest 1. Each sequence's absolute prediction error determines its group. This uses the reconstructed convention with strongest agreement to published retained IDs across classes: 8,913/9,036, not an exact reproduction. The residual 123 differences remain unresolved.

## Recreate

Python 3, standard library only. Supply the reconstructed predictions CSV and original public data archive:

```sh
python3 export_t_sequences.py --predictions reference_predictions.csv --archive data.zip --out data
python3 -m unittest -v test_export.py
```

The output directory must not exist (use a new directory name rather than the supplied `data` directory when regenerating). If a run fails, do not use its partial output. The script validates the decompressed reference entry against its recorded CRC32 and size and checks that every prediction ID maps to exactly one valid window. It can alternatively consume the existing 150-chunk archive cache with `--chunks DIRECTORY`.

## Provenance

Public source archive landing page: https://seafile.lirmm.fr/f/2b74469daf5542149d1a/

Reference entry: `data/hg38.hipstr_reference.cage.500bp.around3end.fa`

Expected decompressed bytes: 3,179,903,896. CRC32: 3,712,841,224.

Reference predictions were reconstructed separately using model code at commit `d265a955fd09889c5b0580e6c740c3546a2b7cab` of https://gite.lirmm.fr/ml4reggen/rstr (see `RunModelPredictionsOnVcf/bin/launch_make_predictions.py`). The prediction CSV hash is recorded in the manifest. This export script does not rerun model inference.

The generated outputs contain only reference sequences and reconstructed reference predictions. No controlled-access individual genotypes are used by this export script. The two reference-sequence FASTA files are included for the sequence-composition comparison. They do not contain individual genotypes or the original GTEx association inputs.
