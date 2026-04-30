# DeepARG v2 Head-to-Head Findings

## Overall Status

PASS semantic equivalence for both LS and SS.

The final files are not byte-identical because NumPy and PyTorch print tiny
floating-point probability differences, but all labels, rows, and non-probability
fields match. All output probability differences are within `1e-5`.

## Method

Both models were compared on the repository example FASTA, `test/ORFs.fa`.
For each model:

1. DIAMOND 2.1.24 generated an alignment TSV through the modern HF-backed CLI.
2. The legacy side loaded original Zenodo `metadata_*.pkl` and `model_*.pkl`.
3. A NumPy forward pass reproduced the original Lasagne dense network using the
   legacy matrices directly.
4. The modern side loaded the converted PyTorch/Transformers model from
   `data/hf/deeparg-v2/<MODEL>`.
5. Prediction rows and final `.ARG` / `.potential.ARG` files were compared.

This avoids needing Python 2/Theano while still testing the exact legacy weights
against the converted weights.

## Model-Level Results

| Model | Status | Raw alignment rows | Samples after filters | Legacy prediction rows | Modern prediction rows | Top-label mismatches | Prediction-row mismatches @1e-5 | Max probability diff |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LS | PASS | 2174 | 141 | 169 | 169 | 0 | 0 | 1.162290573e-06 |
| SS | PASS | 2174 | 241 | 246 | 246 | 0 | 0 | 3.576278687e-07 |

## Output File Results

| Output | Legacy lines | Modern lines | Byte-identical | Non-probability mismatches | Max probability diff | Within 1e-5 |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| LS `.ARG` | 115 | 115 | False | 0 | 2.384185791e-07 | True |
| LS `.potential.ARG` | 56 | 56 | False | 0 | 1.162290573e-06 | True |
| SS `.ARG` | 237 | 237 | False | 0 | 2.384185791e-07 | True |
| SS `.potential.ARG` | 11 | 11 | False | 0 | 3.576278687e-07 | True |

## Findings

* LS and SS both reproduce legacy labels and prediction rows exactly at the
  decision level.
* No ARG class changed in either model.
* No best-hit, query coordinate, identity, alignment length, bitscore, e-value,
  or count field changed.
* Differences are limited to floating-point probability formatting/noise from
  NumPy vs PyTorch execution.
* The largest observed model probability difference is below `1.2e-6`, well
  below the report tolerance of `1e-5`.

## Artifacts

* `reports/head_to_head_v2/combined_summary.json`
* `reports/head_to_head_v2/ls_summary.json`
* `reports/head_to_head_v2/ss_summary.json`
* `reports/head_to_head_v2/ls_legacy_numpy.ARG`
* `reports/head_to_head_v2/ls_modern_transformers.ARG`
* `reports/head_to_head_v2/ss_legacy_numpy.ARG`
* `reports/head_to_head_v2/ss_modern_transformers.ARG`
* matching `.potential.ARG`, label mismatch, and prediction mismatch JSON files

## Conclusion

The converted PyTorch/Transformers DeepARG v2 models are semantically equivalent
to the legacy Zenodo PKL models on the repository example data for both LS and
SS. Based on this test, the modern runtime is suitable to become the default
runtime before removing heavy local artifacts from the repository.
