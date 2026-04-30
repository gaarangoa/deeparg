# Short Reads Pipeline Modernization

The legacy `deeparg short_reads_pipeline` was inspected and reimplemented on top
of the new PyTorch/Transformers runtime.

## Legacy Behavior Preserved

The original paired-end workflow did:

1. Trimmomatic PE trimming.
2. vsearch paired-end merge.
3. `deeparg predict --model SS`.
4. ARG interval merge and per-gene quantification.
5. Optional 16S normalization with Bowtie2, samtools, and the gg13 index.

## Modern Implementation

The new implementation lives in `deeparg/short_reads_pipeline/pipeline.py`.

Changes:

* Calls `deeparg.modern.api.predict(..., model="SS")` instead of shelling out
  to the legacy Python 2 CLI.
* Uses `subprocess.run(..., check=True)` for external tools.
* Replaces the ARG `bedtools merge` step with native Python interval merging.
* Supports `--input-fasta` for the simplest path when reads are already clean.
* Keeps paired FASTQ mode with Trimmomatic and vsearch when those tools are on
  `PATH`.
* Keeps optional 16S normalization and downloads `gg13/dataset.*` from
  `gaarangoa/deeparg` when needed.

## Hugging Face Runtime Assets

The HF repo `gaarangoa/deeparg` now contains:

* LS and SS Transformers model weights.
* LS and SS DeepARG v2 DIAMOND databases.
* DIAMOND 2.1.24 binaries for Linux x86_64 and Linux aarch64.
* `gg13/dataset.*` Bowtie2 index files for 16S normalization.
* Modern pipeline source and reports.

External command-line tools still expected from the system for raw FASTQ and
16S workflows:

* `trimmomatic`
* `vsearch`
* `bowtie2`
* `samtools`

The clean-FASTA short reads path only needs DeepARG/HF assets plus DIAMOND.
