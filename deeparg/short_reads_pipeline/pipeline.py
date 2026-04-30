import logging
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from deeparg.modern.api import DEFAULT_MODEL_REPO, predict as deeparg_predict

logger = logging.getLogger(__name__)


@dataclass
class ShortReadsResult:
    clean_fasta: str
    deeparg_prefix: str
    arg_file: str
    potential_arg_file: str
    arg_quant_file: str
    normalized_subtype_file: str = None
    normalized_type_file: str = None


def _quote_command(command):
    return " ".join(str(part) for part in command)


def _run(command, stdout_path=None):
    logger.info("Running: %s", _quote_command(command))
    if stdout_path:
        with open(stdout_path, "wb") as stdout:
            subprocess.run(command, check=True, stdout=stdout)
    else:
        subprocess.run(command, check=True)


def _require_tool(name):
    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            "{} is required for this step but was not found on PATH. "
            "Install it with bioconda/conda-forge or pass an already-clean FASTA "
            "with --input-fasta.".format(name)
        )
    return path


def _copy_input_fasta(input_fasta, output_file):
    clean_fasta = output_file + ".clean"
    if os.path.abspath(input_fasta) != os.path.abspath(clean_fasta):
        shutil.copyfile(input_fasta, clean_fasta)
    return clean_fasta


def trim_paired_reads(forward_fastq, reverse_fastq):
    trimmomatic = _require_tool("trimmomatic")
    outputs = {
        "forward_paired": forward_fastq + ".paired",
        "forward_unpaired": forward_fastq + ".unpaired",
        "reverse_paired": reverse_fastq + ".paired",
        "reverse_unpaired": reverse_fastq + ".unpaired",
    }
    _run([
        trimmomatic,
        "PE",
        "-phred33",
        forward_fastq,
        reverse_fastq,
        outputs["forward_paired"],
        outputs["forward_unpaired"],
        outputs["reverse_paired"],
        outputs["reverse_unpaired"],
        "LEADING:3",
        "TRAILING:3",
        "SLIDINGWINDOW:4:15",
        "MINLEN:36",
    ])
    return outputs


def merge_paired_reads(forward_fastaq, reverse_fastaq, output_file):
    vsearch = _require_tool("vsearch")
    merged = forward_fastaq + ".merged"
    unmerged_forward = forward_fastaq + ".unmerged"
    unmerged_reverse = reverse_fastaq + ".unmerged"
    clean_fasta = output_file + ".clean"
    _run([
        vsearch,
        "--fastq_mergepairs",
        forward_fastaq,
        "--reverse",
        reverse_fastaq,
        "--fastaout",
        merged,
        "--fastaout_notmerged_fwd",
        unmerged_forward,
        "--fastaout_notmerged_rev",
        unmerged_reverse,
    ])
    with open(clean_fasta, "wb") as output:
        for path in (merged, unmerged_forward, unmerged_reverse):
            with open(path, "rb") as handle:
                shutil.copyfileobj(handle, output)
    return clean_fasta


def _hf_file(model_repo, filename):
    local_base = Path(model_repo)
    if local_base.exists():
        candidate = local_base / filename
        if not candidate.exists():
            raise FileNotFoundError(str(candidate))
        return str(candidate)

    from deeparg.entry import _force_ipv4_if_requested

    _force_ipv4_if_requested()
    from huggingface_hub import hf_hub_download

    return hf_hub_download(repo_id=model_repo, filename=filename, repo_type="model")


def resolve_arg_gene_lengths(data_path=None, model_repo=DEFAULT_MODEL_REPO):
    if data_path:
        candidate = Path(data_path) / "database" / "v2" / "features.gene.length"
        if candidate.exists():
            return str(candidate)
        candidate = Path(data_path) / "SS" / "database" / "features.gene.length"
        if candidate.exists():
            return str(candidate)
    return _hf_file(model_repo, "SS/database/features.gene.length")


def _read_arg_gene_lengths(path):
    gene_lengths = {}
    with open(path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 2:
                continue
            gene = parts[0].split("|")[-1].upper()
            gene_lengths[gene] = float(parts[1])
    return gene_lengths


def _parse_arg_intervals(arg_file):
    intervals = []
    with open(arg_file) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 12:
                continue
            start = int(float(row[1]))
            end = int(float(row[2]))
            if end < start:
                start, end = end, start
            intervals.append({
                "gene": row[0].upper(),
                "start": start,
                "end": end,
                "count": int(float(row[11])),
                "arg_class": row[4],
            })
    return sorted(intervals, key=lambda item: (item["gene"], item["start"], item["end"]))


def merge_arg_intervals(arg_file, gene_lengths_path):
    intervals = _parse_arg_intervals(arg_file)
    merged_path = arg_file + ".merged"
    quant_path = merged_path + ".quant"
    gene_lengths = _read_arg_gene_lengths(gene_lengths_path)

    merged = []
    current = None
    for item in intervals:
        if current and item["gene"] == current["gene"] and item["start"] <= current["end"]:
            current["end"] = max(current["end"], item["end"])
            current["count"] += item["count"]
            current["classes"].add(item["arg_class"])
        else:
            if current:
                merged.append(current)
            current = {
                "gene": item["gene"],
                "start": item["start"],
                "end": item["end"],
                "count": item["count"],
                "classes": {item["arg_class"]},
            }
    if current:
        merged.append(current)

    genes = {}
    with open(merged_path, "w") as handle:
        for item in merged:
            classes = "/".join(sorted(item["classes"]))
            handle.write("\t".join([
                item["gene"],
                str(item["start"]),
                str(item["end"]),
                str(item["count"]),
                classes,
            ]) + "\n")
            gene = genes.setdefault(item["gene"], {"count": 0, "length": 0, "classes": set()})
            gene["count"] += item["count"]
            gene["length"] += abs(item["end"] - item["start"])
            gene["classes"].update(item["classes"])

    with open(quant_path, "w") as handle:
        for gene, item in sorted(genes.items()):
            gene_len = gene_lengths.get(gene)
            if not gene_len:
                logger.warning("No gene length found for %s; skipping quantification row.", gene)
                continue
            cov = item["length"] / float(gene_len)
            handle.write("\t".join([
                gene,
                "/".join(sorted(item["classes"])),
                str(item["count"]),
                str(item["length"]),
                str(int(gene_len) if gene_len.is_integer() else gene_len),
                str(round(cov, 3)),
            ]) + "\n")
    return quant_path


def _resolve_gg13_prefix(data_path=None, gg13_prefix=None, model_repo=DEFAULT_MODEL_REPO):
    if gg13_prefix:
        return gg13_prefix
    if data_path:
        candidate = Path(data_path) / "gg13" / "dataset"
        if (candidate.with_suffix(".len")).exists():
            return str(candidate)

    local_base = Path(model_repo)
    if local_base.exists():
        candidate = local_base / "gg13" / "dataset"
        if (candidate.with_suffix(".len")).exists():
            return str(candidate)
        raise FileNotFoundError(str(candidate.with_suffix(".len")))

    try:
        from deeparg.entry import _force_ipv4_if_requested

        _force_ipv4_if_requested()
        from huggingface_hub import snapshot_download

        snapshot_dir = snapshot_download(
            repo_id=model_repo,
            repo_type="model",
            allow_patterns=["gg13/*"],
        )
    except Exception as exc:
        raise RuntimeError(
            "16S normalization requires the Greengenes bowtie2 index. "
            "Pass --gg13-prefix, provide --deeparg-data-path with gg13/dataset, "
            "or upload gg13/* to the HF model repo. Original error: {}".format(exc)
        )
    return str(Path(snapshot_dir) / "gg13" / "dataset")


def _cigar_reference_length(cigar):
    length = 0
    for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in {"M", "D", "N", "=", "X"}:
            length += int(size)
    return length


def _merge_intervals(rows):
    rows = sorted(rows, key=lambda item: (item[0], item[1], item[2]))
    merged = []
    current = None
    for name, start, end, count in rows:
        if current and name == current[0] and start <= current[2]:
            current[2] = max(current[2], end)
            current[3] += count
        else:
            if current:
                merged.append(tuple(current))
            current = [name, start, end, count]
    if current:
        merged.append(tuple(current))
    return merged


def quantify_16s(clean_fasta, gg13_prefix):
    bowtie2 = _require_tool("bowtie2")
    samtools = _require_tool("samtools")
    sam_path = clean_fasta + ".sam"
    bam_path = clean_fasta + ".bam"
    sorted_bam_path = clean_fasta + ".sorted.bam"
    merged_path = sorted_bam_path + ".merged"
    quant_path = merged_path + ".quant"

    _run([bowtie2, "-f", "--fast-local", "--no-unal", "-x", gg13_prefix, "-U", clean_fasta, "-S", sam_path])
    _run([samtools, "view", "-bS", sam_path], stdout_path=bam_path)
    _run([samtools, "sort", bam_path, "-o", sorted_bam_path])

    view = subprocess.run([samtools, "view", sorted_bam_path], check=True, text=True, stdout=subprocess.PIPE)
    intervals = []
    for line in view.stdout.splitlines():
        row = line.split("\t")
        if len(row) < 6 or row[2] == "*":
            continue
        start = int(row[3]) - 1
        ref_len = _cigar_reference_length(row[5])
        if ref_len <= 0:
            continue
        intervals.append((row[2], start, start + ref_len, 1))

    gene_lengths = {}
    with open(gg13_prefix + ".len") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                gene_lengths[parts[0]] = float(parts[1])

    genes = {}
    with open(merged_path, "w") as handle:
        for name, start, end, count in _merge_intervals(intervals):
            handle.write("\t".join([name, str(start), str(end), str(count)]) + "\n")
            gene = genes.setdefault(name, {"count": 0, "length": 0})
            gene["count"] += count
            gene["length"] += abs(end - start)

    with open(quant_path, "w") as handle:
        for gene, item in sorted(genes.items()):
            gene_len = gene_lengths.get(gene)
            if not gene_len:
                continue
            cov = item["length"] / float(gene_len)
            handle.write("\t".join([
                gene,
                str(item["count"]),
                str(item["length"]),
                str(int(gene_len) if gene_len.is_integer() else gene_len),
                str(round(cov, 3)),
            ]) + "\n")
    return quant_path


def normalize_to_16s(quant_16s_file, arg_quant_file, coverage_fraction):
    n16s = sum(
        int(line.split()[1])
        for line in open(quant_16s_file)
        if line.split() and float(line.split()[2]) >= 100
    )
    logger.info("Total number of 16S reads in the sample: %s", n16s)
    l16s = 1432
    subtype_path = arg_quant_file + ".subtype"
    type_path = arg_quant_file + ".type"
    arg_types = {}

    with open(subtype_path, "w") as subtype_out:
        subtype_out.write("#ARG-group\tReadCount\t16s-NormalizedReadCount\n")
        with open(arg_quant_file) as handle:
            for line in handle:
                parts = line.split()
                if len(parts) < 6:
                    continue
                subtype, arg_type, count, _alg_len, gene_len, cov = parts[:6]
                if float(cov) <= coverage_fraction:
                    continue
                normalized = (int(count) / float(gene_len)) / (float(n16s) / l16s) if n16s > 0 else 0
                subtype_out.write("\t".join([subtype, count, str(normalized)]) + "\n")
                current = arg_types.setdefault(arg_type, [0, gene_len])
                current[0] += int(count)

    with open(type_path, "w") as type_out:
        type_out.write("#ARG-category\tReadCount\t16s-NormalizedReadCount\n")
        for arg_type, (count, gene_len) in sorted(arg_types.items()):
            normalized = (count / float(gene_len)) / (float(n16s) / l16s) if n16s > 0 else 0
            type_out.write("\t".join([arg_type, str(count), str(normalized)]) + "\n")

    return subtype_path, type_path


def run_short_reads_pipeline(
    output_file,
    forward_pe_file=None,
    reverse_pe_file=None,
    input_fasta=None,
    data_path=None,
    model_repo=DEFAULT_MODEL_REPO,
    diamond_path=None,
    identity=80,
    probability=0.8,
    evalue=1e-10,
    gene_coverage=1,
    normalize_16s=False,
    gg13_prefix=None,
    skip_trimming=False,
    device=None,
    chunk_size=10000,
):
    output_parent = Path(output_file).parent
    if str(output_parent) not in ("", "."):
        output_parent.mkdir(parents=True, exist_ok=True)

    if input_fasta:
        clean_fasta = _copy_input_fasta(input_fasta, output_file)
    else:
        if not forward_pe_file or not reverse_pe_file:
            raise ValueError("Provide either --input-fasta or both paired-end FASTQ files.")
        if skip_trimming:
            trimmed = {"forward_paired": forward_pe_file, "reverse_paired": reverse_pe_file}
        else:
            logger.info("Step 1: trimming and QC with Trimmomatic")
            trimmed = trim_paired_reads(forward_pe_file, reverse_pe_file)
        logger.info("Step 2: merging paired-end reads with vsearch")
        clean_fasta = merge_paired_reads(trimmed["forward_paired"], trimmed["reverse_paired"], output_file)

    logger.info("Step 3: running DeepARG-SS")
    deeparg_prefix = clean_fasta + ".deeparg"
    prediction_outputs = deeparg_predict(
        input_file=clean_fasta,
        output_file=deeparg_prefix,
        model="SS",
        input_type="nucl",
        model_repo=model_repo,
        diamond_path=diamond_path,
        data_path=data_path,
        min_prob=probability,
        identity=identity,
        evalue=evalue,
        overlap=0.8,
        device=device,
        chunk_size=chunk_size,
    )

    logger.info("Step 4: quantifying ARG-like reads")
    gene_lengths = resolve_arg_gene_lengths(data_path=data_path, model_repo=model_repo)
    arg_quant_file = merge_arg_intervals(prediction_outputs["arg"], gene_lengths)

    subtype_path = None
    type_path = None
    if normalize_16s:
        logger.info("Step 5: normalizing ARG counts to 16S rRNAs")
        resolved_gg13 = _resolve_gg13_prefix(data_path=data_path, gg13_prefix=gg13_prefix, model_repo=model_repo)
        quant_16s = quantify_16s(clean_fasta, resolved_gg13)
        subtype_path, type_path = normalize_to_16s(
            quant_16s,
            arg_quant_file,
            float(gene_coverage) / 100.0,
        )

    return ShortReadsResult(
        clean_fasta=clean_fasta,
        deeparg_prefix=deeparg_prefix,
        arg_file=prediction_outputs["arg"],
        potential_arg_file=prediction_outputs["potential_arg"],
        arg_quant_file=arg_quant_file,
        normalized_subtype_file=subtype_path,
        normalized_type_file=type_path,
    )


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description="Run the modern DeepARG short reads pipeline.")
    parser.add_argument("--forward-pe-file", "--forward_pe_file", dest="forward_pe_file")
    parser.add_argument("--reverse-pe-file", "--reverse_pe_file", dest="reverse_pe_file")
    parser.add_argument("--input-fasta", dest="input_fasta", help="Already-clean FASTA; skips Trimmomatic and vsearch.")
    parser.add_argument("--output-file", "--output_file", dest="output_file", required=True)
    parser.add_argument("-d", "--deeparg-data-path", "--deeparg_data_path", dest="data_path")
    parser.add_argument("--hf-model-path", default=os.environ.get("DEEPARG_HF_MODEL", DEFAULT_MODEL_REPO))
    parser.add_argument("--diamond-path")
    parser.add_argument("--deeparg-identity", "--deeparg_identity", dest="identity", default=80, type=float)
    parser.add_argument("--deeparg-probability", "--deeparg_probability", dest="probability", default=0.8, type=float)
    parser.add_argument("--deeparg-evalue", "--deeparg_evalue", dest="evalue", default=1e-10, type=float)
    parser.add_argument("--gene-coverage", "--gene_coverage", dest="gene_coverage", default=1, type=float)
    parser.add_argument("--normalize-16s", action="store_true", help="Run bowtie2/samtools 16S normalization.")
    parser.add_argument("--gg13-prefix", help="Bowtie2 index prefix for Greengenes gg13 dataset.")
    parser.add_argument("--skip-trimming", action="store_true", help="Use paired FASTQ files directly as vsearch input.")
    parser.add_argument("--hf-device", default=None)
    parser.add_argument("--hf-chunk-size", default=10000, type=int)
    args = parser.parse_args(argv)

    result = run_short_reads_pipeline(
        output_file=args.output_file,
        forward_pe_file=args.forward_pe_file,
        reverse_pe_file=args.reverse_pe_file,
        input_fasta=args.input_fasta,
        data_path=args.data_path,
        model_repo=args.hf_model_path,
        diamond_path=args.diamond_path,
        identity=args.identity,
        probability=args.probability,
        evalue=args.evalue,
        gene_coverage=args.gene_coverage,
        normalize_16s=args.normalize_16s,
        gg13_prefix=args.gg13_prefix,
        skip_trimming=args.skip_trimming,
        device=args.hf_device,
        chunk_size=args.hf_chunk_size,
    )
    print(result)
    return 0
