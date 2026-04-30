import argparse
import sys
import logging
import os
import shlex
import shutil
import subprocess
import platform
from pathlib import Path

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)


def _quote_cmd(command):
    return " ".join(shlex.quote(str(part)) for part in command)


def _run_command(command):
    logger.info("Running: {}".format(_quote_cmd(command)))
    subprocess.run(command, check=True)


def _test_diamond(candidate):
    try:
        result = subprocess.run(
            [candidate, "version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
        )
    except Exception as exc:
        return False, str(exc)
    if result.returncode != 0:
        return False, (result.stderr or result.stdout).strip()
    return True, (result.stdout or result.stderr).strip()


def resolve_diamond(args):
    candidates = []
    if args.diamond_path:
        candidates.append(args.diamond_path)

    path_diamond = shutil.which("diamond")
    if path_diamond:
        candidates.append(path_diamond)

    if args.data_path:
        bundled = os.path.join(args.data_path, "bin", "diamond")
        candidates.append(bundled)

    if args.hf_model_path:
        hf_diamond = resolve_hf_diamond(args.hf_model_path)
        if hf_diamond:
            candidates.append(hf_diamond)

    seen = set()
    errors = []
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        if not os.path.exists(candidate) and os.path.sep in candidate:
            errors.append("{}: not found".format(candidate))
            continue
        ok, detail = _test_diamond(candidate)
        if ok:
            logger.info("Using DIAMOND: {} ({})".format(candidate, detail))
            return candidate
        errors.append("{}: {}".format(candidate, detail))

    message = [
        "No compatible DIAMOND executable found.",
        "Install a native DIAMOND build, for example:",
        "  conda install -c conda-forge -c bioconda diamond=2.1.24",
        "or pass --diamond-path /path/to/diamond.",
    ]
    if errors:
        message.append("Tried:")
        message.extend("  " + error for error in errors)
    raise RuntimeError("\n".join(message))


def platform_diamond_dir():
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system != "linux":
        return None
    if machine in ("x86_64", "amd64"):
        return "linux-x86_64"
    if machine in ("aarch64", "arm64"):
        return "linux-aarch64"
    return None


def resolve_hf_diamond(model_path):
    diamond_dir = platform_diamond_dir()
    if not diamond_dir:
        return None
    filename = "bin/{}/diamond".format(diamond_dir)
    try:
        path = _hf_download(model_path, None, filename)
    except Exception as exc:
        logger.info("Could not download bundled DIAMOND from HF: {}".format(exc))
        return None
    try:
        os.chmod(path, os.stat(path).st_mode | 0o755)
    except OSError:
        pass
    return path


def _strip_dmnd_suffix(path):
    path = str(path)
    if path.endswith(".dmnd"):
        return path[:-5]
    return path


def _force_ipv4_if_requested():
    if os.environ.get("HF_FORCE_IPV4") != "1":
        return
    import socket

    original_getaddrinfo = socket.getaddrinfo

    def getaddrinfo_ipv4(*args, **kwargs):
        return [
            item for item in original_getaddrinfo(*args, **kwargs)
            if item[0] == socket.AF_INET
        ]

    socket.getaddrinfo = getaddrinfo_ipv4


def _hf_download(model_path, subfolder, filename):
    local_base = Path(model_path)
    if local_base.exists():
        local_path = local_base
        if subfolder:
            local_path = local_path / subfolder.strip("/")
        local_path = local_path / filename
        if not local_path.exists():
            raise FileNotFoundError(str(local_path))
        return str(local_path)

    _force_ipv4_if_requested()
    from huggingface_hub import hf_hub_download

    remote_path = filename
    if subfolder:
        remote_path = "{}/{}".format(subfolder.strip("/"), filename)
    return hf_hub_download(repo_id=model_path, filename=remote_path, repo_type="model")


def resolve_database(args, pipeline):
    if args.data_path:
        version_dir = os.path.join(args.data_path, "database", args.model_version)
        return {
            "diamond_db": os.path.join(version_dir, "features"),
            "gene_lengths": os.path.join(version_dir, "features.gene.length"),
        }

    if not args.hf_model_path:
        raise RuntimeError("--data-path is required unless --hf-model-path points to a DeepARG HF bundle.")

    subfolder = args.hf_model_subfolder or args.model
    dmnd = _hf_download(args.hf_model_path, subfolder, "database/features.dmnd")
    gene_lengths = None
    if pipeline == "genes":
        gene_lengths = _hf_download(args.hf_model_path, subfolder, "database/features.gene.length")
    return {
        "diamond_db": _strip_dmnd_suffix(dmnd),
        "gene_lengths": gene_lengths,
    }


def predict(args):

    if args.model == 'SS':
        mdl = '_SS'
        pipeline = 'reads'
    else:
        mdl = '_LS'
        pipeline = 'genes'

    if args.type == "prot":
        args.aligner = "blastp"
    else:
        args.aligner = "blastx"   

    if args.hf_model_path and not args.hf_model_subfolder:
        args.hf_model_subfolder = args.model

    diamond = resolve_diamond(args)
    database = resolve_database(args, pipeline)

    logger.info("DIAMOND {} alignment".format(args.aligner))
    cmd = [
        diamond,
        args.aligner,
        '-q', args.input_file,
        '-d', database["diamond_db"],
        '-k', str(args.arg_num_alignments_per_entry),
        '--id', str(args.arg_alignment_identity),
        '--sensitive',
        '-e', str(args.arg_alignment_evalue),
        '-a', args.output_file+'.align',
    ]
    _run_command(cmd)

    logger.info("Input:{} output:{} model: deepARG{}, Input type: {}".format(args.input_file, args.output_file, mdl, args.aligner)) 

    cmd = [
        diamond,
        'view',
        '-a', args.output_file+'.align.daa',
        '-o', args.output_file+'.align.daa.tsv'
    ]
    logger.info("parsing output file {}".format(cmd))
    _run_command(cmd)

    from deeparg.modern.predict import predict_file

    predict_file(
        args.output_file+'.align.daa.tsv',
        args.output_file+'.mapping',
        args.hf_model_path,
        args.arg_alignment_identity,
        args.arg_alignment_evalue,
        args.min_prob,
        args.arg_alignment_overlap,
        pipeline,
        gene_lengths_path=database["gene_lengths"],
        chunk_size=args.hf_chunk_size,
        device=args.hf_device,
        subfolder=args.hf_model_subfolder,
    )


def short_reads_pipeline(args):
    from deeparg.short_reads_pipeline import run_short_reads_pipeline

    result = run_short_reads_pipeline(
        output_file=args.output_file,
        forward_pe_file=args.forward_pe_file,
        reverse_pe_file=args.reverse_pe_file,
        input_fasta=args.input_fasta,
        data_path=args.deeparg_data_path,
        model_repo=args.hf_model_path,
        diamond_path=args.diamond_path,
        identity=args.deeparg_identity,
        probability=args.deeparg_probability,
        evalue=args.deeparg_evalue,
        gene_coverage=args.gene_coverage,
        normalize_16s=args.normalize_16s,
        gg13_prefix=args.gg13_prefix,
        skip_trimming=args.skip_trimming,
        device=args.hf_device,
        chunk_size=args.hf_chunk_size,
    )
    logger.info("Short reads outputs: %s", result)
    return 0


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # use deeparg section
    reads = subparsers.add_parser("predict", help="Predict ARG from reads or genes")
    reads.add_argument('--model', required=True,
                       help='Select model to use (short sequences for reads | long sequences for genes) SS|LS [No default]')
    reads.add_argument('-i', '--input-file', required=False,
                       help='Input file (Fasta input file)')
    reads.add_argument('-o', '--output-file', required=True,
                       help='Output file where to store results')
    reads.add_argument('-d', '--data-path', required=False,
                       help="Optional local DeepARG data path. If omitted with --hf-model-path, database files are downloaded from Hugging Face.")
    reads.add_argument('--type', default='nucl',
                       help='Molecular data type prot/nucl [Default: nucl]')
    reads.add_argument('--min-prob', default=0.8, type=float,
                       help='Minimum probability cutoff [Default: 0.8]')
    reads.add_argument('--arg-alignment-identity', default=50, type=float,
                       help='Identity cutoff for sequence alignment in percent [Default: 50]')
    reads.add_argument('--arg-alignment-evalue', default=1e-10, type=float,
                       help='Evalue cutoff [Default: 1e-10]')
    reads.add_argument('--arg-alignment-overlap', default=0.8, type=float,
                       help='Alignment overlap cutoff between read and genes [Default: 0.8]')
    reads.add_argument('--arg-num-alignments-per-entry', default=1000, type=int,
                       help='Diamond, minimum number of alignments per entry [Default: 1000]')
    reads.add_argument('--model-version', default='v2',
                       help='Model deepARG version  [Default: v2]')
    reads.add_argument('--diamond-path', default=None,
                       help='Optional path to a native DIAMOND executable. If omitted, deeparg searches PATH and data-path/bin/diamond.')
    reads.add_argument('--hf-model-path', default=os.environ.get("DEEPARG_HF_MODEL", "gaarangoa/deeparg"),
                       help='Local path or Hugging Face repo id for a Transformers DeepARG model [Default: gaarangoa/deeparg].')
    reads.add_argument('--hf-device', default=None,
                       help='Torch device for --hf-model-path, for example cuda, cuda:0, or cpu [Default: auto]')
    reads.add_argument('--hf-model-subfolder', default=None,
                       help='Optional subfolder inside --hf-model-path, for example SS or LS.')
    reads.add_argument('--hf-chunk-size', default=10000, type=int,
                       help='Prediction chunk size for --hf-model-path [Default: 10000]')
    reads.set_defaults(func=predict)

    short_reads = subparsers.add_parser(
        "short_reads_pipeline",
        help="Run the modern paired-end short reads pipeline with DeepARG-SS",
    )
    short_reads.add_argument('--forward-pe-file', '--forward_pe_file', dest='forward_pe_file',
                             help='Forward mate from a paired-end library.')
    short_reads.add_argument('--reverse-pe-file', '--reverse_pe_file', dest='reverse_pe_file',
                             help='Reverse mate from a paired-end library.')
    short_reads.add_argument('--input-fasta', default=None,
                             help='Already-clean FASTA input. Skips Trimmomatic and vsearch.')
    short_reads.add_argument('--output-file', '--output_file', dest='output_file', required=True,
                             help='Output file prefix.')
    short_reads.add_argument('-d', '--deeparg-data-path', '--deeparg_data_path', dest='deeparg_data_path',
                             default=None,
                             help='Optional local DeepARG data path. If omitted, HF assets are used.')
    short_reads.add_argument('--hf-model-path', default=os.environ.get("DEEPARG_HF_MODEL", "gaarangoa/deeparg"),
                             help='Local path or Hugging Face repo id [Default: gaarangoa/deeparg].')
    short_reads.add_argument('--diamond-path', default=None,
                             help='Optional path to a native DIAMOND executable.')
    short_reads.add_argument('--deeparg-identity', '--deeparg_identity', dest='deeparg_identity',
                             default=80, type=float,
                             help='Minimum identity for ARG alignments [Default: 80].')
    short_reads.add_argument('--deeparg-probability', '--deeparg_probability', dest='deeparg_probability',
                             default=0.8, type=float,
                             help='Minimum probability for ARG-like reads [Default: 0.8].')
    short_reads.add_argument('--deeparg-evalue', '--deeparg_evalue', dest='deeparg_evalue',
                             default=1e-10, type=float,
                             help='Maximum e-value for ARG alignments [Default: 1e-10].')
    short_reads.add_argument('--gene-coverage', '--gene_coverage', dest='gene_coverage',
                             default=1, type=float,
                             help='Minimum gene coverage percentage used by 16S-normalized outputs [Default: 1].')
    short_reads.add_argument('--normalize-16s', action='store_true',
                             help='Run optional bowtie2/samtools 16S normalization against gg13.')
    short_reads.add_argument('--gg13-prefix', default=None,
                             help='Bowtie2 index prefix for Greengenes gg13, e.g. /path/gg13/dataset.')
    short_reads.add_argument('--skip-trimming', action='store_true',
                             help='Use paired FASTQ files directly as vsearch input.')
    short_reads.add_argument('--hf-device', default=None,
                             help='Torch device for model inference [Default: auto].')
    short_reads.add_argument('--hf-chunk-size', default=10000, type=int,
                             help='Prediction chunk size [Default: 10000].')
    short_reads.set_defaults(func=short_reads_pipeline)

    # Get all arguments
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    return args.func(args)


if __name__ == '__main__':
    main()
