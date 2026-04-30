from types import SimpleNamespace


DEFAULT_MODEL_REPO = "gaarangoa/deeparg"


def predict(
    input_file,
    output_file,
    model="LS",
    input_type="nucl",
    model_repo=DEFAULT_MODEL_REPO,
    diamond_path=None,
    data_path=None,
    min_prob=0.8,
    identity=50,
    evalue=1e-10,
    overlap=0.8,
    num_alignments_per_entry=1000,
    model_version="v2",
    device=None,
    chunk_size=10000,
):
    from deeparg.entry import predict as cli_predict

    model = model.upper()
    args = SimpleNamespace(
        model=model,
        input_file=input_file,
        output_file=output_file,
        data_path=data_path,
        type=input_type,
        min_prob=min_prob,
        arg_alignment_identity=identity,
        arg_alignment_evalue=evalue,
        arg_alignment_overlap=overlap,
        arg_num_alignments_per_entry=num_alignments_per_entry,
        model_version=model_version,
        diamond_path=diamond_path,
        hf_model_path=model_repo,
        hf_model_subfolder=model,
        hf_device=device,
        hf_chunk_size=chunk_size,
    )
    cli_predict(args)
    return {
        "arg": output_file + ".mapping.ARG",
        "potential_arg": output_file + ".mapping.potential.ARG",
        "alignment_tsv": output_file + ".align.daa.tsv",
    }

