import operator
import os
import socket

import torch

from deeparg.modern.data import read_alignment_tsv, read_gene_lengths
from deeparg.modern.features import DeepARGFeatureEncoder
from deeparg.modern.modeling_deeparg import DeepARGForSequenceClassification


def force_ipv4_if_requested():
    if os.environ.get("HF_FORCE_IPV4") != "1":
        return
    original_getaddrinfo = socket.getaddrinfo

    def getaddrinfo_ipv4(*args, **kwargs):
        return [
            item for item in original_getaddrinfo(*args, **kwargs)
            if item[0] == socket.AF_INET
        ]

    socket.getaddrinfo = getaddrinfo_ipv4


def chunks(mapping, size=10000):
    keys = list(mapping.keys())
    for index in range(0, len(keys), size):
        selected = keys[index:index + size]
        yield {key: mapping[key] for key in selected}


def load_model_and_encoder(model_path, device=None, subfolder=None):
    force_ipv4_if_requested()
    model = DeepARGForSequenceClassification.from_pretrained(
        model_path,
        subfolder=subfolder,
        trust_remote_code=True,
    )
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)
    model.eval()
    encoder = DeepARGFeatureEncoder(
        features=model.config.features,
        normalization=model.config.feature_normalization,
    )
    return model, encoder, device


@torch.no_grad()
def predict_alignments(model, encoder, alignments, chunk_size=10000, device=None):
    if device is None:
        device = next(model.parameters()).device

    predictions = []
    for chunk in chunks(alignments, size=chunk_size):
        sample_ids, input_features = encoder.alignments_to_torch(chunk, device=device)
        probabilities = model.predict_proba(input_features).detach().cpu()
        for sample_id, row in zip(sample_ids, probabilities):
            order = torch.argsort(row)
            best = int(order[-1])
            best_prob = float(row[best])
            if 0.0 < best_prob < 0.9 and row.numel() > 1:
                second = int(order[-2])
                predictions.append([sample_id, model.config.id2label[second], float(row[second])])
            predictions.append([sample_id, model.config.id2label[best], best_prob])
    return predictions


def predict_file(
    alignment_tsv,
    output_prefix,
    model_path,
    identity,
    evalue,
    probability,
    coverage,
    pipeline,
    gene_lengths_path=None,
    chunk_size=10000,
    device=None,
    subfolder=None,
):
    model, encoder, device = load_model_and_encoder(model_path, device=device, subfolder=subfolder)
    gene_lengths = read_gene_lengths(gene_lengths_path)
    alignments, best_hits = read_alignment_tsv(
        alignment_tsv,
        identity=identity,
        evalue=evalue,
        coverage=coverage,
        features=encoder.features,
        gene_lengths=gene_lengths,
        pipeline=pipeline,
        use_bitscore=True,
    )
    predictions = predict_alignments(model, encoder, alignments, chunk_size=chunk_size, device=device)
    write_prediction_outputs(predictions, alignments, best_hits, output_prefix, probability)


def write_prediction_outputs(predictions, alignments, best_hits, output_prefix, probability):
    header = "#ARG\tquery-start\tquery-end\tread_id\tpredicted_ARG-class\tbest-hit\tprobability\tidentity\talignment-length\talignment-bitscore\talignment-evalue\tcounts\n"
    with open(output_prefix + ".ARG", "w") as confident, open(output_prefix + ".potential.ARG", "w") as potential:
        confident.write(header)
        potential.write(header)
        for sample_id, label, prob in predictions:
            sample_alignments = alignments[sample_id]
            matching_alignments = {
                key: value for key, value in sample_alignments.items() if "|" + label + "|" in key
            }
            if matching_alignments:
                best_hit = max(matching_alignments.items(), key=operator.itemgetter(1))[0]
                output_best_hit = best_hit
            else:
                best_hit = best_hits[sample_id][0]
                output_best_hit = "undefined"

            best_hit_row = best_hits[sample_id][2]
            row = [
                best_hit.split("|")[-1].upper(),
                best_hit_row[8],
                best_hit_row[9],
                sample_id,
                label,
                output_best_hit,
                str(prob),
                best_hit_row[2],
                best_hit_row[3],
                best_hit_row[-1],
                best_hit_row[-2],
                "1",
            ]
            handle = confident if prob >= probability and matching_alignments else potential
            handle.write("\t".join(row) + "\n")
