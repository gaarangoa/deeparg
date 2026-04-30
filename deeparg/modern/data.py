from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

from deeparg.modern.features import DeepARGFeatureEncoder


def read_gene_lengths(path):
    if path is None:
        return {}
    gene_lengths = {}
    with Path(path).open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) >= 2:
                gene_lengths[parts[0]] = float(parts[1])
    return gene_lengths


def read_alignment_tsv(
    path,
    identity=50.0,
    evalue=1e-5,
    coverage=0.8,
    features=None,
    gene_lengths=None,
    pipeline="reads",
    use_bitscore=True,
):
    alignments = {}
    best_hits = {}
    selected_features = set(features or [])
    gene_lengths = gene_lengths or {}
    measure_column = 11 if use_bitscore else 2

    with Path(path).open() as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) <= measure_column:
                continue
            if float(row[2]) < identity:
                continue
            if float(row[10]) > evalue:
                continue

            subject = row[1]
            if selected_features and subject not in selected_features:
                continue

            alignment_length = float(int(row[3]))
            if pipeline == "genes":
                gene_length = gene_lengths.get(subject)
                if gene_length is None:
                    continue
                if alignment_length / gene_length < coverage:
                    continue
            elif pipeline == "reads":
                if alignment_length <= coverage:
                    continue

            query = row[0]
            score = float(row[measure_column])
            alignments.setdefault(query, {})[subject] = score

            best_hit = best_hits.get(query)
            if best_hit is None or best_hit[1] < score:
                best_hits[query] = [subject, score, row]

    return alignments, best_hits


def labels_from_samples(samples):
    labels = []
    for sample in samples:
        parts = sample.split("|")
        if len(parts) < 3:
            raise ValueError(
                "Cannot extract class label from sample id {!r}; expected pipe-separated field 3.".format(sample)
            )
        labels.append(parts[2])
    return labels


def build_label_maps(labels):
    label_names = sorted(set(labels))
    label2id = {label: index for index, label in enumerate(label_names)}
    id2label = {index: label for label, index in label2id.items()}
    return label2id, id2label


def build_training_arrays(alignments, features=None, label2id=None, normalization="batch_minmax"):
    if features is None:
        features = sorted({feature for sample in alignments.values() for feature in sample.keys()})
    encoder = DeepARGFeatureEncoder(features, normalization=normalization)
    samples, inputs = encoder.alignments_to_numpy(alignments)
    labels = labels_from_samples(samples)
    if label2id is None:
        label2id, id2label = build_label_maps(labels)
    else:
        id2label = {index: label for label, index in label2id.items()}
    targets = np.array([label2id[label] for label in labels], dtype=np.int64)
    return samples, inputs, targets, label2id, id2label


class DeepARGDataset(Dataset):
    def __init__(self, input_features, labels=None):
        self.input_features = torch.as_tensor(input_features, dtype=torch.float32)
        self.labels = None if labels is None else torch.as_tensor(labels, dtype=torch.long)

    def __len__(self):
        return self.input_features.shape[0]

    def __getitem__(self, index):
        item = {"input_features": self.input_features[index]}
        if self.labels is not None:
            item["labels"] = self.labels[index]
        return item

