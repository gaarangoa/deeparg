import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.model_selection import train_test_split
from transformers import Trainer, TrainingArguments

from deeparg.modern.configuration_deeparg import DeepARGConfig
from deeparg.modern.data import (
    DeepARGDataset,
    build_training_arrays,
    read_alignment_tsv,
    read_gene_lengths,
)
from deeparg.modern.modeling_deeparg import DeepARGForSequenceClassification


def parse_args():
    parser = argparse.ArgumentParser(description="Train a DeepARG Transformers model from DIAMOND/BLAST TSV.")
    parser.add_argument("--train-alignments", required=True, help="Training alignment TSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for save_pretrained output.")
    parser.add_argument("--pipeline", choices=["reads", "genes"], default="reads")
    parser.add_argument("--gene-lengths", default=None, help="features.gene.length path, required for genes.")
    parser.add_argument("--identity", type=float, default=30.0)
    parser.add_argument("--evalue", type=float, default=1.0)
    parser.add_argument("--coverage", type=float, default=25.0)
    parser.add_argument("--epochs", type=float, default=100.0)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--learning-rate", type=float, default=1e-2)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--validation-split", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--normalization", choices=["batch_minmax", "none"], default="batch_minmax")
    parser.add_argument("--hidden-sizes", default="2000,1000,500,100")
    parser.add_argument("--dropout", type=float, default=0.5)
    parser.add_argument("--use-cpu", action="store_true")
    return parser.parse_args()


def write_model_card(output_dir, metadata):
    readme = Path(output_dir) / "README.md"
    readme.write_text(
        "\n".join(
            [
                "---",
                "library_name: transformers",
                "tags:",
                "- deeparg",
                "- antibiotic-resistance",
                "- pytorch",
                "---",
                "",
                "# DeepARG Transformers Model",
                "",
                "This repository contains a PyTorch/Transformers wrapper for the DeepARG dense classifier.",
                "Inputs are normalized alignment bit-score feature vectors ordered by `config.features`.",
                "",
                "## Training metadata",
                "",
                "```json",
                json.dumps(metadata, indent=2, sort_keys=True),
                "```",
                "",
            ]
        )
    )


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_lengths = read_gene_lengths(args.gene_lengths)
    alignments, _ = read_alignment_tsv(
        args.train_alignments,
        identity=args.identity,
        evalue=args.evalue,
        coverage=args.coverage,
        gene_lengths=gene_lengths,
        pipeline=args.pipeline,
        use_bitscore=True,
    )
    samples, inputs, targets, label2id, id2label = build_training_arrays(
        alignments,
        normalization=args.normalization,
    )
    hidden_sizes = [int(value) for value in args.hidden_sizes.split(",") if value.strip()]

    if args.validation_split > 0:
        train_idx, eval_idx = train_test_split(
            np.arange(len(targets)),
            test_size=args.validation_split,
            random_state=args.seed,
            stratify=targets if len(set(targets.tolist())) > 1 else None,
        )
        train_dataset = DeepARGDataset(inputs[train_idx], targets[train_idx])
        eval_dataset = DeepARGDataset(inputs[eval_idx], targets[eval_idx])
        eval_strategy = "epoch"
    else:
        train_dataset = DeepARGDataset(inputs, targets)
        eval_dataset = None
        eval_strategy = "no"

    config = DeepARGConfig(
        input_size=inputs.shape[1],
        num_labels=len(label2id),
        features=sorted({feature for sample in alignments.values() for feature in sample.keys()}),
        id2label=id2label,
        label2id=label2id,
        hidden_sizes=hidden_sizes,
        dropout=args.dropout,
        feature_normalization=args.normalization,
    )
    model = DeepARGForSequenceClassification(config)

    training_args = TrainingArguments(
        output_dir=str(output_dir / "checkpoints"),
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        num_train_epochs=args.epochs,
        learning_rate=args.learning_rate,
        weight_decay=args.weight_decay,
        eval_strategy=eval_strategy,
        save_strategy="no",
        logging_strategy="epoch",
        report_to="none",
        remove_unused_columns=False,
        seed=args.seed,
        use_cpu=args.use_cpu,
    )
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=eval_dataset,
    )
    trainer.train()

    DeepARGConfig.register_for_auto_class()
    DeepARGForSequenceClassification.register_for_auto_class("AutoModelForSequenceClassification")
    model.save_pretrained(output_dir, safe_serialization=True)
    config.save_pretrained(output_dir)

    metadata = {
        "train_alignments": str(args.train_alignments),
        "num_samples": len(samples),
        "num_features": inputs.shape[1],
        "num_labels": len(label2id),
        "pipeline": args.pipeline,
        "identity": args.identity,
        "evalue": args.evalue,
        "coverage": args.coverage,
        "normalization": args.normalization,
    }
    (output_dir / "training_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True))
    write_model_card(output_dir, metadata)


if __name__ == "__main__":
    main()

