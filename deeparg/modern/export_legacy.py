import argparse
import json
import pickle
import shutil
from pathlib import Path

import torch

from deeparg.modern.configuration_deeparg import DeepARGConfig
from deeparg.modern.modeling_deeparg import DeepARGForSequenceClassification


def parse_args():
    parser = argparse.ArgumentParser(description="Export legacy DeepARG metadata/model PKLs as a Transformers model.")
    parser.add_argument("--metadata-pkl", required=True, help="Legacy metadata_SS.pkl or metadata_LS.pkl.")
    parser.add_argument("--model-pkl", required=True, help="Legacy Lasagne model_SS.pkl or model_LS.pkl.")
    parser.add_argument("--output-dir", required=True, help="Directory for save_pretrained output.")
    parser.add_argument("--database-dir", default=None, help="Optional database directory to copy into output/database.")
    parser.add_argument("--feature-lengths", default=None, help="Optional features.gene.length file to copy.")
    return parser.parse_args()


def load_pickle(path):
    with open(path, "rb") as handle:
        return pickle.load(handle, encoding="latin1")


def label_maps(metadata):
    raw_id2label = metadata["Y_rev"]
    id2label = {int(index): label for index, label in raw_id2label.items()}
    label2id = {label: index for index, label in id2label.items()}
    return label2id, id2label


def copy_legacy_params(model, params):
    if hasattr(params, "values"):
        flattened = []
        for layer_params in params.values():
            if len(layer_params) == 2:
                flattened.extend(layer_params)
        params = flattened

    linear_layers = [module for module in model.modules() if isinstance(module, torch.nn.Linear)]
    expected = len(linear_layers) * 2
    if len(params) != expected:
        raise ValueError("Expected {} Lasagne params, got {}.".format(expected, len(params)))

    with torch.no_grad():
        for layer, weight, bias in zip(linear_layers, params[0::2], params[1::2]):
            layer.weight.copy_(torch.as_tensor(weight, dtype=layer.weight.dtype).t())
            layer.bias.copy_(torch.as_tensor(bias, dtype=layer.bias.dtype))


def write_model_card(output_dir, metadata):
    readme = Path(output_dir) / "README.md"
    card = {
        "source": "legacy-lasagne",
        "input_nodes": metadata.get("input_nodes"),
        "output_nodes": metadata.get("output_nodes"),
        "num_features": len(metadata.get("features", [])),
    }
    readme.write_text(
        "\n".join(
            [
                "---",
                "library_name: transformers",
                "tags:",
                "- deeparg",
                "- antibiotic-resistance",
                "- pytorch",
                "- legacy-conversion",
                "---",
                "",
                "# DeepARG Legacy Conversion",
                "",
                "This model was converted from the original DeepARG Lasagne/Theano dense classifier.",
                "Inputs are normalized alignment bit-score feature vectors ordered by `config.features`.",
                "",
                "```json",
                json.dumps(card, indent=2, sort_keys=True),
                "```",
                "",
            ]
        )
    )


def json_default(value):
    if hasattr(value, "item"):
        return value.item()
    return str(value)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = load_pickle(args.metadata_pkl)
    params = load_pickle(args.model_pkl)
    label2id, id2label = label_maps(metadata)
    config = DeepARGConfig(
        input_size=int(metadata["input_nodes"]),
        num_labels=int(metadata["output_nodes"]),
        features=list(metadata["features"]),
        id2label=id2label,
        label2id=label2id,
        hidden_sizes=[2000, 1000, 500, 100],
        dropout=0.5,
        dropout_after_layers=3,
        feature_normalization="batch_minmax",
    )
    model = DeepARGForSequenceClassification(config)
    copy_legacy_params(model, params)

    DeepARGConfig.register_for_auto_class()
    DeepARGForSequenceClassification.register_for_auto_class("AutoModelForSequenceClassification")
    model.save_pretrained(output_dir, safe_serialization=True)
    config.save_pretrained(output_dir)
    (output_dir / "legacy_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True, default=json_default)
    )
    write_model_card(output_dir, metadata)

    if args.database_dir:
        destination = output_dir / "database"
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(args.database_dir, destination)
    if args.feature_lengths:
        shutil.copy2(args.feature_lengths, output_dir / "features.gene.length")


if __name__ == "__main__":
    main()
