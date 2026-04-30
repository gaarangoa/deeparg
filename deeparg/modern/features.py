import json
from pathlib import Path

import numpy as np
import torch


class DeepARGFeatureEncoder:
    def __init__(self, features, normalization="batch_minmax"):
        self.features = list(features)
        self.feature_to_index = {feature: index for index, feature in enumerate(self.features)}
        self.normalization = normalization

    @classmethod
    def from_pretrained(cls, path_or_model):
        config_path = Path(path_or_model) / "config.json"
        if not config_path.exists():
            raise FileNotFoundError("Expected config.json in {}".format(path_or_model))
        with config_path.open() as handle:
            config = json.load(handle)
        return cls(
            features=config.get("features", []),
            normalization=config.get("feature_normalization", "batch_minmax"),
        )

    def alignments_to_numpy(self, alignments, dtype=np.float32):
        sample_ids = list(alignments.keys())
        matrix = np.zeros((len(sample_ids), len(self.features)), dtype=dtype)
        for row_index, sample_id in enumerate(sample_ids):
            for feature, value in alignments[sample_id].items():
                column_index = self.feature_to_index.get(feature)
                if column_index is not None:
                    matrix[row_index, column_index] = float(value)
        return sample_ids, self.normalize(matrix)

    def alignments_to_torch(self, alignments, device=None):
        sample_ids, matrix = self.alignments_to_numpy(alignments)
        tensor = torch.from_numpy(matrix)
        if device is not None:
            tensor = tensor.to(device)
        return sample_ids, tensor

    def normalize(self, matrix):
        if self.normalization in (None, "none"):
            return matrix
        if self.normalization != "batch_minmax":
            raise ValueError("Unsupported feature normalization: {}".format(self.normalization))
        if matrix.size == 0:
            return matrix

        minimum = matrix.min(axis=0)
        maximum = matrix.max(axis=0)
        scale = maximum - minimum
        normalized = matrix.copy()
        non_constant = scale != 0
        normalized[:, non_constant] = (
            normalized[:, non_constant] - minimum[non_constant]
        ) / scale[non_constant]
        normalized[:, ~non_constant] = 0
        return normalized.astype(np.float32, copy=False)

    def to_dict(self):
        return {
            "features": self.features,
            "normalization": self.normalization,
        }

