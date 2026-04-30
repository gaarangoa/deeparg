from transformers import PretrainedConfig


class DeepARGConfig(PretrainedConfig):
    model_type = "deeparg"

    def __init__(
        self,
        input_size=0,
        num_labels=0,
        features=None,
        id2label=None,
        label2id=None,
        hidden_sizes=None,
        dropout=0.5,
        dropout_after_layers=None,
        feature_normalization="batch_minmax",
        problem_type="single_label_classification",
        **kwargs,
    ):
        features = list(features or [])
        hidden_sizes = list(hidden_sizes or [2000, 1000, 500, 100])
        if dropout_after_layers is None:
            dropout_after_layers = max(0, len(hidden_sizes) - 1)

        if input_size == 0 and features:
            input_size = len(features)
        if num_labels == 0 and id2label:
            num_labels = len(id2label)

        super().__init__(
            num_labels=num_labels,
            id2label=id2label,
            label2id=label2id,
            problem_type=problem_type,
            **kwargs,
        )
        self.input_size = int(input_size)
        self.features = features
        self.hidden_sizes = hidden_sizes
        self.dropout = float(dropout)
        self.dropout_after_layers = int(dropout_after_layers)
        self.feature_normalization = feature_normalization
