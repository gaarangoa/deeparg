import torch
from torch import nn
from transformers import PreTrainedModel
from transformers.modeling_outputs import SequenceClassifierOutput

from deeparg.modern.configuration_deeparg import DeepARGConfig


class DeepARGForSequenceClassification(PreTrainedModel):
    config_class = DeepARGConfig
    base_model_prefix = "deeparg"
    main_input_name = "input_features"

    def __init__(self, config):
        super().__init__(config)
        if config.input_size <= 0:
            raise ValueError("DeepARGConfig.input_size must be greater than zero.")
        if config.num_labels <= 0:
            raise ValueError("DeepARGConfig.num_labels must be greater than zero.")

        layers = []
        previous_size = config.input_size
        for layer_index, hidden_size in enumerate(config.hidden_sizes):
            layers.append(nn.Linear(previous_size, hidden_size))
            layers.append(nn.ReLU())
            if layer_index < config.dropout_after_layers:
                layers.append(nn.Dropout(config.dropout))
            previous_size = hidden_size

        self.deeparg = nn.Sequential(*layers)
        self.classifier = nn.Linear(previous_size, config.num_labels)
        self.post_init()

    def forward(
        self,
        input_features=None,
        labels=None,
        output_attentions=None,
        output_hidden_states=None,
        return_dict=None,
        **kwargs,
    ):
        if input_features is None:
            input_features = kwargs.get("inputs_embeds")
        if input_features is None:
            raise ValueError("DeepARG requires input_features.")

        input_features = input_features.to(dtype=self.dtype)
        hidden_states = self.deeparg(input_features)
        logits = self.classifier(hidden_states)

        loss = None
        if labels is not None:
            loss = nn.CrossEntropyLoss()(logits.view(-1, self.config.num_labels), labels.view(-1))

        if return_dict is False:
            output = (logits,)
            if output_hidden_states:
                output = output + (hidden_states,)
            return ((loss,) + output) if loss is not None else output

        hidden_states_output = (hidden_states,) if output_hidden_states else None
        return SequenceClassifierOutput(
            loss=loss,
            logits=logits,
            hidden_states=hidden_states_output,
            attentions=None,
        )

    @torch.no_grad()
    def predict_proba(self, input_features):
        was_training = self.training
        self.eval()
        outputs = self(input_features=input_features)
        probabilities = torch.softmax(outputs.logits, dim=-1)
        if was_training:
            self.train()
        return probabilities
