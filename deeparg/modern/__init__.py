__all__ = [
    "DeepARGConfig",
    "DeepARGFeatureEncoder",
    "DeepARGForSequenceClassification",
]


def __getattr__(name):
    if name == "DeepARGConfig":
        from deeparg.modern.configuration_deeparg import DeepARGConfig

        return DeepARGConfig
    if name == "DeepARGFeatureEncoder":
        from deeparg.modern.features import DeepARGFeatureEncoder

        return DeepARGFeatureEncoder
    if name == "DeepARGForSequenceClassification":
        from deeparg.modern.modeling_deeparg import DeepARGForSequenceClassification

        return DeepARGForSequenceClassification
    raise AttributeError(name)
