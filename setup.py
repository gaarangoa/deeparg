from setuptools import setup, find_packages

setup(
    name='deeparg-modern',
    version='2.0.0',
    packages=find_packages(
        exclude=(".git", "data")
    ),
    include_package_data=True,
    package_data={},
    install_requires=[
        'accelerate',
        'biopython',
        'huggingface-hub',
        'numpy',
        'requests',
        'safetensors',
        'scikit-learn',
        'torch',
        'tqdm',
        'transformers',
        'wget',
    ],
    python_requires=">=3.10",
    entry_points={
        'console_scripts': [
            'deeparg=deeparg.entry:main',
            'deeparg-short-reads=deeparg.short_reads_pipeline.pipeline:main',
            'deeparg-train-transformers=deeparg.modern.train_transformers:main',
            'deeparg-export-legacy-transformers=deeparg.modern.export_legacy:main',
            'deeparg-upload-hf=deeparg.modern.upload_hf:main',
        ],
    },
)
