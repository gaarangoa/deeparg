from setuptools import setup, find_packages

setup(
    name='DeepARG',
    version='1.0.2',
    packages=find_packages(
        exclude=(".git", "data")
    ),
    include_package_data=True,
    package_data={},
    install_requires=[
        'numpy==1.10.4',
        'scipy==0.16.1',
        'joblib==0.13.0',
        'BioPython==1.73',
        'ete3',
        'tqdm',
        'nolearn==0.6',
        'lasagne==0.1',
        'scikit-learn==0.19.1',
        'theano==0.8.2',
        'requests==2.18.4',
        'wget'
    ],
    python_requires=">=2.7.15,!=3.*",
    entry_points='''
        [console_scripts]
        deeparg=deeparg.entry:main
    ''',
)
