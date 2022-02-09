import inspect
from setuptools import setup

setup(
    name="mfdnres",
    version="0.0.1",
    author="Mark A. Caprio, University of Notre Dame",
    description=("A scripting library for universal results postprocessing"),
    license="MIT",
    packages=['mfdnres'],
    install_requires=[
        "numpy>=1.0.0",
    ],
    extras_require={
        'analysis': ['pandas>=1.2.0', 'matplotlib>=3.4.1', 'more_itertools']
    },
    classifiers=[],
)
