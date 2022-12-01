import inspect
from setuptools import setup, find_packages

setup(
    name="mfdnres",
    version="0.1.0",
    author="Mark A. Caprio, University of Notre Dame",
    description=("A scripting library for universal results postprocessing"),
    license="MIT",
    packages=find_packages(include='mfdnres*'),
    python_requires='>=3.8',
    install_requires=[
        "numpy>=1.0.0",
        "more-itertools>=2.4",
    ],
    extras_require={
        'analysis': ['pandas>=1.2.0', 'matplotlib>=3.4.1', 'more_itertools']
    },
    classifiers=[],
)
