import inspect
import mfdnres
from setuptools import setup

setup(
    name="mfdnres",
    version="0.0.1",
    author="Mark A. Caprio, University of Notre Dame",
    description=("A scripting library for universal results postprocessing"),
    license="MIT",
    packages=['mfdnres'],
    long_description=inspect.getdoc(mfdnres),
    classifiers=[],
)
