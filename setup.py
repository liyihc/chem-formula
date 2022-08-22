from pathlib import Path
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.dist import Distribution

setup(
    name="chem_formula",
    version="0.1.0",
    description="test formula",
    packages=find_packages(),
    python_requires=">=3.6, <4",
    install_requires=[
        # "Cython>=0.29.1, <0.30",
        "pyteomics>=4.0"
    ],
    # distclass=ChemFormulaDist,
    # https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
)
