from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.dist import Distribution

setup(
    name="chem_formula",
    version="0.1.0",
    description="test formula",
    packages=find_packages(),
    # package_dir={"": "chem_formula"},
    package_data={"chem_formula": ["*.csv", "*.pxd", "*.pyx", "*.cpp"]},
    python_requires=">=3.6, <4",
    install_requires=[
        # "Cython>=0.29.1, <0.30",
        "numpy>=1.17",
        "pyteomics>=4.0"
    ]
)
