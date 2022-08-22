from setuptools import setup, find_packages, Extension

# pyinstaller
# 可参考https://github.com/pyinstaller/pyinstaller-hooks-contrib
# 其实就是entrypoint，在哪里定义都行
setup(
    name="chem_formula",
    version="0.1.0",
    description="test formula",
    packages=find_packages(),
    package_data={"chem_formula": ["*.csv", "*.pxd", "*.pyx", "*.cpp"]},
    python_requires=">=3.6, <4",
    install_requires=[
        # "Cython>=0.29.1, <0.30",
        "numpy>=1.17",
        "pyteomics>=4.0"
    ]
)
