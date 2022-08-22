import pyteomics.mass  # just for pyinstaller package
import pyteomics  # just for pyinstaller package
try:
    from . import _element 
except:
    print("Building")
    from setuptools import setup, Extension
    import numpy as np
    from pathlib import Path

    extensions = []
    root = Path(__file__).parent
    for file in root.glob("*.cpp"):
        name = f"chem_formula.{file.stem}"
        extensions.append(Extension(name, [str(file)]))
    tmp_dir = root / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    setup(
        script_args=["build_ext"],
        include_dirs=np.get_include(),
        options={
            "build_ext": {
                "build_lib": str(root.parent),
                "build_temp": str(tmp_dir)
            }},
        ext_modules=extensions
    )
from . import _element
from . import _formula
from . import _restrictedCalc
from ._formula import Formula
from ._restrictedCalc import Calculator as RestrictedCalc
from ._forceCalc import Calculator as ForceCalc
