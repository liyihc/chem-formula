from Cython.Build import cythonize
import numpy as np

cythonize("chem_formula/*.pyx", include_path=[np.get_include()], force=True)
