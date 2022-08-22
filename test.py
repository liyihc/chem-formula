from Cython.Build import cythonize
from distutils.extension import Extension

ext:Extension
for ext in cythonize("chem_formula/*.pyx"):
    print(ext.name, ext.sources)