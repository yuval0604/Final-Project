from setuptools import setup, Extension
import numpy as np

module = Extension(
    'symnmf',
    sources=['symnmfmodule.c', 'symnmf.c'],
    include_dirs=[np.get_include()],
    extra_compile_args=['-std=c99','-Wall','-Wextra','-Werror','-pedantic-errors', "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION"],
)

setup(
    name='symnmf',
    version='1.0',
    description='SymNMF extension',
    ext_modules=[module]
)
