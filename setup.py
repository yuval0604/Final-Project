from setuptools import setup, Extension
import numpy as np

module = Extension(
    'symnmf',
    sources=['symnmfmodule.c', 'symnmf.c'],
    include_dirs=[np.get_include()],
    extra_compile_args=['-ansi','-Wall','-Wextra','-Werror','-pedantic-errors'],
)

setup(
    name='symnmf',
    version='1.0',
    description='SymNMF extension',
    ext_modules=[module]
)
