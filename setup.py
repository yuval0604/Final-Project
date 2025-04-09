from setuptools import setup, Extension

module = Extension(
    'symnmf',
    sources=['symnmfmodule.c', 'symnmf.c'],
    extra_compile_args=['-ansi','-Wall','-Wextra','-Werror','-pedantic-errors'],
)

setup(
    name='symnmf',
    version='1.0',
    description='SymNMF extension',
    ext_modules=[module]
)
