from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
Plumedplugin_header_dir = '@PlumedPLUGIN_HEADER_DIR@'
Plumedplugin_library_dir = '@PlumedPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = []
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_PLUMEDplugin',
                      sources=['PlumedPluginWrapper.cpp'],
                      libraries=['OpenMM', 'PlumedPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), Plumedplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), Plumedplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='Plumedplugin',
      version='1.0',
      py_modules=['Plumedplugin'],
      ext_modules=[extension],
     )
