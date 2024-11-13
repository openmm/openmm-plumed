from distutils.core import setup
from distutils.extension import Extension
import os
import platform
import numpy

openmm_dir = '@OPENMM_DIR@'
openmmplumed_header_dir = '@OPENMMPLUMED_HEADER_DIR@'
openmmplumed_library_dir = '@OPENMMPLUMED_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args=['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_openmmplumed',
                      sources=['PlumedPluginWrapper.cpp'],
                      libraries=['OpenMM', 'OpenMMPlumed'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), openmmplumed_header_dir, numpy.get_include()],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), openmmplumed_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='OpenMMPlumed',
      version='2.1',
      py_modules=['openmmplumed'],
      ext_modules=[extension],
      install_requires=['openmm', 'numpy']
     )
