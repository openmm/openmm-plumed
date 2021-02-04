[![GH Actions Status](https://github.com/openmm/openmm-plumed/workflows/CI/badge.svg)](https://github.com/openmm/openmm-plumed/actions?query=branch%3Amaster+workflow%3ACI)
[![Conda](https://img.shields.io/conda/v/conda-forge/openmm-plumed.svg)](https://anaconda.org/conda-forge/openmm-plumed)
[![Anaconda Cloud Badge](https://anaconda.org/conda-forge/openmm-plumed/badges/downloads.svg)](https://anaconda.org/conda-forge/openmm-plumed)

OpenMM PLUMED Plugin
=====================

This project provides a connection between [OpenMM](http://openmm.org) and [PLUMED](http://www.plumed.org).
It allows you to bias or analyze an OpenMM simulation based on collective variables.

This plugin requires PLUMED version 2.3b or greater.

Installing The Plugin
=====================

We provide [conda](https://docs.conda.io/) packages for Linux and MacOS via [`conda-forge`](https://conda-forge.org/), which can be installed from the [conda-forge channel](https://anaconda.org/conda-forge/openmm-plumed):

```bash
conda install -c conda-forge openmm-plumed
```

If you don't have `conda` available, we recommend installing [Miniconda for Python 3](https://docs.conda.io/en/latest/miniconda.html) to provide the `conda` package manager.  

Building The Plugin
===================

This project uses [CMake](http://www.cmake.org) for its build system.  To build it, follow these
steps:

1. Create a directory in which to build the plugin.

2. Run the CMake GUI or ccmake, specifying your new directory as the build directory and the top
level directory of this project as the source directory.

3. Press "Configure".

4. Set OPENMM_DIR to point to the directory where OpenMM is installed.  This is needed to locate
the OpenMM header files and libraries.

5. Set PLUMED_INCLUDE_DIR and PLUMED_LIBRARY_DIR to point to the directories where the PLUMED header
files and libraries are installed.

6. Set CMAKE_INSTALL_PREFIX to the directory where the plugin should be installed.  Usually,
this will be the same as OPENMM_DIR, so the plugin will be added to your OpenMM installation.

7. If you plan to build the OpenCL platform, make sure that OPENCL_INCLUDE_DIR and
OPENCL_LIBRARY are set correctly, and that PLUMED_BUILD_OPENCL_LIB is selected.

8. If you plan to build the CUDA platform, make sure that CUDA_TOOLKIT_ROOT_DIR is set correctly
and that PLUMED_BUILD_CUDA_LIB is selected.

9. Press "Configure" again if necessary, then press "Generate".

10. Use the build system you selected to build and install the plugin.  For example, if you
selected Unix Makefiles, type `make install` to install the plugin, and `make PythonInstall` to
install the Python wrapper.

Using The Plugin
================

Simply create a `PlumedForce` object, passing the PLUMED control script as an argument to the
constructor, then add it to your `System`.  For example,

```Python
script = """
d: DISTANCE ATOMS=1,10
METAD ARG=d SIGMA=0.2 HEIGHT=0.3 PACE=500"""
system.addForce(PlumedForce(script))
```

Be aware that PLUMED numbers atoms starting from 1, whereas OpenMM numbers them starting from 0.
The example above performs metadynamics based on the distance between atoms 0 and 9 (in OpenMM's
numbering system).


License
=======

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2016 Stanford University and the Authors.

Authors: Peter Eastman

Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
