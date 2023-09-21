import numpy as np
import os
import platform
from setuptools import setup, Extension, find_packages
from os.path import join as pjoin


# IFT_DIR = "../" # changed to make installation not directory dependent
IFT_DIR = os.environ["NEWIFT_DIR"]
CUDA_DIR1 = "/usr/local/cuda"
CUDA_DIR2 = "/opt/cuda"

libraries = ['ift', 'svm', 'lapack', 'blas', 'm']

# adding the CUDA libraries
if os.path.exists(CUDA_DIR1) or os.path.exists(CUDA_DIR2):
    libraries.extend(['cublas', 'cudart'])

include_dirs = [dirpath for dirpath, _, _ in os.walk(pjoin(IFT_DIR, "include"))] +\
               [pjoin(IFT_DIR, 'externals/libsvm/include'),
                pjoin(IFT_DIR, 'externals/libnifti/include'),
                pjoin(IFT_DIR, 'externals/libjpeg/include'),
                pjoin(IFT_DIR, 'externals/libpng/include'),
                pjoin(IFT_DIR, 'externals/zlib/include'),
                pjoin(IFT_DIR, 'externals/tsne/include'),
                CUDA_DIR1,
                CUDA_DIR2,
                '/usr/local/opt/openblas/include',  # added this directory for MAC OS
                np.get_include(),
                ]

library_dirs = [pjoin(IFT_DIR, 'lib'),
                pjoin(IFT_DIR, 'externals/libsvm/lib'),
                pjoin(CUDA_DIR1, 'lib64'),
                pjoin(CUDA_DIR2, 'lib64'),
                ]

extra_compile_args = ['-w', '-fPIC', '-O3', '-pthread', '-std=gnu11', '-pedantic']

extra_link_args = []

# MAC OS
if platform.system() == "Darwin":
    extra_compile_args.append("-Xpreprocessor -fopenmp")
    extra_link_args.append("-lomp")
else:
    libraries.append("stdc++")
    extra_compile_args.append("-fopenmp")
    extra_link_args.append("-lgomp")

extensions = [Extension('_pyift', ['pyift/pyift.i'],
                        swig_opts=['-py3'],
                        include_dirs=include_dirs,
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_compile_args=extra_compile_args,
                        extra_link_args=extra_link_args,
                        )]

setup(name='pyift',
      version='0.1',
      author='LIDS',
      maintainer='Jordao Bragantini',
      maintainer_email='jordao.bragantini@gmail.com',
      description="""Python  wrapper of libift""",
      ext_modules=extensions,
      packages=find_packages(),
      package_data={'': ['_pyift.so']},
      )

