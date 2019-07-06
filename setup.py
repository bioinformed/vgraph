# Copyright 2015 Kevin B Jacobs
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License.  You may obtain
# a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations
# under the License.

"""Setup script for vgraph."""

import sys

from setuptools import setup, find_packages
from Cython.Distutils import build_ext
from distutils.extension import Extension

if sys.version_info < (3, 4):
    sys.exit('Sorry, Python 3.4 or newer is required to install and run vgraph')

install_requires = ['pysam']
setup_requires   = ['pysam', 'setuptools_scm==1.15.0', 'setuptools_scm_git_archive==1.0']
tests_require    = ['pytest-runner', 'pytest', 'coverage']


ext_modules = [
    Extension('vgraph.norm',      ['vgraph/norm.pyx']),
    Extension('vgraph.intervals', ['vgraph/intervals.pyx']),
]


classifiers = """
Development Status :: 2 - Alpha
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python :: 3.7
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""


if __name__ == '__main__':
    setup(
        name='vgraph',
        description='Graph-based variant normalization and comparison tools',
        url='https://github.com/bioinformed/vgraph',
        author='Kevin Jacobs',
        maintainer='Kevin Jacobs',
        author_email='jacobs@bioinformed.com',
        maintainer_email='jacobs@bioinformed.com',
        license='APACHE-2.0',
        classifiers=classifiers.split('\n'),
        use_scm_version=True,
        zip_safe=False,
        tests_require=tests_require,
        packages=find_packages(),
        install_requires=install_requires,
        setup_requires=setup_requires,
        cmdclass={'build_ext': build_ext},
        ext_modules=ext_modules,
        entry_points={'console_scripts': ['vgraph=vgraph.vgraph:main']},
    )
