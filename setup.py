from setuptools import setup, find_packages
from Cython.Distutils import build_ext
from distutils.extension import Extension

install_requires = ['Cython>=0.22', 'nose', 'pysam>=0.8.3']
setup_requires   = ['pysam>=0.8.3']
tests_require    = ['coverage']


ext_modules = [Extension('vgraph.norm', ['vgraph/norm.pyx'])]


classifiers = """
Development Status :: 2 - Alpha
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""


if __name__ == '__main__':
    setup(
        name = 'vgraph',
        version = '0.2',
        description = 'Graph-based variant normalization and comparison tools',
        url = 'https://github.com/bioinformed/vgraph',
        author = 'Kevin Jacobs',
        maintainer = 'Kevin Jacobs',
        author_email = 'jacobs@bioinformed.com',
        maintainer_email = 'jacobs@bioinformed.com',
        license = 'APACHE-2.0',
        classifiers = classifiers,
        zip_safe = False,
        test_suite = 'nose.collector',
        tests_require = tests_require,
        packages = find_packages(),
        install_requires = install_requires,
        setup_requires = setup_requires,
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules,
        scripts=['bin/vgraph'],
    )
