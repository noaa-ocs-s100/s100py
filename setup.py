import setuptools
import os
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# Get the long description from the README file
with open(os.path.join(HERE, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='s100py',
    license='CC0 1.0 Universal',
    author='Barry Gallagher, Erin Nagel',
    author_email='barry.gallagher@noaa.gov, erin.nagel@noaa.gov',
    description='This python package provides api and utilities for encoding hydrographic datasets in the International Hydrographic Organization (IHO) S-100 format',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/noaa-ocs-s100/s100py',
    packages=setuptools.find_packages(),
    use_scm_version=True,
    setup_requires=['numpy', 'setuptools_scm'],
    install_requires=['numpy', 'h5py'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
        'Operating System :: OS Independent',
    ],
)
