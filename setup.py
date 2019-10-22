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
    version='0.4.0',
    license='BSD-2-Clause',
    author='Erin Nagel, Jason Greenlaw',
    author_email='erin.nagel@noaa.gov, jason.greenlaw@noaa.gov',
    description='This python package provides utilities for encoding hydrographic datasets in the International Hydrographic Organization (IHO) S-100 format',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/noaa-ocs-s100/s100py',
    packages=setuptools.find_packages(),
    install_requires=['thyme(>=0.4.0)', 'numpy', 'h5py'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
)
