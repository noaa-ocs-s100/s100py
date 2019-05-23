import setuptools
import os
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# Get the long description from the README file
with open(os.path.join(HERE, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='s100py',
    version='0.2.0',
    license='',
    author='Erin Nagel, Jason Greenlaw',
    author_email='erin.nagel@noaa.gov, jason.greenlaw@noaa.gov',
    description='Convert Oceanographic Operational Forecast System Datasets into IHO S-100/S-111 Format',
    long_description=long_description,
    url='',
    packages=setuptools.find_packages(),
    install_requires=['h5py', 'numpy', 'thyme(>=0.2.0)'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved',
        'Operating System :: OS Independent',
    ],
)
