from .v2_2.utils import *
"""
Filenaming supposed to be 102PPPPØØØØØØØØØØØØ.H5
102 - the first 3 characters identify the dataset as an S-102 dataset (mandatory).

PPPP - the fourth to seventh characters identify the producer code of the issuing agency (mandatory for S-102).   
Where the producer code is derived from a 2- or 3-character format (for instance when converting S-57 ENCs), 
the missing characters of the producer code must be populated with zeros (“00” or “0” respectively) 
for the sixth and seventh characters of the dataset file name, as required.

ØØØØØØØØØØØØ - the eighth to the maximum nineteenth characters are optional and may be used in any way by the producer to provide the unique file name.   
The following characters are allowed in the dataset name: A to Z, 0 to 9 and the special character _ (underscore).
H5 - denotes and HDF5 file.
"""