import pathlib

from .v4_0 import s100 as v4_0
from .v4_0.s100 import *
import s100py.s102

def open(filename: (str, pathlib.Path), mode: str = "r") -> S1XXFile:
    # sends back the S100 version if it can't figure out which exact type to use
    file_object = S100File(filename, mode)
    # @todo use the new case statement when python 3.10 is required
    spec = file_object.root.product_specification
    if spec == b'INT.IHO.S-102.2.0':
        file_object = s100py.s102.v2_0.api.S102File(filename, mode)
    elif spec == b'INT.IHO.S-102.2.1':
        file_object = s100py.s102.v2_1.api.S102File(filename, mode)
    else:
        print("Warning: unrecognized s100 product specification, using generic S100File object")
    return file_object
