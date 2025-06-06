from .v1_0 import api as v1_0
from .v1_2 import api as v1_2
from .v2_0 import api as v2_0
from .v2_0.api import *


def open(name, *args, version=EDITION, **kwargs):
    obj = None
    if version == 2.0:
        obj = v2_0.S111File(name, *args, **kwargs)
    elif version == 1.2:
        obj = v1_2.S111File(name, *args, **kwargs)
    elif version == 1.0:
        obj = v1_0.S111File(name, *args, **kwargs)
    elif version == 1.1:
        raise NotImplementedError("Version 1.1 of S111 is not supported")
    elif version < 1.0:
        raise NotImplementedError("Version 1.x of S111 is not supported")
    else:
        raise ValueError(f"Version {version} is not supported in this version of the S100py module")
    return obj
