from .v1_0 import api as v1_0
from .v1_1 import api as v1_1
from .v2_0 import api as v2_0
from .v2_0.api import *


def open(name, *args, version=EDITION, **kwargs):
    obj = None
    if version == 2.0:
        obj = v2_0.S104File(name, *args, **kwargs)
    if version == 1.1:
        obj = v1_1.S104File(name, *args, **kwargs)
    if version == 1.0:
        obj = v1_0.S104File(name, *args, **kwargs)
    elif version < 1.0:
        raise NotImplementedError("Version x.x of S104 is not supported")
    else:
        raise ValueError(f"Version {version} is not supported in this version of the S100py module")
    return obj
