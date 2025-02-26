import shutil

from .v2_0 import api as v2_0
from .v2_1 import api as v2_1
from .v2_2 import api as v2_2
from .v3_0 import api as v3_0
from .v3_0.api import *


def open(name, *args, version=EDITION, **kwargs):
    obj = None
    if version == 3.0:
        obj = v3_0.S102File(name, *args, **kwargs)
    elif version == 2.2:
        obj = v2_2.S102File(name, *args, **kwargs)
    elif version == 2.1:
        obj = v2_1.S102File(name, *args, **kwargs)
    elif version == 2.0:
        obj = v2_0.S102File(name, *args, **kwargs)
    elif version < 2.0:
        raise NotImplementedError("Version 1.x of S102 is not supported")
    else:
        raise ValueError(f"Version {version} is not supported in this version of the S100py module")
    return obj


# update a file in place
def update(name, *args, version=EDITION, **kwargs):
    raise NotImplementedError("not done yet")


def copy_update(src_name, dest_name, *args, version=EDITION, **kwargs):
    shutil.copy(src_name, dest_name)
    return update(dest_name, *args, version=version, **kwargs)
