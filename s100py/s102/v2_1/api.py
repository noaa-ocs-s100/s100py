try:
    from ... import s1xx
except:  # fake out sphinx and autodoc which are loading the module directly and losing the namespace
    __package__ = "s100py.s102"

from ..v2_0 import api as v2_0
# Anything not overridden in this module will use whatever was available in the previous version
from ..v2_0.api import *

EDITION = 2.1

CHANGELOG = """
Removed TrackingList  --  4.2.1.1.8 TrackingListCoverage
Removed min/max display scale -- 4.2.1.1.1.2 and 4.2.1.1.1.5 BathymetryCoverage semantics
Added flip_z parameters in utils since z orientation is going from positive up to positive down
 
"""

# removed the tracking list mixin
class FeatureCodes(v2_0.FeatureCodesBase):
    pass


# removed min/max display scale mixin
class BathymetryCoverage(v2_0.BathymetryCoverageBase):
    pass


class BathymetryGroupList(v2_0.BathymetryGroupList):
    @property
    def metadata_type(self) -> type:
        return BathymetryCoverage


class BathymetryFeatureInstance(v2_0.BathymetryFeatureInstance):
    @property
    def __bathymetry_group_type__(self):
        return BathymetryGroupList


class BathymetryCoveragesList(v2_0.BathymetryCoveragesList):
    @property
    def metadata_type(self) -> Type[BathymetryFeatureInstance]:
        return BathymetryFeatureInstance


class BathymetryContainer(v2_0.BathymetryContainer):
    """ This is the BathymetryCoverage right off the root of the HDF5 which has possible attributes from S100 spec table 10c-10
    This will hold child groups named BathymetryCoverage.NN
    """
    #: attribute name will be automatically determined based on the containing list's index
    __bathymetry_coverage_hdf_name__ = BATHY_COVERAGE + r"[\._]\d+"

    @property
    def __version__(self) -> int:
        return 1

    @property
    def __bathymetry_coverage_type__(self):
        return BathymetryCoveragesList


class S102Root(v2_0.S102RootBase):
    @property
    def __feature_information_type__(self):
        return FeatureCodes

    @property
    def __bathymetry_coverage_type__(self):
        return BathymetryContainer


class S102File(v2_0.S102File):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-102.2.1')
    def update(self, s102_obj):
        raise NotImplementedError(f"Haven't implemented the upgrade of existing data yet")

