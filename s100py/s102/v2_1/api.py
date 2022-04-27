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
Change FeatureInformation datatype to H5T_FLOAT from H5T_NATIVE_FLOAT - per table 10-3 
"""

del DisplayScaleMixin
del TrackingListValues
del TrackingListValuesList
del TrackingListSetList
del TrackingListCoverage
del TrackingListGroupList
del TrackingListCoveragesList
del TrackingListCoverageDataset
del FeatureCodesTrackingMixin
del S102RootTrackingMixin


# removed the tracking list mixin
class FeatureCodes(v2_0.FeatureCodesBase):
    def feature_name_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_name = self.__feature_name_type__([BATHY_COVERAGE, ], dtype=h5py_string_dtype)

    def feature_code_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.feature_code = self.__feature_code_type__([BATHY_COVERAGE, ], dtype=h5py_string_dtype)


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


class S102FeatureInformation(v2_0.FeatureInformation):

    def datatype_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.datatype = self.__datatype_type__("H5T_FLOAT")


class S102FeatureInformationDataset(v2_0.S102FeatureInformationDataset):
    @property
    def metadata_type(self) -> Type[S102FeatureInformation]:
        return S102FeatureInformation


class BathymetryCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return BATHY_COVERAGE


class S102Root(v2_0.S102RootBase):
    @property
    def __feature_information_type__(self):
        return FeatureCodes

    @property
    def __bathymetry_coverage_type__(self):
        return BathymetryContainer


class S102File(v2_0.S102File):
    PRODUCT_SPECIFICATION = numpy.string_('INT.IHO.S-102.2.1')

    def __init__(self, name, *args, **kywrds):
        if 'root' not in kywrds:
            kywrds['root'] = S102Root  # inherited classes will specify their own root type
        super().__init__(name, *args, **kywrds)

    @property
    def z_down(self) -> bool:  # reverse Z direction
        return True

    def upgrade(self, s102_obj):
        raise NotImplementedError(f"Haven't implemented the upgrade of existing data yet")

    def set_defaults(self, overwrite=True):  # remove tracking list
        self.create_empty_metadata()  # init the root with a fully filled out empty metadata set
        self._set_bathy_defaults()
