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


class S102FeatureInformation(v2_0.FeatureInformation):

    def datatype_create(self):
        # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.datatype = self.__datatype_type__("H5T_FLOAT")


class S102FeatureInformationDataset(v2_0.S102FeatureInformationDataset):
    @property
    def metadata_type(self) -> Type[S102FeatureInformation]:
        return S102FeatureInformation


class TrackingListCoverageDataset(S102FeatureInformationDataset):
    @property
    def metadata_name(self) -> str:
        return TRACKING_COVERAGE


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
    def update(self, s102_obj):
        raise NotImplementedError(f"Haven't implemented the upgrade of existing data yet")

    def set_defaults(self, overwrite=True) -> S102File:
        """ Creates or updates an S102File object.
        Default values are set for any data that don't have options or are mandatory to be filled in the S102 spec.

        Parameters
        ----------
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        overwrite
            If updating an existing file then set this option to False in order to retain data (not sure this is needed).

        Returns
        -------
        S102File
            The object created or updated by this function.


        """
        # @fixme @todo -- I think this will overwrite no matter what, need to look into that
        self.create_empty_metadata()  # init the root with a fully filled out empty metadata set
        root = self.root
        bathy_cov_dset = root.feature_information.bathymetry_coverage_dataset
        bathy_depth_info = bathy_cov_dset.append_new_item()  # bathy_cov_dset.append(bathy_cov_dset.metadata_type())
        bathy_depth_info.initialize_properties(True, overwrite=overwrite)
        bathy_depth_info.code = DEPTH
        bathy_depth_info.name = DEPTH
        # these are auto-filled by the api:
        # unit_of_measure
        # fill_value
        # datatype
        # lower
        # upper
        # closure

        bathy_uncertainty_info = bathy_cov_dset.append_new_item()
        bathy_uncertainty_info.initialize_properties(True, overwrite=overwrite)
        bathy_uncertainty_info.code = UNCERTAINTY
        bathy_uncertainty_info.name = UNCERTAINTY

        root.bathymetry_coverage.axis_names = numpy.array(["Longitude", "Latitude"])  # row major order means X/longitude first
        root.bathymetry_coverage.sequencing_rule_scan_direction = "Longitude, Latitude"
        root.bathymetry_coverage.common_point_rule = 1  # average
        # root.bathymetry_coverage.data_coding_format = 2  # default
        # root.bathymetry_coverage.dimension = 2  # default value
        root.bathymetry_coverage.interpolation_type = 1  # nearest neighbor
        root.bathymetry_coverage.num_instances = 1  # how many Bathycoverages
        root.bathymetry_coverage.sequencing_rule_type = 1  # linear
        del root.bathymetry_coverage.time_uncertainty

    def load_bag(self, bagfile, output_file, metadata: dict = None) -> S102File:
        """
        Parameters
        ----------
        bagfile
            Either a path to a raster file that GDAL can open or a gdal.Dataset object.
        output_file
            Can be an S102File object or anything the h5py.File would accept, e.g. string file path, tempfile obect, BytesIO etc.
        metadata
            Supports the metadata options in :any:`from_from_arrays_with_metadata`.
            In addition, 'resample_resolution' can supplied to use a particular resolution using gdal "MODE=RESAMPLED_GRID"
        Returns
        -------

        """
        # @todo update method docstring for possible metadata fields
        if metadata is None:
            metadata = {}
        else:
            metadata = metadata.copy()

        if isinstance(bagfile, gdal.Dataset):
            bag = bagfile
        else:
            bag = gdal.Open(bagfile)

        # check for and resample variable resolution BAG if able
        gdal_metadata = bag.GetMetadata()
        if 'HAS_SUPERGRIDS' in gdal_metadata and gdal_metadata['HAS_SUPERGRIDS'] == 'TRUE':
            bag_filename = bag.GetFileList()[0]
            if "resample_resolution" in metadata:
                res = metadata["resample_resolution"]
                bag = None
                bag = gdal.OpenEx(bag_filename, open_options=['MODE=RESAMPLED_GRID', f'RESX={res}', f'RESY={res}'])
            else:
                warnings.warn(f'No resampling resolution provided for variable resolution bag {bag_filename}.  Using overview resolution.',
                              category=RuntimeWarning)

        # populate the issueDate if possible from a simple string search
        if 'issueDate' not in metadata:
            xml_str = bag.GetMetadata('xml:BAG')[0]
            root = et.fromstring(xml_str)
            elem = root.find(".//" + gco + "Date")
            if elem is not None and elem.text:
                metadata['issueDate'] = elem.text

        self.load_gdal(bag, metadata=metadata, flip_z=True)
