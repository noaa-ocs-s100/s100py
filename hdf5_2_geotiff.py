#!/usr/bin/env python
"""Convert S-111 compliant HDF5 Files to 2-Band GeoTIFFS."""
import h5py
import gdal
import osr
import argparse


class H5FILE:
    """Read and Manage S111 HDF5 File."""

    def __init__(self, input_file, ofs_model):
        """Initializes h5 file and output file path.

        Args:
            input_file: Filename and path of target hdf5 file.
            ofs_model: Model abbreviation.
        """
        self.input_file = input_file
        self.model = ofs_model
        self.open()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5_file.close()

    def open(self):
        self.h5_file = h5py.File(self.input_file, "r")
        
    def close(self):
        self.h5_file.close()

    def toGeotiff(self, output_path):
        """Create a 2-Band GeoTIFF for every HDF5 speed and direction compound dataset
        
        Args:
            output_path: Path to a directory where GeoTIFF file(s) will be
                generated.
        """

        # Read S111 HDF5 feature instance, attributes and values
        feature_instance = self.h5_file["/SurfaceCurrent/SurfaceCurrent.01/"]
        num_grp = feature_instance.attrs['numGRP']

        for idx in range(1, num_grp + 1):
            values = feature_instance["Group_{:03d}/values".format(idx)]
            speed = values['SurfaceCurrentSpeed']
            direction = values['SurfaceCurrentDirection']
            datetime = feature_instance["Group_{:03d}".format(idx)].attrs['timePoint'][0:16]
            print (datetime)

            # Set image size
            x_dim = speed.shape[1]
            y_dim = speed.shape[0]

            # Set Geospatial Information
            geoTransform = []
            for i in range(6):
                geoTransform.append(0.0)
            geoTransform[0] = feature_instance.attrs['gridOriginLongitude']
            geoTransform[1] = feature_instance.attrs['gridSpacingLongitudinal']
            geoTransform[2] = 0
            geoTransform[3] = feature_instance.attrs['gridOriginLatitude']
            geoTransform[4] = 0
            geoTransform[5] = feature_instance.attrs['gridSpacingLatitudinal']

            srs = osr.SpatialReference()
            srs.SetWellKnownGeogCS("WGS84")

            num_bands = 2
            name = "{}S111US_{}_TYP2_{}T{}_{}.tif".format(output_path, str(self.model),
                                                          self.h5_file.attrs['issueDate'],
                                                          self.h5_file.attrs['issueTime'], datetime)
            new_dataset = gdal.GetDriverByName('GTiff').Create(name, x_dim, y_dim, num_bands, gdal.GDT_Float32)
            new_dataset.SetGeoTransform(geoTransform)
            new_dataset.SetProjection(srs.ExportToWkt())

            new_dataset.GetRasterBand(1).WriteArray(speed)
            new_dataset.GetRasterBand(1).SetDescription('speed')
            new_dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
            new_dataset.GetRasterBand(2).WriteArray(direction)
            new_dataset.GetRasterBand(2).SetDescription('direction')
            new_dataset.GetRasterBand(2).SetNoDataValue(-9999.0)


def main():
    """Parse command line arguments and execute target functions."""
    parser = argparse.ArgumentParser(description="Convert S-111 HDF5 to 2-Band GeoTIFFs")
    parser.add_argument("-i", "--input_file", help="Input HDF5 file and directory path.")
    parser.add_argument("-o", "--output_path", help="Path to a directory where GeoTIFF file(s) will be generated.")
    parser.add_argument("-m", "--model", help="OFS Model.")
    args = parser.parse_args()

    if not args.input_file:
        parser.error("HDF5 file --input file must be specified.")
        return 1
    if not args.output_path:
        parser.error("GeoTIFF output path --output path must be specified.")
        return 1
    if not args.model:
        parser.error("Model --model must be specified.")
        return 1
        
    with H5FILE(args.input_file, args.model) as h5_file:
        h5_file.toGeotiff(args.output_path)
    return 0


if __name__ == "__main__":
    main()