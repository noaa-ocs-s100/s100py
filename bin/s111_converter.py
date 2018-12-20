#!/usr/bin/env python
"""Convert S-111 compliant HDF5 File(s).

S-111 is an IHO standard outlining formats for storing and sending surface
water current data and metadata.
"""
import gdal
import osr
import argparse
from glob import glob
import os

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class S111Converter:
    """Convert S111 HDF5 File(s)."""

    def toGeotiff(input_path, output_path):
        """Create a 2-Band GeoTIFF for every speed and direction compound dataset
           within each HDF5 file(s).
        
        Args:
            input_path: Path to a single S-111 HDF5 file or a directory containing
                one or more.
            output_path: Path to a directory where GeoTIFF file(s) will be
                generated.
        """
        # Creates a list of all HDF5 file(s) specified in input
        if input_path.endswith(".h5"):
            hdf5_files = [input_path]
        else:
            hdf5_files = glob("{}/*.h5".format(input_path))

        for file in hdf5_files:
            with h5py.File(file, "r") as h5_file:
                # Read S111 HDF5 feature instance, attributes and values
                feature_instance = h5_file["/SurfaceCurrent/SurfaceCurrent.01/"]
                num_grp = feature_instance.attrs['numGRP']
                split_path = os.path.split(file)
                filename = os.path.splitext(split_path[1])

                for idx in range(1, num_grp + 1):
                    values = feature_instance["Group_{:03d}/values".format(idx)]
                    speed = values['surfaceCurrentSpeed']
                    direction = values['surfaceCurrentDirection']
                    datetime = feature_instance["Group_{:03d}".format(idx)].attrs['timePoint'][0:16]

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
                    name = "{}{}_{}.tif".format(output_path, filename[0], datetime)
                    new_dataset = gdal.GetDriverByName('GTiff').Create(name, x_dim, y_dim, num_bands, gdal.GDT_Float32)
                    new_dataset.SetGeoTransform(geoTransform)
                    new_dataset.SetProjection(srs.ExportToWkt())

                    new_dataset.GetRasterBand(1).WriteArray(speed)
                    new_dataset.GetRasterBand(1).SetDescription('speed')
                    new_dataset.GetRasterBand(1).SetNoDataValue(-9999.0)
                    new_dataset.GetRasterBand(2).WriteArray(direction)
                    new_dataset.GetRasterBand(2).SetDescription('direction')
                    new_dataset.GetRasterBand(2).SetNoDataValue(-9999.0)
        print ("Conversion Complete")


def main():
    """Parse command line arguments and execute target functions."""
    parser = argparse.ArgumentParser(description="Convert S-111 HDF5(s) to 2-Band GeoTIFFs")
    parser.add_argument("-i", "--input_path", help="Path to a single HDF5 (*.h5) file or a directory where multiple HDF5 file(s) are located.")
    parser.add_argument("-o", "--output_path", help="Path to a directory where GeoTIFF file(s) will be generated.")
    args = parser.parse_args()

    if not args.input_path:
        parser.error("Path to input file(s) (--input_path) must be specified.")
        return 1
    if not args.output_path:
        parser.error("GeoTIFF output path (--output path) must be specified.")
        return 1

    S111Converter.toGeotiff(args.input_path, args.output_path)

    return 0


if __name__ == "__main__":
    main()
