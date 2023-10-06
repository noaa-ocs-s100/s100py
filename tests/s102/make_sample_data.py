import os
import pathlib
import logging

from osgeo import gdal, osr
from s100py import s100
from s100py.s102 import v2_0
from s100py.s102 import v2_1
from s100py.s102 import v2_2

def remove_file(pth, quiet=True):
    """Remove a file if it exists"""
    try:
        os.remove(pth)
    except (FileNotFoundError, PermissionError):  #
        if not quiet:
            logging.warning(f"{pth} not found")


# test_rat(str(local_path.joinpath("F00788_SR_8m.tif")), output_path)
# test_rat(tiffname, output_path)
if 0:
    metadata = {"horizontalDatumReference": "EPSG", "horizontalDatumValue": 32610}
    out_path = r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102b.2_2.h5"
    remove_file(out_path)
    new_s102_20 = v2_2.utils.from_gdal(r"C:\Data\BlueTopo\RATs\BlueTopo_BC25M26L_20221102b.tiff", out_path, metadata=metadata)
if 1:
    root_path = pathlib.Path("C:\\Data\\S102\\S102_v2.2\\")
    bags = [r"NBS_US4NY1CQ_20230221.bag",
            r"NBS_US5NYCII_20230602.bag",
            r"NBS_US5NYCJH_20230110.bag",
            r"NBS_US5NYCJI_20230221.bag",
            r"NBS_US5NYCKG_20230110.bag",
            r"NBS_US5NYCKH_20230110.bag",
            r"NBS_US5NYCLG_20230110.bag",
            r"NBS_US5NYCLH_20230221.bag",
            r"NBS_US5NYCMG_20230110.bag",
            r"NBS_US5NYCMH_20230110.bag",
            r"NBS_US5NYCNG_20230110.bag",
            r"NBS_US5NYCNH_20230221.bag",
            r"NBS_US5NYCOG_20230110.bag",
            r"NBS_US5NYCGH_20230602.bag",
            r"NBS_US5NYCHH_20230602.bag",
            r"NBS_US5NYCIH_20230602.bag",]

    from urllib.request import urlretrieve

    # download the savannah files, if needed
    # for line in savanah_files.split("\n"):
    #     try:
    #         tif_url, aux_url = line.split(", ")
    #     except ValueError:
    #         continue
    #     else:
    #         for url in (tif_url, aux_url):
    #             filename = os.path.split(url)[-1]
    #             target = root_path.joinpath(filename)
    #             tiffs.append(target)
    #             if not target.exists():
    #                 urlretrieve(url, target)
    # Used the bluetopo downloader since the dates had changed invalidating the links for the savannah files
    for tif in root_path.glob("*.tif?"):
        print('processing', tif)
        # Bluetopop TIFs were using NAD83 (epsg~26900) instead of WGS84 (epsg~32600)
        dataset = gdal.Open(str(tif))
        sr = osr.SpatialReference(dataset.GetProjection())
        epsg = sr.GetAuthorityCode(None)
        # It also was using compound CRS with custom vertical
        if epsg is None and sr.IsProjected():
            crs_2d = osr.SpatialReference(dataset.GetProjection())
            # in GDAL 3.2 this was added which may make getting the horizontal CRS more obvious
            crs_2d.DemoteTo2D()
            epsg = crs_2d.GetAuthorityCode(None)
        if 26900 <= int(epsg) < 26999:
            epsg = 32600 + int(epsg) - 26900
        metadata = {'geographicIdentifier': "Sample Data", "horizontalDatumReference": "EPSG", "horizontalDatumValue": epsg}
        fname_2_2 = tif.with_suffix(".2_2.h5")
        remove_file(fname_2_2)
        s102_22 = v2_2.utils.from_gdal(tif, fname_2_2, metadata=metadata, flip_z=True)  # bluetopo tiff is in elevation instead of depth
        fname_2_1 = tif.with_suffix(".2_1.h5")
        remove_file(fname_2_1)
        s102_21 = v2_1.utils.from_gdal(tif, fname_2_1, metadata=metadata, flip_z=True)  # bluetopo tiff is in elevation instead of depth

if 0:
    metadata = {'geographicIdentifier': "Sample Data"}
    for bag in root_path.glob("*.bag"):
        fname_2_2 = bag.with_suffix(".bag.2_2.h5")
        remove_file(fname_2_2)
        s102_22 = v2_2.utils.from_bag(bag, fname_2_2, metadata=metadata)
        fname_2_1 = bag.with_suffix(".bag.2_1.h5")
        remove_file(fname_2_1)
        s102_21 = v2_1.utils.from_bag(bag, fname_2_1, metadata=metadata)

if 0:
    # test upgrade of v2.1 to v2.2
    remove_file(r"C:\data\S102\S102_v2.2\NBS_US4NY1CQ_20230221_162221_base.tiff.2_1_to_2_2.h5")
    sf22 = v2_2.api.S102File.upgrade(r"C:\data\S102\S102_v2.2\NBS_US4NY1CQ_20230221_162221_base.tiff.2_1.h5", r"C:\data\S102\S102_v2.2\NBS_US4NY1CQ_20230221_162221_base.tiff.2_1_to_2_2.h5")
