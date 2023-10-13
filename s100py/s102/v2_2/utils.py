""" Functions to create S102 data from other sources

If the utils module is run as __main__ then it will run  :func:`from_bag` or :func:`from_gdal` on the filename given in the command line arguments.

"""
# this also works but doesn't put parenthesis in the html   :any:`from_bag`
# fully qualified would also work   :func:`s100py.s102.utils.from_gdal`

import os
import sys

if getattr(sys, 'frozen', False):
    # in a frozen exe os.pyc is at the root level
    proj_db_path = os.path.join(os.path.dirname(os.__file__), "library\\share\\proj")
    # print("frozen exe", proj_db_path)
    os.environ["PROJ_LIB"] = proj_db_path
    # print(os.listdir(proj_db_path))  # os.listdir(os.path.dirname(os.__file__)))
else:
    pass
    # print("running as script")

import logging
import warnings
import argparse
from xml.etree import ElementTree as et
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy
from osgeo import gdal, osr
import h5py
import re

try:
    from matplotlib import pyplot
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
    from matplotlib.colorbar import ColorbarBase
except:
    if not getattr(sys, 'frozen', False):  # we expect the frozen exe to not have matplotlib
        print("matplotlib.pyplot failed to import, plotting will not work")

from s100py.s102.v2_2.api import DEPTH, UNCERTAINTY, S102File, S102Exception

__all__ = ['plot_depth_using_h5py', 'create_s102', 'from_arrays', 'from_arrays_with_metadata',
           'from_gdal', 'from_bag', 'get_valid_epsg', 'to_geotiff']

gco = "{http://www.isotc211.org/2005/gco}"

# @todo create a friendly name mapping to s102 nested location, then add s102 functions for "to dictionary" and "from dictionary" to api
#   that would make these functions easily invertable

r"""
from s100py.s102 import utils
navo_name = r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag.navo_%d.h5"
utils.plot_depth_using_h5py(navo_name)
bag_names = [r"G:\Data\S102 Data\NBS_US5NYCAH_20200430.bag", r"G:\Data\S102 Data\NBS_US5NYCBH_20200429.bag", r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag"]
fout = utils.from_bag(bag_name, bag_name+r".noaa.h5")
for bag_name in bag_names:
    fout = utils.from_bag(bag_name, bag_name+r".noaa.h5")
    fout.close()
    utils.plot_depth_using_h5py(bag_name + r".noaa.h5")
"""

def plot_depth_using_h5py(filename, enc_color=False):
    # filename = r"G:\Data\S102 Data\GlenS102Test\102USA15NYCAH200430.H5"
    # h5py.File(r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag_%d.h5", mode="r", driver="family", memb_size=681574400)
    try:
        f = h5py.File(filename, mode="r", driver="family")
    except OSError as e:
        # see if the member size isn't right and then fall back to standard opening
        # OSError: Unable to open file (Family member size should be 681574400.  But the size from file access property is 2147483647)
        try:
            error_string = str(e)
            m = re.search(r'Family member size should be\s+(\d+)', error_string)
            if m:
                sz = int(m[1])
            f = h5py.File(filename, mode="r", driver="family", memb_size=int(sz))
        except:
            f = h5py.File(filename, mode="r")
    fill_val = 1000000
    try:
        d = f["BathymetryCoverage/BathymetryCoverage.01/Group.001/values"]['depth']
    except KeyError:
        try:
            d = f["BathymetryCoverage/BathymetryCoverage.001/Group.001/values"]['depth']
        except KeyError:
            d = f["SurfaceCurrent/SurfaceCurrent.01/Group_001/values"]['surfaceCurrentSpeed']
            fill_val = -9999
    d[d==fill_val] = numpy.nan

    # ud = numpy.flipud(d)
    if enc_color:
        colors = numpy.array([(255, 255, 255), (201, 237, 252), (167, 217, 251), (130, 202, 255), (97, 180, 255)]) / 255
        # bins = [-1000, -9.1, -5.4, -3.6, -1.8, ]
        bins = [-1000, -30, -22, -15, -12]
        norm = BoundaryNorm(bins, colors.shape[0])
        cmap = ListedColormap(colors)
        pyplot.imshow(d, interpolation='nearest', cmap=cmap, norm=norm)
    else:
        pyplot.imshow(d)
    pyplot.gca().invert_yaxis()
    pyplot.show()
    # cm = LinearSegmentedColormap.from_list("enc", colors, N=4)
    # im = pyplot.imshow(ud, cmap=cm, interpolation='nearest')
    # pyplot.colorbar(im)


def browse_files(question):
    # using tkinter since it is built in to python and smaller to distribute than PySide2 or wxPython in an executable
    root = tk.Tk()
    root.withdraw()
    # root.filename = tkFileDialog.askopenfilename(initialdir="/", title="Select file",
    #                                              filetypes=(("jpeg files", "*.jpg"), ("all files", "*.*")))
    file_path = filedialog.askopenfilename(title=question)
    return file_path


def bool_question(question, title="", icon="warning"):
    root = tk.Tk()
    root.withdraw()
    result = messagebox.askquestion(title, question, icon=icon)
    return result == 'yes'


def make_parser():
    parser = argparse.ArgumentParser(description='Convert a georeferenced file to S102')
    parser.add_argument("-?", "--show_help", action="store_true", help="show this help message and exit")
    parser.add_argument("-i", "--input_filename", help="full path to the file to be processed")
    parser.add_argument("-o", "--output_filename", help="output filename, default is same name as input with .h5 appended")
    parser.add_argument("-r", "--res", help="Resolution.  If the input file is a BAG then use attempt to use the given resolution" )
    return parser


def to_geotiff(input_path, output_path):
    s102_data = S102File(input_path)
    s102_data.to_geotiff(output_path)


create_s102 = S102File.create_s102
from_arrays = S102File.from_arrays
from_arrays_with_metadata = S102File.from_arrays_with_metadata
from_gdal = S102File.from_raster
from_raster = S102File.from_raster
from_bag = S102File.from_bag
get_valid_epsg = S102File.get_valid_epsg

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    if args.show_help:
        parser.print_help()
        sys.exit()

    if not args.input_filename:
        path = browse_files('Choose file to convert to S102')
        if path:
            args.input_filename = path

    if args.input_filename:
        # get the output file path from the command line or the user if it wasn't specified
        if args.output_filename:
            output_name = args.output_filename
        else:
            output_name = args.input_filename + ".h5"
            if os.path.exists(output_name):
                if bool_question(f"{output_name} already exists, overwrite?", "Overwrite File"):
                    os.remove(output_name)

        # check if the data is a bag and should be sent to the from_bag function or just raster and send to from_gdal
        ds = gdal.Open(args.input_filename)
        drv = ds.GetDriver()
        if drv.GetDescription() == "BAG":
            from_bag(ds, output_name)
        else:
            from_gdal(ds, output_name)
