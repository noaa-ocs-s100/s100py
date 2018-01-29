#!/usr/bin/env python

import argparse
import h5py
import numpy
import os

"""
OFS S-100/S-111 HDF5 File 

+ Create HDF5 File
+ Create Empty Global attributes

"""
def add_metadata(hdf_file):
    """ Add metadata to HDF5 file."""

    hdf_file.attrs.create('horizDatumValue', 0 , dtype=numpy.int32) 
    hdf_file.attrs.create('timeRecordInterval', 0 , dtype=numpy.int32)
    hdf_file.attrs.create('numberOfTimes', 0 , dtype=numpy.int32)
    hdf_file.attrs.create('numberOfStations', 0 , dtype=numpy.int32)
    hdf_file.attrs.create('numPointsLongitudinal', 0 , dtype=numpy.int32)
    hdf_file.attrs.create('numPointsLatitudinal', 0 , dtype=numpy.int32)
    hdf_file.attrs.create('numberOfNodes', 0 , dtype=numpy.int32)
   
    #Real types
    hdf_file.attrs.create('surfaceCurrentDepth', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('gridOriginLongitude', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('gridOriginLatitude', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('gridSpacingLongitudinal', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('gridSpacingLatitudinal', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('minGridPointLongitudinal', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('minGridPointLatitudinal', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('gridLandMaskValue', -9999.0 , dtype=numpy.float32)
    hdf_file.attrs.create('speedUncertainty', -1.0 , dtype=numpy.float32)
    hdf_file.attrs.create('directionUncertainty', -1.0 , dtype=numpy.float32)
    hdf_file.attrs.create('horizontalPositionUncertainty', -1.0 , dtype=numpy.float32)
    hdf_file.attrs.create('verticalUncertainty', -1.0 , dtype=numpy.float32)
    hdf_file.attrs.create('timeUncertainty', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('minDatasetCurrentSpeed', 0 , dtype=numpy.float32)
    hdf_file.attrs.create('maxDatasetCurrentSpeed', 0 , dtype=numpy.float32)

    #String types
    #Changed dtype = numpy.bytes_ to dt (variable length ascii) can't store null bytes
    dt = h5py.special_dtype(vlen=bytes)
    hdf_file.attrs.create('productSpecification', 'S-111_v1.11.0' , dtype=dt)
    hdf_file.attrs.create('dateTimeOfIssue', 'None' , dtype=dt)
    hdf_file.attrs.create('nameRegion', 'US_East_Coast' , dtype=dt)
    hdf_file.attrs.create('nameSubregion', 'Cheaspeake_Bay' , dtype=dt)
    hdf_file.attrs.create('horizDatumReference', 'EPSG' , dtype=dt)
    hdf_file.attrs.create('dateTimeOfFirstRecord', 'None' , dtype=dt)
    hdf_file.attrs.create('dateTimeOfLastRecord', 'None' , dtype=dt)
    hdf_file.attrs.create('methodCurrentsProduct', 'ROMS_Hydrodynamic_Model' , dtype=dt)

    #Enumeration types
    hdf_file.attrs.create('typeOfCurrentData', 6 , dtype=numpy.int32)
    hdf_file.attrs.create('dataCodingFormat', 2 , dtype=numpy.int32)
    hdf_file.attrs.create('depthTypeIndex', 0 , dtype=numpy.int32)
    hdf_file.attrs.create('verticalDatum', 0 , dtype=numpy.int32)

    for att in hdf_file.attrs:
        print(att)
    
    hdf_file.close()
#****************************************************************************** 

def create_dataset(output_file):
    """ Create a new S-111 dataset.
        output_file: The name of the file to be created.
    """

    #Make sure the output file has the correct extension.
    filename, file_extension = os.path.splitext(output_file)
    output_file_with_extension = filename + ".h5"

    #Create the new HDF5 file.
    with h5py.File(output_file_with_extension, "w") as hdf_file:
    
        #Add the metadata to the file.
        add_metadata(hdf_file)
        
#******************************************************************************        
def create_command_line():
    """Create and initialize the command line parser. """

    parser = argparse.ArgumentParser(description='Create S-111 File')

    parser.add_argument("outputFile", nargs=1)

    return parser

#******************************************************************************        
def main():

    #Create the command line parser.
    parser = create_command_line()

    #Parse the command line.
    results = parser.parse_args()
    
    create_dataset(results.outputFile[0])



if __name__ == "__main__":
    main()