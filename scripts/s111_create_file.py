#!/usr/bin/env python

import argparse
import h5py
import numpy
import csv
import os

def clear_metadata_value(attributes, attribute_name):
    """ Clear the specified attribute value.

    :param attributes: The list of attributes containing the value to be cleared.
    :param attribute_name: The name of the attribute to be cleared.
    """

    if attribute_name in attributes:
        print("Information: The value for", attribute_name, "has been ignored.")
        del attributes[attribute_name]

#******************************************************************************
def get_metadata_type(attribute_name):
    """ Retrieve the specified attribute's type.

    :param attribute_name: The name of the attribute to retrive the type for.
    :returns: The attribute's type, None if not found.
    """

    typeMap = dict()
        
    """
         Carrier Metadata
    """
    #Integer types
    typeMap['horizDatumValue'] = numpy.int64 
    typeMap['timeRecordInterval'] = numpy.int64
    typeMap['numberOfTimes'] = numpy.int64
    typeMap['numberOfStations'] = numpy.int64
    typeMap['numPointsLongitudinal'] = numpy.int64
    typeMap['numPointsLatitudinal'] = numpy.int64
    typeMap['minGridPointLongitudinal'] = numpy.int64
    typeMap['minGridPointLatitudinal'] = numpy.int64
    typeMap['numberOfNodes'] = numpy.int64
    #Real types
    typeMap['surfaceCurrentDepth'] = numpy.float64
    typeMap['gridOriginLongitude'] = numpy.float64
    typeMap['gridOriginLatitude'] = numpy.float64
    typeMap['gridSpacingLongitudinal'] = numpy.float64
    typeMap['gridSpacingLatitudinal'] = numpy.float64
    typeMap['gridLandMaskValue'] = numpy.float64
    typeMap['speedUncertainty'] = numpy.float64
    typeMap['directionUncertainty'] = numpy.float64
    typeMap['horizontalPositionUncertainty'] = numpy.float64
    typeMap['verticalUncertainty'] = numpy.float64
    typeMap['timeUncertainty'] = numpy.float64
    typeMap['minDatasetCurrentSpeed'] = numpy.float64
    typeMap['maxDatasetCurrentSpeed'] = numpy.float64

    #String types
    typeMap['productSpecification'] = numpy.bytes_
    typeMap['dateTimeOfIssue'] = numpy.bytes_
    typeMap['nameRegion'] = numpy.bytes_
    typeMap['nameSubregion'] = numpy.bytes_
    typeMap['horizDatumReference'] = numpy.bytes_
    typeMap['dateTimeOfFirstRecord'] = numpy.bytes_
    typeMap['dateTimeOfLastRecord'] = numpy.bytes_ 
    typeMap['methodCurrentsProduct'] = numpy.bytes_

    #Enumeration types
    typeMap['typeOfCurrentData'] = numpy.int64
    typeMap['dataCodingFormat'] = numpy.int64
    typeMap['depthTypeIndex'] = numpy.int64
    typeMap['verticalDatum'] = numpy.int64


    if attribute_name not in typeMap:
        return None
        
    return typeMap[attribute_name]
    
#******************************************************************************
def add_metadata(attributes, metadata_file):
    """ Add metadata values to the S-111 attributes.

    :param attributes: The S-111 attributes to be populated.
    :param metadata_file: The ASCII CSV file to retrieve the metadata values from.
    """

    with open(metadata_file) as csvfile:
        reader = csv.reader(csvfile)
        
        #Grab the header and data rows.
        header = next(reader)
        data = next(reader)
        
        colnum = 0
                
        #For each column in the data row...
        for col in data:
            attribute_name = header[colnum].strip()
            attribute_value = col.strip().encode()
            attribute_type = get_metadata_type(attribute_name)
            
            #If we don't know what this attribute is, just report it to the user.
            if attribute_type == None:
                print("Warning: Unknown metadata value", attribute_name)
            #Else if this is a string type...
            elif attribute_type == numpy.bytes_:
                attributes.create(attribute_name, attribute_value)
            #Else use the type returned.
            else:
                attributes.create(attribute_name, attribute_value, dtype=attribute_type)
                
            colnum += 1


    #Attribute information is derived from netcdf files 
    clear_metadata_value(attributes, 'dateTimeOfFirstRecord')
    clear_metadata_value(attributes, 'dateTimeOfLastRecord')
    clear_metadata_value(attributes, 'numberOfStations')
    clear_metadata_value(attributes, 'numberOfTimes')
    clear_metadata_value(attributes, 'dataCodingFormat')
    clear_metadata_value(attributes, 'timeRecordInterval')
    clear_metadata_value(attributes, 'minDatasetCurrentSpeed')
    clear_metadata_value(attributes, 'maxDatasetCurrentSpeed')
    

    attributes.create('numberOfStations', 0, dtype=numpy.int64)
    attributes.create('numberOfTimes', 0, dtype=numpy.int64)

#******************************************************************************    
def create_dataset(output_file, metadata_file):
    """ Create a new S-111 dataset.

    :param output_file: The name of the file to be created.
    :param metadata_file: The ASCII CSV file to retrieve the metadata values from.
    """

    #Make sure the output file has the correct extension.
    filename, file_extension = os.path.splitext(output_file)
    output_file_with_extension = filename + ".h5"

    #Create the new HDF5 file.
    with h5py.File(output_file_with_extension, "w") as hdf_file:
    
        #Add the metadata to the file.
        add_metadata(hdf_file.attrs, metadata_file)
        
#******************************************************************************        
def create_command_line():
    """Create and initialize the command line parser.
    
    :returns: The command line parser.
    """

    parser = argparse.ArgumentParser(description='Create S-111 File')

    parser.add_argument('-m', '--metadata-file', help='The text file containing the file metadata.', required=True)
    parser.add_argument("outputFile", nargs=1)

    return parser

#******************************************************************************        
def main():

    #Create the command line parser.
    parser = create_command_line()

    #Parse the command line.
    results = parser.parse_args()
    
    create_dataset(results.outputFile[0], results.metadata_file)



if __name__ == "__main__":
    main()
