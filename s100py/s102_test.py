import os
import logging

from HSTB.drivers.s100 import s102, bag_to_s102

if __name__ == "__main__":

    logging.basicConfig(level=logging.DEBUG)
    if 0:  # read and print an existing file using plain HDF5
        convertor_path = r"C:\Git_Repos\BagToS102\x64\Release\BAG_to_S102.exe"
        bag_path = r"C:\Git_Repos\BagToS102\x64\Debug\LA_LB_Area_GEO.bag"
        bag_path = r"C:\downloads\S102\S102__linux_from_NAVO\BAG_to_S102_converter\sample_data\LA_LB_Area_UTM_original.bag"
        output_path = r"C:\Git_Repos\BagToS102\x64\Release\test_output"

        convert = False
        if convert:
            print("converting", bag_path, "to", output_path)
            S102File.convert_bag(bag_path, output_path, convertor_path)
            print("finished")

        paths = [  # r"C:\downloads\S102\S102__linux_from_NAVO\BAG_to_S102_converter\sample_data\LA_LB_Area_UTM_original.bag_%d.h5",
            # r"C:\downloads\S102\s102_first_from_george\sample_output\LA_LB_AREA_UTM_S102_%d.h5",
            # r"C:\downloads\S102\s102_third_from_george\x64\Debug\LA_LB_Area_GEO_reprojected.bag_%d.h5",
            # r"C:\downloads\S102\s102_third_from_george\x64\Debug\BATHY_GEN_resolution_10m.bag_%d.h5",
            # output_path + "_%d.h5",
            # r"C:\downloads\S102\sample_hdf5_data\h5ex_t_bitatt.h5",
            # r"C:\Data\S102 Data\BATHY_GEN_resolution_10m_%d.h5",
            r"C:\Data\S102 Data\reproj_s102_%d.h5",
            # r"C:\Git_Repos\BagToS102\x64\Release\test_output_%d.h5",
        ]
        show_hdf5 = True
        if show_hdf5:
            # r"C:\downloads\S102\S102__linux_from_NAVO\BAG_to_S102_converter\sample_data\LA_LB_Area_UTM_original.bag_%d.h5"
            for fname in paths:
                print(fname)
                f = S102File(fname)
                f.print_overview()
                f.print_depth_attributes()
                f.show_keys(f)
                # for k in S102File.top_level_keys:
                #    if k in f:
                #        d = f[k]
                # d = f["Group_F"]["BathymetryCoverage"]
        if 1:  # read the test file using the python driver
            print("Python Driver not tested yet")
    if 0:  # create a new test file
        # test read
        try:
            fname  # use the last fname from above?
        except:
            fname = r"C:\Git_Repos\BagToS102\x64\Release\test_output_%d.h5"  # r"C:\Data\S102 Data\reproj_s102_%d.h5"
        f = S102File(fname)
        f.read_s102_metadata()

        # create a second file and write the info from the file we just read into it.
        # add _0.h5 on the end since it's writing as a 'family' file currently
        f2 = S102File(fname + ".write.test_0.h5", "w")
        f.root.write(f2)

        # test write
        root = Root()
        feat_codes = FeatureCodes()
        feat_codes.feature_name = "Test external creation and setter"
        if 0:  # haven't done write yet
            root.feature_information = feat_codes
            root.write(None)
            print()
            root.feature_information_create()
            root.feature_information.feature_name = "Test via _create and getter"
            root.write(None)
        print()
    if 1:
        try:
            bagname
            output_path = ""
        except:
            bagname = r"C:\Data\S57_and_S100_testData\S102\Band5_LALB_4m.bag"
            if not os.path.exists(bagname):
                bagname = r"C:\Data\S102 Data\Band5_LALB_4m.bag"
            output_path = bagname + ".s102_output_test.h5"

        convert = True
        if convert:
            try:
                os.remove(output_path)
            except FileNotFoundError:
                pass
            bag_to_s102.sr_bag_to_s102(bagname, output_path=output_path)

        s102_read_test = s102.S102File(output_path, "r", driver=None)
        s102_copy_root_test = s102.S102File(output_path + "copyroot.h5", "w", driver=None)
        s102_copy_root_test.root = s102_read_test.root
        s102_copy_root_test.write_s102_data()

        # s102_read_test.read_s102_metadata()

    raise Exception("""Questions to answer/raise ----

    Does data have to be float32?
    GeographicExtent is a parent class for GeographicBoundingBox, does NAVO treat it that way
    In BathymetryCoverage (the one with values) the NAVO is missing min,max depth+uncertainty.
        Not sure what to put for min/max displayscale (see table 11.1)
        NAVO also missing origin, offsetvectors, dimension, extent, sequenceRule, startSequence
        NAVO does have axisName (mismatch between figure 4.4 and section 4.2.1.1.1.1 -- plural or not) -- in S100 Fig 8-23 and following text it is axisNames (plural)
        axisNames says "Sequence<CharacterString>"  this means dataset or attribute?
    origin in NAVO is a string but should be DirectPosition class (group)
    CV_SequenceRule 4.2.1.1.7 should be a sub group, per Figure 4.4 of S102?
    scanDirection is missing documentation or is supposed to by under SequenceRule 4.2.1.1.8
    S100 10c-9.3 says   Note that numROWS, numCOLS, numZ, and numPOS
    Annex B in the S102 v2.0.0 is based on NOAA's data - not from spec but an output from the NAVO program.
    originalDepth is not listed in the S102 or S100.  I think NAVO is interpretting wrong, should leave original depths in the array 
        and override with the tracking list values specified - per my read of the spec.
    BathyCoverage has gridMatrix containing Value arrays - not sure NAVO has that container Figure 4.4
    Table 10.1 is showing BathyCoverages as under level 2 or 3 but the UML diagram in Figure 4.4 shows a gridMatrix with 'values' under it -
        where does that fit?
    Inconsistent use of UML, text, tables is confusing and leads to error.
    Figure 4.3 shows gridMatrix/values for Bathy and set/point for TrackingList coverages in the Coverage types - implies extra level compared to NAVO and maybe other tables (like 10.1) in doc
    rangeType in the TrackingListCoverage is of type "RecordType" but I don't understand exactly what that is.  Table 4.1 and sec 4.2.1.1.9.1 reference this.  Is this the "set"?
    BathymetricCoverage uses "Sequence" but need to find it's definition.  HDF5 dataset?
    In S100 squence is referred to from a GML perspective, which is XML.  Does that apply to HDF5 for S102?
    4.2.1.1.6.2 has an angle bracket issue at sequence<integer> -- using non-standard font?
    section 4.2.1.1.2.2 references S102_DataIdentification class which is not listed in the docs (at least not in the main sections)
        there is a MD_DataIdentification in chapter 12 about Metadata
        It appears that the docs are saying the S102_DataIdentification in the overall metadata object supplies the meaning of the depths (datum? etc)
        so there is no place that it is used in the BathymetryValues classes.
    Why is NAVO making a nested array for FeatureInformation, is it one per bathycoverage? I think it should be shape = (2,) not (2,1)
    In FeatureInformation, does fillValue need to be a string or can it be a float?
    Do datasets (e.g. FeatureInformation) need to be ordered -- we should support writing in the order of the S102 doc but read in any order.
    Section 10.2.1 and Table 10.1 have a mismatch.  The table shows "featureCode" but the text says "featureName" three times.  NAVO used featureCode.
    featureName (or featureCode depending on answer from above) only lists BathymetryCoverage and TrackingListCoverage (sec 10.2.1) but NAVO has another array of Latitude, Longitude in the same dataset
    Why does 10.2.4 talk about two dimensional arrays while section 5.2 talks about a one dimension read using C columns so sample C+1 is row=2
    S102_BathymetryValueRecord is referenced but not defined.  
    multiple items are referred to in text but not explicitly listed.  verticalUncertaintyType, depthCorrectionType, numPointsLongitudinal and numPointsLatitudinal as examples (some listed in S100)
    Why does figure 4.4 call things gridMatrix 1  and values 1..*  and set 0..* and point 1..*  and never mention 'Group" then chapter 10 talks about Group.NNN
    and 4.2.1.1.2.2 - values talks like it's a sequential list of values and not a NxM matrix dataset?
    Is there a mismatch in fig 4.4 and the 4.2.1.1.1.1 text, one says maximumElevation and the other says maximumDepth.
    sequencignRule.type and sequencingRule.scanDirection -- literal or are they saying use the sequencingRule class.  S111 does literal, NAVO uses just after the dot.
    """)
    raise NotImplementedError("S102 says it doesn't expect TilingScheme inside the data currently, implment after cleaning up initial data tests")
