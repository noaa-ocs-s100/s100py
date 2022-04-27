from s100py.s102.v2_0 import utils as v2_0
from s100py.s102.v2_0.utils import *
from s100py.s102.v2_1.api import DEPTH, UNCERTAINTY, S102File, S102Exception

create_s102 = S102File.create_s102
from_arrays = S102File.from_arrays
from_arrays_with_metadata = S102File.from_arrays_with_metadata
from_gdal = S102File.from_gdal
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

