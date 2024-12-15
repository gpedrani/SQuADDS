"""
=====================================================================
Imports
=====================================================================
"""
import pathlib 




"""
=====================================================================
Helper Functions
=====================================================================
"""

def transform_to_dxf(path_object, output_path):
    
    with path_object.open(mode="rb") as file:
        




"""
=====================================================================
Public Functions
=====================================================================
"""
def convert_circuit_files(path_to_design: str, target_format: str, output_path: str) -> None: 
"""
    Convert_circuit_files() transforms a given circuit desing, identifies its current format, and converts it 
        into a target format specified.

    Paramters:
        path_to_design (str) [Required]: a string which holds the filepath to the ciruit design you wish to convert.
            Windows Example: path_to_design = ".\path\to\my\design.file" 
            UNIX/MAC-OS Example: path_to_design = "./path/to/my/design.file" 

        output_path (str) [Optional]: The directory where final output design will be saved; if no input is given,
            the process defaults to the same directory as the input file.

        target_format (str) [Required]: the file type which you wish to convert your file to; must be one of the specified formats below.
            Supported formats:
            GDS --> (DXF, A, B, C) 
"""

    #Make an array of supported formats
    supported_formats = ["DXF"]

    #Sainity Checks for path_to_design 
    assert path_to_design != None, "Error: path_to_design was not assigned a value - please specify an input file"
    #Turn string(s) to Path objects for easy assignments and checks
    path_object = Path(path_to_design)
    assert path_object.is_dir == False, "Error: path_to_design does not point to a file - please specify a target file."
    

    #Sanity Checks for target_format
    assert target_format  != None, "Error: target_format was not assigned a value - please specify a target format!" 
    assert target_format in supported_formats, "Error: target_format is not one of the following supported formats {}".format(supported_formats)

    #Define output_path
    if output_path == None:
        output_path = path_object.parent
    else: 
        output_path = Path(output_path)

    filename, extension = path_object.stem, path_object.suffix
    assert extension not in ['gds', "gds"], "Error: the input file given is not a GDS file - this function converts GDS into the target_format." 

    #Call required fucntion for the transformation into target format
    match target_format:
        case "DXF":
            transform_to_dxf(path_object, output_path)
        case "GDS":
            transform_to_gds(path_object, output_path)
        case _:
            continue 


def fix_routemeanders(cpw_objects):


