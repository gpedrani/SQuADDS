import gdspy
import klayout.db as kdb
import numpy as np
from os.path import basename, splitext, abspath, join


def merge_shapes_in_layer(gds_file, output_gds_file, layer_number):
    """
    Selects all shapes in the given layer number from the input GDS file,
    merges them together, and saves the result in the output GDS file.

    Args:
        gds_file (str): Path to the input GDS file.
        output_gds_file (str): Path to the output GDS file.
        layer_number (int): The layer number whose shapes should be merged.

    Returns:
        None
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(gds_file)

    # Find the index of the specified layer with datatype 0 (common in GDS files)
    layer_index = layout.layer(layer_number, 0)

    # Iterate over each cell in the layout
    for cell in layout.each_cell():
        # Get the shapes in the specified layer
        shapes = cell.shapes(layer_index)

        # If there are shapes, merge them together using boolean union
        if not shapes.is_empty():
            print(f"Merging shapes in layer {layer_number} in cell: {cell.name}")
            
            # Create a region for the shapes on this layer
            region = kdb.Region(shapes)
            
            # Perform the boolean union to merge shapes
            merged_shapes = region.merged()

            # Clear the original shapes
            shapes.clear()

            # Insert the merged shapes back into the layer
            shapes.insert(merged_shapes)

    # Save the modified layout to the output GDS file
    layout.write(output_gds_file)

    print(f"Merged shapes in layer {layer_number} and saved to {output_gds_file}")

def crop_top_left_rectangle(gds_file, width=300, height=100, layer_number=5, datatype=0):
    """
    Removes a 300 x 100 um rectangle from the top left of the layer_number/datatype rectangle in the GDS file.

    Args:
        gds_file (str): Path to the input GDS file.
        width (int, optional): Width of the rectangle to remove in micrometers. Defaults to 300.
        height (int, optional): Height of the rectangle to remove in micrometers. Defaults to 100.
        layer_number (int, optional): The layer number of the rectangle to crop. Defaults to 5.
        datatype (int, optional): The datatype of the rectangle to crop. Defaults to 0.

    Returns:
        None
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(gds_file)

    # Iterate over each cell in the layout
    for cell in layout.each_cell():
        # Get the shapes in the 5/0 layer
        layer_index_5_0 = layout.layer(layer_number, datatype)
        shapes_5_0 = cell.shapes(layer_index_5_0)

        # Create a region from the 5/0 layer shapes
        region_5_0 = kdb.Region(shapes_5_0)

        # Get the bounding box of the 5/0 layer shapes
        bbox_5_0 = region_5_0.bbox()

        # Calculate the coordinates for the 300 x 100 um rectangle
        top_left_x = bbox_5_0.left
        top_left_y = bbox_5_0.top
        rect_width = width / layout.dbu  # Convert um to database units
        rect_height = height / layout.dbu  # Convert um to database units

        # Create the 300 x 100 um rectangle
        rect_top_left = kdb.Box(top_left_x, top_left_y - rect_height,
                                top_left_x + rect_width, top_left_y)

        # Remove the 300 x 100 um rectangle from the 5/0 layer
        cropped_shapes = region_5_0 - kdb.Region(rect_top_left)

        # Clear the original 5/0 layer shapes
        cell.shapes(layer_index_5_0).clear()

        # Insert the cropped shapes into the 5/0 layer
        cell.shapes(layer_index_5_0).insert(cropped_shapes)

    # Write the modified layout to the GDS file
    layout.write(gds_file)

    print(f"{width} x {height} um rectangle removed from the top left of the {layer_number}/{datatype} rectangle successfully.")
    print(f"Output file: {gds_file}")

def apply_fixes(gds_file, datatype=0):
    """
    Applies the required fixes to the GDS file.

    Args:
        gds_file (str): Path to the input GDS file.
        datatype (int, optional): The new datatype value for all layers. Defaults to 0.

    Returns:
        None
    """
    # Call the method to crop a 300 x 100 um rectangle on the top left of the 5/0 rectangle
    # crop_top_left_rectangle(gds_file)

    # Call the method to add a 703 layer datatype 0 rectangle covering the 5/0 layer
    add_703_layer(gds_file)

    # Call the method to modify the GDS file datatypes
    modify_gds_datatypes(gds_file, datatype)

    # Call the method to delete layers with non-zero datatypes
    delete_non_zero_datatype_layers(gds_file.replace('.gds', '_final.gds'))

    # Merge all metal layers
    merge_shapes_in_layer(gds_file.replace('.gds', '_final.gds'), gds_file.replace('.gds', '_final.gds'), 5)


def modify_gds_datatypes(gds_file, datatype=0):
    """
    Modifies the datatype of all layers in a GDS file.

    Args:
        gds_file (str): Path to the input GDS file.
        datatype (int, optional): The new datatype value for all layers. Defaults to 0.

    Returns:
        None
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(gds_file)

    # Dictionary to store the merged layers
    merged_layers = {}

    # Iterate over each cell in the layout
    for cell in layout.each_cell():
        # Dictionary to store the layer numbers and their corresponding shapes
        layer_shapes = {}

        # Iterate over each layer in the layout
        for layer_index in layout.layer_indices():
            # Get the layer number and datatype
            layer_info = layout.get_info(layer_index)
            layer_number = layer_info.layer
            layer_datatype = layer_info.datatype

            # Retrieve the shapes in the current layer for the cell
            shapes = cell.shapes(layer_index)

            # Add the shapes to the layer_shapes dictionary
            if layer_number not in layer_shapes:
                layer_shapes[layer_number] = shapes
            else:
                layer_shapes[layer_number].insert(shapes)

        # Merge layers with the same layer number using boolean union
        for layer_number, shapes in layer_shapes.items():
            if layer_number not in merged_layers:
                # Create a new layer with the modified datatype
                merged_layer = layout.layer(layer_number, datatype)
                merged_layers[layer_number] = merged_layer
            else:
                # Use the existing merged layer
                merged_layer = merged_layers[layer_number]

            # Perform boolean union on the shapes
            merged_shapes = kdb.Region(shapes)
            cell.shapes(merged_layer).insert(merged_shapes)

        # Delete the redundant layers
        for layer_index in layout.layer_indices():
            layer_info = layout.get_info(layer_index)
            if layer_info.layer not in merged_layers:
                cell.clear(layer_index)

    # Write the modified layout to the GDS file
    file_name = gds_file.replace('.gds', '_final.gds')
    layout.write(file_name)

    print("Merged all redundant layers successfully.")
    print(f"Output file: {file_name}")

    # delete_non_zero_datatype_layers(file_name)

def delete_non_zero_datatype_layers(gds_file):
    """
    Deletes all layers with datatypes not equal to 0 from the GDS file.

    Args:
        gds_file (str): Path to the input GDS file.

    Returns:
        None
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(gds_file)

    # Iterate over each cell in the layout
    for cell in layout.each_cell():
        # Iterate over each layer in the layout
        for layer_index in layout.layer_indices():
            # Get the layer datatype
            layer_info = layout.get_info(layer_index)
            layer_datatype = layer_info.datatype

            # If the layer datatype is not 0, delete the layer
            if layer_datatype != 0:
                cell.clear(layer_index)

    # Write the modified layout to the GDS file
    layout.write(gds_file)

    print("Layers with non-zero datatypes deleted successfully.")
    print(f"Output file: {gds_file}")

    # Call the method to add a 703 layer datatype 0 rectangle covering the 5/0 layer
    # add_703_layer(gds_file.replace('_modified.gds', '_final.gds'))

def add_703_layer(gds_file):
    """
    Adds a 703 layer datatype 0 rectangle covering the 5/0 layer in the GDS file.

    Args:
        gds_file (str): Path to the input GDS file.

    Returns:
        None
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(gds_file)

    # Iterate over each cell in the layout
    for cell in layout.each_cell():
        # Get the shapes in the 5/0 layer
        layer_index_5_0 = layout.layer(5, 0)
        shapes_5_0 = cell.shapes(layer_index_5_0)

        # Create a region from the 5/0 layer shapes
        region_5_0 = kdb.Region(shapes_5_0)

        # Create a bounding box covering the 5/0 layer shapes
        bbox_5_0 = region_5_0.bbox()

        # Create a 703 layer with datatype 0
        layer_index_703_0 = layout.layer(703, 0)

        # Add a rectangle covering the 5/0 layer to the 703/0 layer
        rect_703_0 = kdb.Box(bbox_5_0)
        cell.shapes(layer_index_703_0).insert(rect_703_0)

    # Write the modified layout to the GDS file
    layout.write(gds_file)

    print("703 layer datatype 0 rectangle added successfully.")
    print(f"Output file: {gds_file}")

def get_all_layer_numbers(gds_file):
    """
    Retrieves all unique layer numbers present in the GDS file.

    Args:
        gds_file (str): Path to the input GDS file.

    Returns:
        list: A list of tuples, where each tuple contains (layer_number, datatype).
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(gds_file)

    # Create a set to store unique (layer, datatype) pairs
    layers = set()

    # Iterate over all layers in the layout
    for layer_index in layout.layer_indices():
        layer_info = layout.get_info(layer_index)
        layers.add((layer_info.layer, layer_info.datatype))

    # Convert the set to a sorted list
    sorted_layers = sorted(layers)

    return sorted_layers



def add_squares_to_layer(input_gds, output_gds, selected_layer, selected_datatype, square_size=5, spacing=10, keepout=5):
    """
    Adds squares to a specific layer in a GDS file.

    Parameters:
        input_gds (str): The path to the input GDS file.
        output_gds (str): The path to the output GDS file.
        selected_layer (int): The layer number to add squares to.
        selected_datatype (int): The datatype number to add squares to.
        square_size (int, optional): The size of the squares to be added. Defaults to 5.
        spacing (int, optional): The spacing between squares. Defaults to 10.
        keepout (int, optional): The keepout area around the existing shapes. Defaults to 5.
    """
    # Read the GDS file
    gdsii = gdspy.GdsLibrary(infile=input_gds)

    # Create a new GDS file for the output
    gdsii_new = gdspy.GdsLibrary()

    # Adjust spacing
    spacing = spacing*2

    # Copy all cells from the input GDS to the output GDS
    for cell_name, cell in gdsii.cells.items():
        gdsii_new.add(cell)

        # Select the layer and datatype
        layer_to_process = (selected_layer, selected_datatype)

        # Get the shapes in the selected layer
        polygons = cell.get_polygons(by_spec=True)

        if layer_to_process in polygons:
            # Create a new layer and datatype for the squares
            new_layer = selected_layer
            new_datatype = selected_datatype + 2

            # Calculate the bounding box of the existing shapes
            shapes = polygons[layer_to_process]
            for shape in shapes:
                shape = shape.tolist()  # Ensure shape is a list of tuples

                min_x = min(shape, key=lambda p: p[0])[0] + keepout
                max_x = max(shape, key=lambda p: p[0])[0] - keepout
                min_y = min(shape, key=lambda p: p[1])[1] + keepout
                max_y = max(shape, key=lambda p: p[1])[1] - keepout

                # Add squares within the geometry boundaries with keepout spacing
                x = min_x
                while x <= max_x:
                    y = min_y
                    while y <= max_y:
                        square = gdspy.Rectangle(
                            (x, y), (x + square_size, y + square_size),
                            layer=new_layer, datatype=new_datatype
                        )
                        # Check if the square points and the keepout area are within the existing polygon
                        square_points = [(x, y), (x + square_size, y), (x + square_size, y + square_size), (x, y + square_size)]
                        keepout_points = [
                            (x - keepout, y - keepout), (x + square_size + keepout, y - keepout),
                            (x + square_size + keepout, y + square_size + keepout), (x - keepout, y + square_size + keepout)
                        ]
                        if np.all(gdspy.inside(square_points, [shape])) and np.all(gdspy.inside(keepout_points, [shape])):
                            cell.add(square)
                        y += spacing
                    x += spacing

    # Write the new GDS file
    gdsii_new.write_gds(output_gds)

def create_cheesing_effect(input_gds, output_gds, selected_layer, selected_datatype):
    """
    Creates a cheesing effect on a GDS file by subtracting a square layer from an original layer.

    Parameters:
    - input_gds (str): The path to the input GDS file.
    - output_gds (str): The path to save the modified GDS file.
    - selected_layer (int): The layer number of the original and square layers.
    - selected_datatype (int): The datatype number of the original and square layers.

    Returns:
    None
    """
    # Load the GDS file
    layout = kdb.Layout()
    layout.read(input_gds)

    # Get the top cell
    top_cell = layout.top_cell()

    # Define the layers
    original_layer = layout.layer(selected_layer, selected_datatype)
    square_layer = layout.layer(selected_layer, selected_datatype + 1)

    # Boolean operation: A not B
    original_shapes = kdb.Region(top_cell.shapes(original_layer))
    square_shapes = kdb.Region(top_cell.shapes(square_layer))
    result_shapes = original_shapes - square_shapes

    # Clear original layer and add result shapes
    top_cell.shapes(original_layer).clear()
    top_cell.shapes(original_layer).insert(result_shapes)

    # Remove the squares layer
    top_cell.shapes(square_layer).clear()

    # Write the new GDS file
    layout.write(output_gds)

def bias_gds_features(input_gds, output_gds, bias, layer_number, datatype_number=None):
    """
    Biases features on a specific layer in a GDS file by expanding or contracting them by the specified amount using KLayout.

    Parameters:
        input_gds (str): The path to the input GDS file.
        output_gds (str): The path to the output GDS file.
        bias (float): The amount by which to bias the features (in microns). Positive values expand features, negative values contract them.
        layer_number (int): The layer number of the features to bias.
        datatype_number (int, optional): The datatype number of the features to bias. If None, all datatypes on the layer are biased.
    """
    # Load the GDS file using KLayout
    layout = kdb.Layout()
    layout.read(input_gds)

    # Convert bias from microns to database units (nanometers)
    bias_db = bias * 1000  # KLayout uses nanometers as the unit

    # Get the layer indices for the specified layer and datatype
    if datatype_number is None:
        # If datatype_number is None, get all datatypes on the layer
        layer_indices = []
        for layer_info in layout.layer_infos():
            if layer_info.layer == layer_number:
                layer_indices.append(layout.layer(layer_info.layer, layer_info.datatype))
    else:
        # Get the specific layer index
        layer_indices = [layout.layer(layer_number, datatype_number)]

    # Iterate through all cells
    for cell in layout.each_cell():
        # Iterate through the specified layers
        for layer_index in layer_indices:
            shapes = cell.shapes(layer_index)
            region = kdb.Region(shapes)
            if not region.is_empty():
                # Perform the bias (grow or shrink)
                region.size(bias_db)
                # Clear the original shapes
                shapes.clear()
                # Add the biased shapes back
                shapes.insert(region)

    # Write the biased GDS file
    layout.write(output_gds)

def convert_gds_to_dxf(input_gds:str, output_destination:str = None, output_name:str = None,  database_unit:float = 0.001, scale_factor:float = 0.001, polyline_mode:int = 0): 
    """
    Converts a given gds layout into dxf using Klayout

    Use tip: make sure to wrap strings and string literals with the 'r' prefix to help Python
        interpreate string as a filepath -- example: r"my/file/path" or r"my\file\path"

    Parameters: 
        input_gds (str, required): the path to gds file which will be converted 

        output_destination (str, optional): the path to the directory you wish to send gds to; 
            defaults to the same directory as the input if none is specified

        output_name (str, optional): the name of the final gds file you wish to produce; 
            defaults to the same name as the input_gds file

        database_unit (float, optional): set the database unit (ie: the unit the design is rendered at) in Klayout; default value is 0.001 micrometers

        scale_factor (float,optional): set scale factor for transforming a gds file into the dxf format; default value is 0.001

        polyline_mode (int,optional): set the polyline mode for a dxf file in Klayout - which determines how polyline and lwpolyline entities are converted treated
            in Klayout; Default value is 0. Possible values for this parameter are: 
            0 [automatic], 
            1 [keep's 0-width line geometry constant], 
            2 [create polygons from closed polulines with width =0], 
            3 [merge all lines with width = 0 into polygons], and
            4 [merge all lines with width = 0 into polygons plus auto-close open contours].
            If some geometry on a given input circuit is anomolous, consider selecting a different value for this parameter!
    """ 
    
    #Sanity Checks
    assert type(input_gds) == str, "Error: 'input_gds' given was not given as a str, please input path to gds file" 
    assert polyline_mode in [0,1,2,3,4], "Error: invalid polyline_mode given! Supported values are 0,1,2,3, and 4; see docstring for details"

    filename, extension = splitext(basename(input_gds))
    assert extension in ['.gds', ".gds"], "Error: input file is not a .gds file, it is %s" % (extension) 

    #Process output_destination: 
    if output_destination == None:
        output_destination = abspath(input_gds)
    #turn into a path-object so forward and backslashes are
    #handled properly -- even if the user errorneously puts in the wrong slash 
    elif output_destination != None:
        output_destination = abspath(output_destination)

    
    #Process output_name:
    if output_name == None:
        output_name =  filename + ".dxf"
    #need to check if the user specified an extension or not, and process appropriately  
    if output_name != None:
        fname, ext = splitext(output_name)
        if ext in ['.dxf', ".dxf"]:
            pass  
        else:
           output_name = fname + ".dxf" 

    final_output = join(output_destination, output_name)

    print("Input recieved: processing %s%s " % (filename,extension)) 

    #Use Klayout to load in GDS design
    layout = kdb.Layout()
    layout.read(input_gds)

    #Load in top layer(s) of GDS design to let user know to double check their design(s)
    top_cells = layout.top_cells()
    top_cell_names = [cell.name for cell in top_cells] 

    if len(top_cell_names) > 1:
        print("Current gds design has an ambigiously defined top layer -- might cause issues in conversion to DXF!")
    else:
        print("Current gds design is compatible -- performing conversion")


    #Conversion process:
    if len(top_cell_names) == 1:
        new_output_opts = kdb.SaveLayoutOptions()   
        new_output_opts.format = "DXF"
        new_output_opts.dbu = database_unit 
        new_output_opts.scale_factor = scale_factor 
        new_output_opts.dxf_polyline__mode = polyline_mode 
        layout.write(final_output, new_output_opts)
    #For multiple top layers
    else:
        for t_cell in top_cells:
            new_output_opts = kdb.SaveLayoutOptions()
            new_output_opts.format = "DXF"
            new_output_opts.dbu = database_unit 
            new_output_opts.scale_factor = scale_factor 
            new_output_opts.dxf_polyline_mode = polyline_mode 
            new_output_opts.select_cell( t_cell.cell_index() )
            tcdxf = "__%s__%s" % (t_cell.name, final_output)
            layout.write(tcdxf, new_output_opts)

    print("Processing complete: created new .DXF file")

def convert_dxf_to_gds(input_dxf:str, output_destination:str = None, output_name:str = None, database_unit:float = 0.001, dxf_unit:int = 1, polyline_mode:int = 0): 
    """
    Converts a given dxf layout into gds using Klayout

    Use tip: make sure to wrap strings and string literals with the 'r' prefix to help Python
        interpreate string as a filepath -- example: r"my/file/path" or r"my\file\path"

    Parameters: 
        input_dxf (str, required): the path to gds file which will be converted

        output_destination (str, optional): the path to the directory you wish to send gds to; 
            defaults to the same directory as the input if none is specified

        output_name (str, optional): the name of the final gds file you wish to produce, with proper extension; 
            defaults to the same name as the input_dxf file

        database_unit (float, optional): set the database unit (ie: the unit the design is rendered at) in Klayout; default value is 0.001 micrometers

        dxf_unit (int,optional): set the unit for value for a dxf file in Klayout; default value is 1 mircometer

        polyline_mode (int,optional): set the polyline mode for a dxf file in Klayout - which determines how polyline and lwpolyline entities are converted treated
            in Klayout; Default value is 0. Possible values for this parameter are: 
            0 [automatic], 
            1 [keep's 0-width line geometry constant], 
            2 [create polygons from closed polulines with width =0], 
            3 [merge all lines with width = 0 into polygons], and
            4 [merge all lines with width = 0 into polygons plus auto-close open contours].
            If some geometry on a given input circuit is anomolous, consider selecting a different value for this parameter! 
    """ 
    
    #Sanity Checks
    assert type(input_dxf) == str, "Error: 'input_dxf' given was not given as a str, please input path to DXF file" 
    assert polyline_mode in [0,1,2,3,4], "Error: invalid polyline_mode given! Supported values are 0,1,2,3, and 4; see docstring for details"

    #The empty string is included here, because on Linux Systems .dxf files are sometimes valid without a propper extension name
    filename, extension = splitext(basename(input_dxf))
    assert extension in ['.dxf', ".dxf", ""], "Error: input file is not a .dxf file, it is %s" % (extension) 

    #Process output_destination: 
    if output_destination == None:
        output_destination = abspath(input_dxf)
    #turn into a path-object so forward and backslashes are
    #handled properly -- even if the user errorneously puts in the wrong slash 
    elif output_destination != None:
        output_destination = abspath(output_destination)

    #Process output_name:
    if output_name == None:
        output_name =  filename + ".gds"
    #need to check if the user specified an extension or not, and process appropriately  
    if output_name != None:
        fname, ext = splitext(output_name)
        if ext in ['.gds', ".gds"]:
           pass 
        else:
           output_name = fname + ".gds" 

    final_output = join(output_destination, output_name)


    print("Input recieved: processing %s%s " % (filename,extension)) 

    #Use Klayout to load in DXF design
    layout = kdb.Layout()

    #Define layout options to properly load in DXF
    load_options = kdb.LoadLayoutOptions()
    load_options.dxf_dbu   = database_unit 
    load_options.dxf_unit  = dxf_unit 
    load_options.dxf_polyline_mode = polyline_mode 
    layout.read(input_dxf, load_options)

    #Conversion process:
    new_output_opts = kdb.SaveLayoutOptions()
    new_output_opts.format = "GDS2"
    layout.write(final_output, new_output_opts)


    print("Processing complete: created new .gds file")
