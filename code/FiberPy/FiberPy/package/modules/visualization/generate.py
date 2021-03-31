
# Get location of blender modules.
import os
import time
import json
import sys
import argparse  # to parse options for us and print a nice help message

# Blender packages that are located in Blender's Python Path.
import bpy
import bmesh
import mathutils

ROOT = os.path.realpath(os.path.dirname(__file__))
SARC_STRUCTS_ROOT = os.path.join(ROOT, "..", "util")
sys.path.append(ROOT)
sys.path.append(SARC_STRUCTS_ROOT)

import objects as b_obj
import sarc_structs


###################################################################################################
### Classes
###################################################################################################

class BlenderParams():
  """Container for holding parameters for Blender visualization"""
  def __init__(self, args, xscale=1/10., yzscale=44.):
    print ("I'd like to refactor this and have a flag that creates different dictionaries based on "
      "mechanical or biochemical rendering.")
    ## Create the containers to hold settings.
    # For creating filaments and filament features.
    self.radii = dict()
    self.lengths = {"lambda":40} # Where lambda represents the length of the bare zone.
    self.heights = {}
    self.angles = {}
    self.angles = {
      "tropomyosin_off_angle":0,
      "tropomyosin_on_angle":-45,
      "tropomyosin_random_angle":15,
      "cb_separation":10
    }
    # For creating M and Z lines.
    self.boxes = {
      "thickness":10, 
      "padding":10
    }
    # Extensions in x direction for cross bridge states.
    self.extensions = {1:0, 2:0, 3:5, 4:10, 5:0}
    # Blender doesn't handle large scales very well, so I'm adding in a scale parameter.
    self.xscale = xscale
    self.yzscale = xscale * yzscale

    # These are the biological settings.
    if args["render_mode"] == 0:
      raise RuntimeError("Biochemical rendering is not supported yet.")
      # self.radii["thin"] = 5.5 # 5.5 nm actin monomer, two monomers making up filament.

    # These are the mechanical settings.
    else:
      self.radii["thick"] = 7.5
      self.radii["thin"] = 4.5
      self.radii["titin"] = 1.5
      self.radii["actin_bs_biological"] = 2.25
      self.radii["actin"] = 2
      self.radii["tropomyosin"] = 2
      self.radii["cb_stalk"] = 1.75 
      self.radii["cb_lever_arm"] = 0.5
      self.radii["cb_binding_site"] = 1.5
      self.radii["cb_spring"] = 0.5
      self.radii["actin_bs_mechanical"] = 1      # ad hoc.
      self.radii["thick_spring"] = 0.1 * self.radii["thick"]
      self.radii["thin_spring"] = 0.1 * self.radii["thin"]

      self.radii["cb_spring"] = 0.2 * self.radii["cb_stalk"]
      
      self.lengths["actin_bs"] = 4
      self.lengths["tropomyosin"] = 5
      self.lengths["m_node"] = 3.55
      self.lengths["cb_spring"] = 8        # totally ad hoc
      self.lengths["cb_rand_x"] = 2 
      self.lengths["a_node"] = 2           # This should be 2 * actin_bs_mechanical if it's the 
                                           # mechanical rendering.
      self.lengths["thick_spring"] = 13.5  # Taken from the myosin resting length.
      self.lengths["thin_spring"] = 5.375  # Taken from the actin resting length.
      
      self.heights["cb_stalk"] = 5
      self.heights["cb_lever_arm"] = 3.5
      self.heights["cb_binding_site"] = 1
      self.heights["actin_bs_mechanical"] = 1.5

      self.num_spring_vertices = 12
      self.spring_steps = 32

    self.apply_scale(self.xscale)
  
  def apply_scale(self, scale):
    """Applies the scaling factor defined within the class to the parameter dictionaries."""
    for dictionary in self.get_dicts().values(): # It's absurd I don't use `vars()` here.
      for key, value in dictionary.items():
        dictionary[key] = value * scale
    return

  def get_dicts(self):
    """Returns a dict of the parameter dictionaries in this class."""
    return {key:value for key,value in self.__dict__.items() if isinstance(value,dict)}

###################################################################################################
### Convenience Functions
###################################################################################################

def create_object_lists_and_collection(half_sarcomere, b_params):
  """Returns the dictionary containing lists where objects will be stored and the collections
  
  Note: Collection is a Blender term and is how they manage large collections of objects 
    "efficiently."
  """
  obj_lists = dict()
 
  # Create thick and thin collections to hold the individual filament collections.
  thick = bpy.data.collections.new(name="Thick Collection")
  bpy.context.scene.collection.children.link(thick)
  thin = bpy.data.collections.new(name="Thin Collection")
  bpy.context.scene.collection.children.link(thin)

  for thick_obj in half_sarcomere.thick:
    thick_id = "m_" + str(thick_obj.id)
    obj_lists[thick_id] = []
    col = bpy.data.collections.new(thick_id)
    thick.children.link(col)

  for thin_obj in half_sarcomere.thin:
    thin_id = "a_" + str(thin_obj.id)
    obj_lists[thin_id] = []
    col = bpy.data.collections.new(thin_id)
    thin.children.link(col)
  
  # Create the collections for the primitives.
  prim = bpy.data.collections.new(name="Primitives")
  bpy.context.scene.collection.children.link(prim)
  # Set the primitive collection to invisible.
  prim.hide_viewport = True
  prim.hide_render = True

  return obj_lists

def link_all_objs(obj_lists):
  """Links all objects in `obj_lists` to the current collection
  
  Note: Collection is a Blender term and is how they manage large collections of objects 
    "efficiently."
  """
  start = time.time()

  for key, o_list in obj_lists.items():
    link = bpy.data.collections[key].objects.link
    for o in o_list:
      link(o)
  end = time.time()
  print ("Time for object linking to execute: {} seconds".format(end-start))

  return

def setup_blender_script(args):
  """Adds boilerplate blender scripting instructions."""
  # Set mode to object mode since we have to do most of the operations there.
  if bpy.context.active_object != None and bpy.context.active_object.mode != 'OBJECT':
      bpy.ops.object.mode_set(mode='OBJECT')

  # Delete objects if they currently exist.
  if bpy.data.objects != []:
      bpy.ops.object.select_all(action='SELECT')
      bpy.ops.object.delete()
  
  # Set the viewport background color to white.
  bpy.context.preferences.themes['Default'].view_3d.space.gradients.high_gradient.hsv = (0.0, 0.0, 
    1.0)
  bpy.context.scene.world.use_nodes = False
  bpy.context.scene.world.color = (1, 1, 1)

  setup_render_engine(args)

  return

def setup_render_engine(args):
  """Sets up the render engine that we want to use."""
  # Setting the render samples to 8 instead of 64 to cut render time.
  # Note: I don't think this affects much since the objects aren't very reflective.
  # bpy.context.scene.eevee.taa_render_samples = 8

  if args["render_quality"] == "high":
    bpy.context.scene.render.engine = 'CYCLES'
  elif args["render_quality"] == "medium":
    # Opting for Blender's workbench render engine because it does a good job and it's fast.
    bpy.context.scene.render.engine = 'BLENDER_WORKBENCH'
    bpy.context.scene.display.shading.light = 'STUDIO'
    bpy.context.scene.view_settings.view_transform = 'Filmic'
    bpy.context.scene.display.shading.show_object_outline = True
    bpy.context.scene.view_settings.exposure = 1.5
    bpy.context.scene.view_settings.gamma = 1.1
    bpy.context.scene.display.shading.show_specular_highlight = True
    bpy.context.scene.display.shading.shadow_intensity = 0.1

def read_and_scale_json_dump(json_dump_file_string, b_params):
  """Reads in and scales the parameters from the JSON dump file and returns result in dump_dict."""
  # Read in the JSON dump file.
  with open(json_dump_file_string, 'r') as f:
    dump_dict = json.load(f)
  
  half_sarcomere = sarc_structs.HalfSarcomere(dump_dict)

  half_sarcomere.hs_length *= b_params.xscale

  for thick_obj in half_sarcomere.thick:
    thick_obj.m_y *= b_params.yzscale
    thick_obj.m_z *= b_params.yzscale
    thick_obj.m_inter_crown_rest_length *= b_params.xscale
    thick_obj.m_lambda *= b_params.xscale
    for cb_i in range(thick_obj.m_no_of_cbs):
      thick_obj.cb_x[cb_i] *= b_params.xscale
  for thin_obj in half_sarcomere.thin:
    thin_obj.a_y *= b_params.yzscale
    thin_obj.a_z *= b_params.yzscale
    thin_obj.a_inter_bs_rest_length *= b_params.xscale
    for bs_i in range(thin_obj.a_no_of_bs):
      thin_obj.bs_x[bs_i] *= b_params.xscale
  
  # Get the cross bridge extensions.
  # extension_keys = [key for key in dump_dict["header"]["state_extensions"].keys()]
  extension_keys = [str(i) for i in range(10)]
  for key in extension_keys:
    # Get the state number.
    state = int(key)
    # b_params.extensions[state] = dump_dict["header"]["state_extensions"][key] * b_params.scale
    b_params.extensions[state] = 0

  return half_sarcomere

def render_screenshot(args, dump_file_index):
  """Renders a screen capture of the current geometry."""
  bpy.context.scene.render.filepath = os.path.join(args["output_file_directory"], 
    args["output_file_root"] + "_" + str(dump_file_index))
  bpy.context.scene.frame_set(1)
  bpy.ops.render.render(write_still=True)
  return

###################################################################################################
### Main Functions and Command Line Execution
###################################################################################################

def generate_model(args, output_file_name="generated_script"):
  """Reads FiberSim dump files and generates and saves a blender model of the simulation."""
  start = time.time()

  # Create an instance of the class holding the Blender parameters.
  b_params = BlenderParams(args)
  
  # Setup boilerplate blender code.
  setup_blender_script(args)

  # Get dictionary of materials in the script.
  materials = b_obj.create_materials(args)

  # Form object dictionary.
  object_dict = dict()

  # Read in the first JSON file and create an instance of the HalfSarcomere class.
  dump_file_template = args["dump_file_root"] + "{}.json"
  half_sarcomere = read_and_scale_json_dump(dump_file_template.format(args["dump_file_start"]),
    b_params)
  
  # Place the cameras.
  b_obj.create_cameras(args, half_sarcomere)

  # Place the lights.
  b_obj.create_lights(half_sarcomere)
  
  # Eventually replace with dump file parameter.
  num_half_sarcomeres = 1

  # Make the Z-lines.
  for i in range(num_half_sarcomeres):
    b_obj.create_z_line(half_sarcomere, b_params, object_dict, materials)

  # Create M lines.
  for i in range(num_half_sarcomeres):
    b_obj.create_m_line(half_sarcomere, b_params, object_dict, materials)
  
  # Create the lists and collections for the objects to be held in while creating geometry.
  obj_lists = create_object_lists_and_collection(half_sarcomere, b_params)

  # Create mesh primitives for binding sites, cross bridges, and nodes.
  b_primitives = b_obj.create_blender_primitives(b_params, args, materials)

  # Create the initial model from the first file.
  # Create the thick filaments.
  for i, thick_obj in enumerate(half_sarcomere.thick):
    fil_start_time = time.time()
    thick_key = "m_" + str(thick_obj.id)
    if thick_key not in args["filaments_to_hide"] and args["render_level"] in (0, 1):
      b_obj.create_thick_filament(half_sarcomere, i, b_params, args, materials, obj_lists, 
        b_primitives)
      print (thick_key, "time to draw:" ,time.time()-fil_start_time, 's')

  # Create the thin filaments.
  for i, thin_obj in enumerate(half_sarcomere.thin):
    fil_start_time = time.time()
    thin_key = "a_" + str(thin_obj.id)
    if thin_key not in args["filaments_to_hide"] and args["render_level"] in (0, 1):
      b_obj.create_thin_filament(half_sarcomere, i, b_params, args, materials, obj_lists, 
        b_primitives)
      print (thin_key, "time to draw:" ,time.time()-fil_start_time, 's')
  
  # Draw tropomyosin.
  # TODO
  # if args.render_level in (0,2):
  #   b_script = b_obj.create_tropomyosin(dump_dict, b_script, b_params)
  
  # Link all of the objects together.
  link_all_objs(obj_lists)

  # Take a screen capture of the first scene.
  if not args["no_render"]:
    render_screenshot(args, dump_file_index=args["dump_file_start"])
  
  print ("Took {} seconds to draw initial geometry.".format(time.time() - start))
  
  # Loop through the other files to render and adjust the model as needed.
  for i in range(int(args["dump_file_start"]) + 1, int(args["dump_file_end"])):
    update_model_from_json_dump(dump_file_template.format(i), b_params, args)

    # Take a screenshot of the updated model.
    if not args["no_render"]:
      render_screenshot(args, dump_file_index=i)
  
  # Add in the print statement for how long it took.
  print ("Took {} seconds to create entire animation.".format(time.time() - start))

  return

def update_model_from_json_dump(dump_file_path, b_params, args):
  """Updates an existing Blender model of half-sarcomere(s) based on JSON dump file(s)."""
  # Read the JSON file to form the half_sarcomere object.
  half_sarcomere = read_and_scale_json_dump(dump_file_path, b_params)
  
  # Loop through the filaments and update them.
  draw_start = time.time()
  for i, thick_obj in enumerate(half_sarcomere.thick):
    fil_start_time = time.time()
    thick_key = "m_" + str(thick_obj.id)
    if thick_key not in args["filaments_to_hide"] and args["render_level"] in (0, 1):
      b_obj.update_thick_filament(half_sarcomere, i, b_params, args)
      print (thick_key, "time to draw:" ,time.time()-fil_start_time, 's')
  
  for i, thin_obj in enumerate(half_sarcomere.thin):
    fil_start = time.time()
    thin_key = "a_" + str(thin_obj.id)
    if thin_key not in args["filaments_to_hide"] and args["render_level"] in (0, 1):
      b_obj.update_thin_filament(half_sarcomere, i, b_params, args)

  print ("Took {} seconds to update entire model.".format(time.time() - draw_start))

  return

def setup_parser():
  """Sets up and returns the argument parser."""
  # When --help or no args are given, print this help
  usage_text = (
    "Run blender in background mode with this script:"
    "  blender --background --python " + __file__ + " -- [options]"
  )

  parser = argparse.ArgumentParser(description=usage_text)

  parser.add_argument(
    "-j",
    "--json-file",
    help="The path to the JSON instruction file."
  )

  return parser

def recursive_param_check(d, type_dict, required_check=False, key=None):
  """Checks that the d[key] is of type type_dict[key] recursively.

  See check_required_params for a usage example.
  """
  if key == None or isinstance(type_dict[key], dict):
    if key != None:
      d = d[key]
      type_dict = type_dict[key]

    # Traverse deeper into the dictionary tree.
    for sub_key in type_dict.keys():
      recursive_param_check(d, type_dict, required_check, sub_key)

  elif not isinstance(type_dict[key], dict):
    if required_check:
      assert (key in d.keys()), "The parameter \""+key+"\" is not specified and it is required."
    # Check that the value in d matches with type in type_dict.
    assert isinstance(d[key], type_dict[key]), "The parameter \""+key+"\" is not of the right type."

def check_required_params(instruction_dict):
  """Checks that the required parameters are present in the instruction dictionary."""
  # Set up the required params type dictionary.
  # This is just a dictionary that is the exact same structure as the parameter dictionary, but
  # instead of having values in the place they would be when we read them in, we have a tuple of 
  # the possible types the values can have. Not to toot my own horn, but I think this is pretty 
  # clever.
  required_params = {
    "dump_file_root": str,
    "dump_file_start": int,
    "dump_file_end": int,
    "camera": {
      "location": {
        "x": (float, int),
        "y": (float, int),
        "z": (float, int)
      },
      "rotation": {
        "x": (float, int),
        "y": (float, int),
        "z": (float, int)
      }
    }
  }

  # Recursively check the parameters.
  recursive_param_check(instruction_dict, required_params, required_check=True)

def check_optional_params(instruction_dict):
  """Checks optional parameters in instruction_dict.

  Virtually the same as check_required_params but it's only checking the optional parameters.
  """
  optional_params = {
    "output_file_directory": str,
    "output_file_root": str,
    "render_level": int,
    "render_mode": int,
    "no_render": bool,
    "filaments_to_hide": list,
    "draw_mirrored_cb_connections": bool
  }

  recursive_param_check(instruction_dict, optional_params)

def dfs_set_param(dict_1, dict_default, key=None):
  """Depth-first default parameter setting for use with default parameter dictionaries."""
  if key == None or isinstance(dict_default[key], dict):
    if key != None:
      dict_1 = dict_1[key]
      dict_default = dict_default[key]

    # Recursively call the depth-first parameter setting.
    for sub_key in dict_default.keys():
      dfs_set_param(dict_1, dict_default, sub_key)
  
  elif key not in dict_1.keys():
    # Set the parameter to the default value if it isn't present in dict_1.
    dict_1[key] = dict_default[key]

def default_visualization_param_dict():
  """Returns a dictionary that contains the default visualization parameters."""
  default_dict = {
    "output_file_directory": os.path.join(ROOT, "blender_renders"),
    "output_file_root": "model",
    "render_level": 0,
    "render_mode": 1,
    "no_render": False,
    "render_quality": "medium",
    "filaments_to_hide": [],
    "draw_mirrored_cb_connections": False
  }

  return default_dict

def get_parameters_from_json_file(argv):
  """Gets the visualization parameters from a JSON file."""
  # Setup default parameter dictionary.
  default_dict = default_visualization_param_dict()

  # Get the argument parser.
  parser = setup_parser()

  # Grab only the arguments after the "--" flag.
  argv = argv[argv.index("--") + 1:]

  # Parse the arguments and get the dictionary form.
  args = vars(parser.parse_args(argv))

  # Read in the instruction file parameters.
  with open(args["json_file"], 'r') as f:
    instruction_dict = json.load(f)
  
  # Do some error checking to make sure the non-optional parameters are set and the right type.
  check_required_params(instruction_dict)

  # Set the default parameters if none are specified in the JSON instruction dict.
  dfs_set_param(instruction_dict, default_dict)

  # Check that the optional parameters are of the right type.
  check_optional_params(instruction_dict)

  # Set the paths as absolute paths so nothing unexpected happens with relative referencing.
  instruction_dict["output_file_directory"] = os.path.realpath(
    instruction_dict["output_file_directory"])

  return instruction_dict

def main():
  import sys       # to get command line args

  # Get the parameters from the JSON file.
  args = get_parameters_from_json_file(sys.argv)
  
  generate_model(args)

if __name__ == "__main__":
    main()