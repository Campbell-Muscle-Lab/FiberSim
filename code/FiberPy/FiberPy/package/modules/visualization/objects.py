# These are essentially useless since this is always used as a module, but it stops my linting 
# from blowing up so I'm keeping it.
import math

import bpy
import bmesh
import mathutils
import numpy as np


###################################################################################################
### Convenience Functions
###################################################################################################

def give_circle_vertices_and_edges(radius, num_vertices, offset=(0,0,0)):
    """Returns a list of 3D vertices on a circle aligned with the xz plane."""
    angle_separation = 2 * np.pi / num_vertices
    vertices = []
    for i in range(num_vertices):
        angle = angle_separation * i
        x = np.cos(angle) * radius
        z = np.sin(angle) * radius
        vertices.append((x+offset[0], 0+offset[1], z+offset[2]))
    edges = []
    for j in range(1,num_vertices):
        edges.append((j,j-1))
    edges.append((j,0))
    return vertices, edges

def get_rotation_angle_and_vector(vect, vect_orig=np.asarray([0,0,1])):
  """Returns the rotation angle and vector of rotation between vect and vect_orig."""
  assert isinstance(vect,np.ndarray), "vect must be an ndarray!!"
  assert vect.shape[0] == 3, "vect must be of length = 3!!"
  # Calculate the vector of rotation via cross product.
  v_of_rot = np.cross(vect, vect_orig)
  # Calculate the angle of rotation.
  ang = np.arccos(np.dot(vect,vect_orig) / (np.linalg.norm(vect) * np.linalg.norm(vect_orig)))

  return (ang, v_of_rot)

def position_scale_orient_cylinder(cyl_obj, base_vec, end_vec):
  """Positions, scales, and orients a cylinder primitive s.t. it spans base_vec to end_vec."""
  # Compute the arm rotation angle from the two coordinates using the law of cosines.
  b = end_vec - base_vec
  # NOTE: Following assumes the cylinder primitive is aligned with the Z-axis originally.
  rot_angle, rot_axis = get_rotation_angle_and_vector(b)

  # Find the cylinder location (in the middle of the two positions).
  loc = 0.5 * (base_vec + end_vec)
  scale = np.sqrt(np.sum(np.power(base_vec - end_vec, 2)))

  cyl_obj.location.x = loc[0]
  cyl_obj.location.y = loc[1]
  cyl_obj.location.z = loc[2]
  cyl_obj.scale[2] = scale
  cyl_obj.rotation_mode = 'AXIS_ANGLE'
  cyl_obj.rotation_axis_angle[0] = -rot_angle
  cyl_obj.rotation_axis_angle[1] = rot_axis[0]
  cyl_obj.rotation_axis_angle[2] = rot_axis[1]
  cyl_obj.rotation_axis_angle[3] = rot_axis[2]
  return

def find_adjacent_thin_fils_to_thick_fils(half_sarcomere, thick_index):
  """Returns list of filament IDs that are adjacent to the filament ID specified.

  Inputs:
    half_sarcomere -> Instance of util.sarc_structs.HalfSarcomere object. The half_sarcomere 
                 latttice object that is being drawn.
    thick_index -> int. Index of the thick filament you would like to find the adjacent actin 
                   filaments of.
  """
  # Pull out the relevant thick filament coordinates.
  thick_obj = half_sarcomere.thick[thick_index]
  thick_y = thick_obj.m_y
  thick_z = thick_obj.m_z

  # Pull out the coordinates for all of the thin filaments.
  thin_coords = [(thin_obj.a_y, thin_obj.a_z) for thin_obj in half_sarcomere.thin]

  # Compute the distances between the thick filament and all thin filaments.
  distances = [((thick_y - t_c[0])**2 + (thick_z - t_c[1])**2)**0.5 for t_c in thin_coords]

  # Find the minimum distance and find all thin filaments that fall within the distance threshold.
  error_margin = 1.05
  distance_threshold = np.sqrt(1.25) * error_margin * min(distances)
  thin_fil_indexes_sub_thresh = [i for i, dist in enumerate(distances) if dist < distance_threshold]
  # thin_fil_IDs = [half_sarcomere.thin[i].id for i in thin_fil_indexes_sub_thresh]
      
  # Find all filament ID's with distance below distance cutoff.
  return thin_fil_indexes_sub_thresh

def return_all_filament_y_coords(half_sarcomere):
  """Returns a list of all filament y coordinates"""
  m_ys = [thick_obj.m_y for thick_obj in half_sarcomere.thick]
  a_ys = [thin_obj.a_y for thin_obj in half_sarcomere.thin]
  return m_ys + a_ys

def return_all_filament_z_coords(half_sarcomere):
  """Returns a list of all filament z coordinates"""
  m_zs = [thick_obj.m_z for thick_obj in half_sarcomere.thick]
  a_zs = [thin_obj.a_z for thin_obj in half_sarcomere.thin]
  return m_zs + a_zs

def find_bounding_box_yz_positions(half_sarcomere):
  """Returns min/max yz positions of the bounding box of the geometry."""
  filament_y_positions = return_all_filament_y_coords(half_sarcomere)
  filament_z_positions = return_all_filament_z_coords(half_sarcomere)
  min_y = np.min(filament_y_positions)
  max_y = np.max(filament_y_positions)
  min_z = np.min(filament_z_positions)
  max_z = np.max(filament_z_positions)
  return ((min_y,max_y),(min_z,max_z))

def set_active_object_shade_smooth(obj, auto_smooth_angle = 45):
  """Sets the object to active to use smooth shading with autosmoothe which is a magic bullet."""
  obj.select_set(True)
  bpy.ops.object.shade_smooth()
  bpy.context.active_object.data.use_auto_smooth = True
  bpy.context.active_object.data.auto_smooth_angle= auto_smooth_angle * np.pi / 180
  return

def give_thick_key(thick_obj):
  return "m_" + str(thick_obj.id)

def give_thin_key(thin_obj):
  return "a_" + str(thin_obj.id)

###################################################################################################
### Functions for Creating Individual Primitives
###################################################################################################

def create_cameras(args, half_sarcomere):
  """Creates cameras in positions for rednering."""
  # Add a camera that will look down the length of the filaments.
  if "camera" in args.keys():
    print("ken_camera")
    print(args)
    
    rot = (
      args["camera"]["rotation"]["x"] * np.pi / 180,
      args["camera"]["rotation"]["y"] * np.pi / 180,
      args["camera"]["rotation"]["z"] * np.pi / 180
    )
    bpy.ops.object.camera_add(
      enter_editmode=False,
      align='VIEW',
      location=(args["camera"]["location"]["x"],
                args["camera"]["location"]["y"], 
                args["camera"]["location"]["z"]),
      rotation=rot
    )
    bpy.context.object.data.lens = 35

  else:
    bpy.ops.object.camera_add(
      enter_editmode=False,
      align='VIEW',
      #location=(0.6*dump_dict["header"]["hsl"], 0, 75),
      location=(103.22, -11.144, 1.7829),
      rotation=(np.pi/2, 0, np.pi/2)
    )
  bpy.context.scene.camera = bpy.context.object
  return

def create_lights(half_sarcomere):
  """Creates lights to brighten the rendered scene."""
  print ("Disabling lights in lieu of global lighting.")
  return

  # Space separating the lights and the objects.
  padding = 11
  plane_padding = 10
  # light_strength = 5000
  light_strength = 60

  # Find the x,y,z positions of the lights (centered at the half sarcomere).
  x_pos = half_sarcomere.hs_length / 2
  (min_y, max_y), (min_z, max_z) = find_bounding_box_yz_positions(half_sarcomere)
  light_width = 1.5 * (abs(min_y) + abs(max_y))
  light_height = 1.5 * half_sarcomere.hs_length
  y_z_positions = (
    (0.5 * (max_y + min_y), max_z + padding),
    (min_y - padding,       0.5 * (max_z + min_z)),
    (0.5 * (max_y + min_y), min_z - padding),
    (max_y + padding,       0.5 * (max_z + min_z))
  )

  # Add one light for each side of the half sarcomere.
  for i, (y_pos, z_pos) in enumerate(y_z_positions):
    rotation = (i * 90 * np.pi / 180, 0, 0)
    bpy.ops.object.light_add(
      type='AREA', 
      location=(x_pos, y_pos, z_pos),
      rotation=rotation
    )
    bpy.context.object.data.shape = 'RECTANGLE'
    bpy.context.object.data.size = light_height
    bpy.context.object.data.size_y = light_width
    bpy.context.object.data.energy = light_strength
    bpy.context.object.data.use_shadow = False
    bpy.context.object.hide_viewport = True

  # Need to add the end cap lights for the muscle fiber too. With lights pointing in both directions.
  y_pos = (max_y + min_y) / 2
  z_pos = (max_z + min_z) / 2
  # At start of half-sarcomere point towards muscle.
  bpy.ops.object.light_add(
    type = 'AREA',
    location=(-padding, y_pos, z_pos), 
    rotation=(0, -0.5 * np.pi, 0)
  )
  bpy.context.object.data.shape = 'RECTANGLE'
  bpy.context.object.data.size = light_height
  bpy.context.object.data.size_y = light_width
  bpy.context.object.data.energy = light_strength
  bpy.context.object.data.use_shadow = False
  bpy.context.object.hide_viewport = True
  # At start of half-sarcomere pointing away from muscle.
  bpy.ops.object.light_add(
    type = 'AREA',
    location=(1, y_pos, z_pos), 
    rotation=(0, 0.5 * np.pi, 0)
  )
  bpy.context.object.data.shape = 'RECTANGLE'
  bpy.context.object.data.size = max_z - min_z + 2 * plane_padding
  bpy.context.object.data.size_y = max_y - min_y + 2 * plane_padding
  bpy.context.object.data.energy = light_strength
  bpy.context.object.data.use_shadow = False
  bpy.context.object.hide_viewport = True
  # At end of half-sarcomere pointing towards muscle.
  bpy.ops.object.light_add(
    type = 'AREA',
    location=(half_sarcomere.hs_length + padding, y_pos, z_pos), 
    rotation=(0, 0.5 * np.pi, 0)
  )
  bpy.context.object.data.shape = 'RECTANGLE'
  bpy.context.object.data.size = light_height
  bpy.context.object.data.size_y = light_width
  bpy.context.object.data.energy = light_strength
  bpy.context.object.data.use_shadow = False
  bpy.context.object.hide_viewport = True
  # At end of half-sarcomere pointing away from muscle.
  bpy.ops.object.light_add(
    type = 'AREA',
    location=(half_sarcomere.hs_length - 1, y_pos, z_pos), 
    rotation=(0, -0.5 * np.pi, 0)
  )
  bpy.context.object.data.shape = 'RECTANGLE'
  bpy.context.object.data.size = max_z - min_z + 2 * plane_padding
  bpy.context.object.data.size_y = max_y - min_y + 2 * plane_padding
  bpy.context.object.data.energy = light_strength
  bpy.context.object.data.use_shadow = False
  bpy.context.object.hide_viewport = True

  return

def create_planes(dump_dict, materials_dictionary):
  """Creates and places the planes that are necessary for the render to have a white background."""
  # Space separating the planes and the objects
  padding = 10

  # Find the x,y,z coordinates of the planes.
  x_pos = dump_dict["header"]["hsl"] / 2
  (min_y, max_y), (min_z, max_z) = find_bounding_box_yz_positions(dump_dict)
  y_z_positions = (
    (0.5 * (max_y + min_y), max_z + padding),
    (min_y - padding,       0.5 * (max_z + min_z)),
    (0.5 * (max_y + min_y), min_z - padding),
    (max_y + padding,       0.5 * (max_z + min_z))
  )

  for i, (y_pos,z_pos) in enumerate(y_z_positions):
    rotation = (i * 90 * np.pi / 180, 0, 0)
    bpy.ops.mesh.primitive_plane_add(
      size=1.5 * dump_dict["header"]["hsl"], 
      enter_editmode=False, 
      location=(x_pos, y_pos, z_pos), 
      rotation=rotation
    )
    bpy.context.active_object.data.materials.append(materials_dictionary["background_plane"])
    bpy.context.object.hide_viewport = True

  # Need to add the end caps for the muscle fiber too.
  end_plane_size = max(((max_y - min_y), (max_z - min_z))) + 2 * padding
  y_pos = (max_y + min_y) / 2
  z_pos = (max_z + min_z) / 2
  bpy.ops.mesh.primitive_plane_add(
    size=1.25 * end_plane_size, 
    enter_editmode=False, 
    location=(-padding, y_pos, z_pos), 
    rotation=(0, 0.5 * np.pi, 0)
  )
  bpy.context.active_object.data.materials.append(materials_dictionary["background_plane"])
  bpy.context.object.hide_viewport = True
  
  bpy.ops.mesh.primitive_plane_add(
    size=1.25 * end_plane_size, 
    enter_editmode=False, 
    location=(dump_dict["header"]["hsl"] + padding, y_pos, z_pos), 
    rotation=(0, 0.5 * np.pi, 0)
  )
  bpy.context.active_object.data.materials.append(materials_dictionary["background_plane"])
  bpy.context.object.hide_viewport = True

  return

def create_basic_myosin_head_primitives(b_params, args, b_primitives, materials_dictionary):  
  """Create a myosin head mesh primitives."""
  m_head_primitives = dict()

  # Add the stalk that is the same morphology across all cross bridge states.
  bpy.ops.mesh.primitive_cylinder_add(
    radius=b_params.radii["cb_stalk"],
    depth=b_params.heights["cb_stalk"], 
    enter_editmode=False, 
    location=(0, 0, 0), 
    vertices=20
  )
  # stalk = bpy.ops.object
  stalk_obj = bpy.context.active_object

  # Add a sphere onto the stalk to make the spring seem like it's a part of the stalk
  bpy.ops.mesh.primitive_ico_sphere_add(
    radius=b_params.radii["cb_stalk"],
    subdivisions=3,
    enter_editmode=False,
    location=(0, 0, 0.5 * b_params.heights["cb_stalk"])
  )
  sphere = bpy.context.active_object
  sphere.name = "cb_stalk"
  m_head_primitives["stalk"] = bpy.context.active_object
  stalk_obj.select_set(True)
  bpy.ops.object.join()

  # Set the custom shading options for the arm to make it look a whole lot better.
  set_active_object_shade_smooth(m_head_primitives["stalk"])
  change_object_collection(m_head_primitives["stalk"], bpy.data.collections["Primitives"])

  # Add the lever arm for each cross bridge state.
  bpy.ops.mesh.primitive_cylinder_add(
    radius=b_params.radii["cb_lever_arm"],
    depth=1,
    enter_editmode=False, 
    location=(0, 0, 0),
    vertices=20
  )
  # lever_arm = bpy.ops.object
  m_head_primitives["lever_arm"] = bpy.context.active_object
  bpy.ops.object.editmode_toggle()
  me = m_head_primitives["lever_arm"].data
  bm = bmesh.from_edit_mesh(me)

  # Have to refresh the lookup table for some reason.
  bm.faces.ensure_lookup_table()

  # Scale the bottom face to make a 'flask' shape.
  scale_factor = 1.75
  for f in bm.faces:
      if all([v.co[2] < 0 - 1e-2 for v in f.verts]):
          for v in f.verts:
              v.co[0] *= scale_factor
              v.co[1] *= scale_factor
  me.update()
  bm.free()
  bpy.ops.object.editmode_toggle()

  # Create sphere at the bottom of the lever arm to round it off.
  bpy.ops.mesh.primitive_ico_sphere_add(
    radius=scale_factor * b_params.radii["cb_lever_arm"],
    subdivisions=2,
    enter_editmode=False,
    location=(0, 0, -0.5)
  )
  sph = bpy.context.active_object
  m_head_primitives["lever_arm"].select_set(True)
  sph.select_set(True)
  bpy.context.view_layer.objects.active = m_head_primitives["lever_arm"]
  bpy.ops.object.join()

  set_active_object_shade_smooth(m_head_primitives["lever_arm"])
  change_object_collection(m_head_primitives["lever_arm"], bpy.data.collections["Primitives"])

  create_m_head_spring(b_params, m_head_primitives)
  create_m_head_catalytic_domain(b_params, m_head_primitives)
  copy_materials_to_m_heads(b_primitives, m_head_primitives, materials_dictionary)

  return

def create_m_head_spring(b_params, m_head_primitives):
  """Creates the spring that connects the lever arm and the stalk of the myosin head."""
  # Make the spring connecting the stalk and lever arm.
  spring_data = bpy.data.meshes.new("stalk_connecting_spring")
  line_data = bpy.data.meshes.new("line")

  # Get the vertexes and edges from external function in generate_blender_script.py
  verts,edges = give_circle_vertices_and_edges(
    radius=b_params.radii["cb_spring"],
    num_vertices=b_params.num_spring_vertices,
    offset=(-2 * b_params.radii["cb_spring"], 0, 0)
  )
  face = [tuple([i for i in range(len(verts))])]

  # height vector.
  line_verts = [(0,0,1), (0,0,-1)]
  line_edge = [(0,1)]

  spring_data.from_pydata(verts, edges, face)
  spring_data.update()
  line_data.from_pydata(line_verts, line_edge, [])
  line_data.update()

  # Create the object
  spring = bpy.data.objects.new("stalk_connecting_spring", spring_data)
  line_obj = bpy.data.objects.new("line", line_data)
  bpy.data.collections[0].objects.link(spring)
  bpy.data.collections[0].objects.link(line_obj)

  # Add the screw modifier
  spring.select_set(True)
  spring.modifiers.new("stalk_connecting_spring",'SCREW')
  spring.modifiers["stalk_connecting_spring"].object = bpy.data.objects["line"]
  spring.modifiers["stalk_connecting_spring"].screw_offset = 5.5 * b_params.radii["cb_spring"]
  spring.modifiers["stalk_connecting_spring"].iterations = b_params.lengths["cb_spring"] / 5.5 /  b_params.radii["cb_spring"]
  spring.select_set(True)

  bpy.context.view_layer.objects.active = spring
  bpy.ops.object.modifier_apply(modifier="stalk_connecting_spring")
  spring.rotation_euler[1] = 3.14159 / 2
  spring.location.z = 0.5 * b_params.heights["cb_stalk"] + 0.25 * b_params.radii["cb_stalk"]

  # Delete the unnecessary line.
  bpy.ops.object.select_all(action="DESELECT")
  line_obj.select_set(True)
  bpy.ops.object.delete()

  m_head_primitives["spring"] = spring
  set_active_object_shade_smooth(m_head_primitives["spring"])
  change_object_collection(m_head_primitives["spring"], bpy.data.collections["Primitives"])

  # Store the original length for later scaling
  b_params.lengths["original_cb_spring"] = spring.dimensions[2] - b_params.radii["cb_spring"]

  return

def create_m_head_catalytic_domain(b_params, m_head_primitives):
  """Creates the catalytic domain of the myosin head."""
  # Add the catalytic domain.
  bpy.ops.mesh.primitive_cylinder_add(
    radius=b_params.radii["cb_lever_arm"],
    depth=b_params.heights["cb_binding_site"],
    enter_editmode=True, 
    location=(0,0,0),
    vertices=15
  )
  catalytic_domain = bpy.context.active_object
  catalytic_domain.name = "catalytic_domain"

  me = catalytic_domain.data
  bm = bmesh.from_edit_mesh(me)

  # Have to refresh the lookup table for some reason.
  bm.faces.ensure_lookup_table()

  # Scale the top face to make a 'flask' shape.
  scale_factor = b_params.radii["actin_bs_mechanical"] / b_params.radii["cb_lever_arm"]
  for f in bm.faces:
      if all([v.co[2] > 0 + 1e-2 for v in f.verts]):
          for v in f.verts:
              v.co[0] *= scale_factor
              v.co[1] *= scale_factor

  # Now we set the origin of the binding site to the bottom of it so that rotation
  #   rotates it while remaining in contact with the lever arm.
  for f in bm.faces:
      if all([v.co[2] < 0 - 1e-2 for v in f.verts]):
          o = f.calc_center_median()
          bmesh.ops.translate(bm,
                verts = bm.verts,
                vec = -o,
                )
          bmesh.update_edit_mesh(me)
          me.update()           
          # move the object globally to reflect
          mw = catalytic_domain.matrix_world
          t = mw @ o - mw @ mathutils.Vector()
          mw.translation += t

  me.update()
  bm.free()

  bpy.ops.object.editmode_toggle()

  bpy.ops.mesh.primitive_ico_sphere_add(
    radius=b_params.radii["cb_lever_arm"],
    enter_editmode=False,
    location=(0, 0, -0.5 * b_params.heights["cb_binding_site"]),
    subdivisions=2
  )
  sph = bpy.context.active_object

  # Join the objects together.
  catalytic_domain.select_set(True)
  sph.select_set(True)
  bpy.context.view_layer.objects.active = catalytic_domain
  bpy.ops.object.join()

  # Shade with autosmooth
  set_active_object_shade_smooth(catalytic_domain)

  # Deselect everything.
  bpy.ops.object.select_all(action='DESELECT')

  m_head_primitives["catalytic_domain"] = catalytic_domain
  change_object_collection(m_head_primitives["catalytic_domain"], bpy.data.collections["Primitives"])

  return

def copy_materials_to_m_heads(b_primitives, m_head_primitives, materials):
  """Creates copies of each myosin head primitive for each cross bridge state available."""
  num_states = 5
  for i in range(num_states+1):
    cb_k = "cb_state_"+str(i)
    b_primitives[cb_k] = dict()

    b_primitives[cb_k]["stalk"] = m_head_primitives["stalk"].copy()
    b_primitives[cb_k]["stalk"].data = m_head_primitives["stalk"].data.copy()
    b_primitives[cb_k]["stalk"].data.materials.append(materials[cb_k])
    b_primitives[cb_k]["stalk"].name = "stalk_" + cb_k
    bpy.data.collections["Primitives"].objects.link(b_primitives[cb_k]["stalk"])

    b_primitives[cb_k]["catalytic_domain"] = m_head_primitives["catalytic_domain"].copy()
    b_primitives[cb_k]["catalytic_domain"].data = m_head_primitives["catalytic_domain"].data.copy()
    b_primitives[cb_k]["catalytic_domain"].data.materials.append(materials[cb_k])
    b_primitives[cb_k]["catalytic_domain"].name = "cat_" + cb_k
    bpy.data.collections["Primitives"].objects.link(b_primitives[cb_k]["catalytic_domain"])

    b_primitives[cb_k]["lever_arm"] = m_head_primitives["lever_arm"].copy()
    b_primitives[cb_k]["lever_arm"].data = m_head_primitives["lever_arm"].data.copy()
    b_primitives[cb_k]["lever_arm"].data.materials.append(materials[cb_k])
    b_primitives[cb_k]["lever_arm"].name = "arm_" + cb_k
    bpy.data.collections["Primitives"].objects.link(b_primitives[cb_k]["lever_arm"])
    
    b_primitives[cb_k]["spring"] = m_head_primitives["spring"].copy()
    b_primitives[cb_k]["spring"].data = m_head_primitives["spring"].data.copy()
    b_primitives[cb_k]["spring"].data.materials.append(materials[cb_k])
    b_primitives[cb_k]["spring"].name = "spring_" + cb_k
    bpy.data.collections["Primitives"].objects.link(b_primitives[cb_k]["spring"])

  return

def create_thick_spring_primitive(b_params, b_primitives, materials):
  """Appends instructions to create a spring primitive to b_script."""
  spring_core_gap = 0.6 * b_params.radii["thick"]
  offset = (-(0.5 * spring_core_gap + b_params.radii["thick_spring"]),0,0)
  spring_height = 5.5 * b_params.radii["thick_spring"]

  spring_data = bpy.data.meshes.new("node_connecting_spring")
  line_data = bpy.data.meshes.new("line")

  verts,edges = give_circle_vertices_and_edges(
    radius=b_params.radii["thick_spring"], 
    num_vertices=b_params.num_spring_vertices, 
    offset=offset
  )

  # height vector.
  line_verts = [(0,0,0.5), (0,0,-0.5)]
  line_edge = [(0,1)]

  spring_data.from_pydata(verts, edges, [])
  spring_data.update()
  line_data.from_pydata(line_verts, line_edge, [])
  line_data.update()

  # Create the object
  obj = bpy.data.objects.new("thick_spring", spring_data)
  line_obj = bpy.data.objects.new("line", line_data)
  bpy.data.collections[0].objects.link(obj)
  bpy.data.collections[0].objects.link(line_obj)

  # Add the screw modifier
  obj.select_set(True)

  obj.modifiers.new("Spring",'SCREW')
  obj.modifiers["Spring"].object = bpy.data.objects["line"]
  obj.modifiers["Spring"].screw_offset = spring_height
  obj.modifiers["Spring"].steps = b_params.spring_steps
  obj.modifiers["Spring"].iterations = b_params.lengths["thick_spring"] / spring_height
  obj.select_set(True)
  bpy.context.view_layer.objects.active = obj
  bpy.ops.object.modifier_apply(modifier="Spring")
  obj.data.materials.append(materials["thick_spring"])
  obj.rotation_euler[1] = 3.14159 / 2
  obj.scale[0] = 1.5
  obj.scale[1] = 1.5

  bpy.ops.object.select_all(action="DESELECT")
  line_obj.select_set(True)
  bpy.ops.object.delete()

  b_primitives["thick_spring"] = obj
  set_active_object_shade_smooth(obj)
  change_object_collection(obj, bpy.data.collections["Primitives"])
  return

def create_thin_spring_primitive(b_params, b_primitives, materials):
  """Appends instructions to create a spring primitive to b_script."""
  spring_core_gap = 0.6 * b_params.radii["thin"]
  offset = (-(0.5 * spring_core_gap + b_params.radii["thin_spring"]),0,0)
  spring_height = 7.0 * b_params.radii["thin_spring"]

  spring_data = bpy.data.meshes.new("node_connecting_spring")
  line_data = bpy.data.meshes.new("line")

  verts,edges = give_circle_vertices_and_edges(
    radius=b_params.radii["thin_spring"], 
    num_vertices=b_params.num_spring_vertices, 
    offset=offset
  )

  # height vector.
  line_verts = [(0,0,0.5), (0,0,-0.5)]
  line_edge = [(0,1)]

  spring_data.from_pydata(verts, edges, [])
  spring_data.update()
  line_data.from_pydata(line_verts, line_edge, [])
  line_data.update()

  # Create the object
  obj = bpy.data.objects.new("thin_spring", spring_data)
  line_obj = bpy.data.objects.new("line", line_data)
  bpy.data.collections[0].objects.link(obj)
  bpy.data.collections[0].objects.link(line_obj)

  # Add the screw modifier
  obj.select_set(True)

  obj.modifiers.new("Spring",'SCREW')
  obj.modifiers["Spring"].object = bpy.data.objects["line"]
  obj.modifiers["Spring"].screw_offset = spring_height
  obj.modifiers["Spring"].steps = b_params.spring_steps
  obj.modifiers["Spring"].iterations = b_params.lengths["thin_spring"] / spring_height
  obj.select_set(True)
  bpy.context.view_layer.objects.active = obj
  bpy.ops.object.modifier_apply(modifier="Spring")
  obj.data.materials.append(materials["thin_spring"])
  obj.rotation_euler[1] = 3.14159 / 2
  obj.scale[0] = 1.5
  obj.scale[1] = 1.5

  bpy.ops.object.select_all(action="DESELECT")
  line_obj.select_set(True)
  bpy.ops.object.delete()

  b_primitives["thin_spring"] = obj
  set_active_object_shade_smooth(obj)
  change_object_collection(obj, bpy.data.collections["Primitives"])
  
  return

def create_a_binding_site_primitives(b_params, args, b_primitives, materials):
  """Appends instructions to create an actin binding site primitive to b_script."""
  if args["render_mode"] == 0:
    # Indicates we should make biological binding sites (spheres).
    # Unavailable binding site.
    bpy.ops.mesh.primitive_ico_sphere_add(
      radius=b_params.radii["actin_bs_biological"],
      enter_editmode=False,
      location=(0,0,0),
      subdivisions=1
    )
    set_active_object_shade_smooth(bpy.context.object)
    b_primitives["unavail_bs"] = bpy.context.object
    b_primitives["unavail_bs"].name = "unavail_bs"
    b_primitives["unavail_bs"].data.materials.append(materials["bs_unavailable"])

    # Available binding site.
    bpy.ops.mesh.primitive_ico_sphere_add(
      radius=b_params.radii["actin_bs_biological"],
      enter_editmode=False,
      location=(0,0,0),
      subdivisions=1
    )
    set_active_object_shade_smooth(bpy.context.object)
    b_primitives["avail_bs"] = bpy.context.object
    b_primitives["avail_bs"].name = "avail_bs"
    b_primitives["avail_bs"].data.materials.append(materials["bs_available"])

  else:
    # Indicates we should make mechanical binding sites (cylinders).
    # Unavailable binding site.
    bpy.ops.mesh.primitive_cylinder_add(
      vertices=15,
      radius=b_params.radii["actin_bs_mechanical"], 
      depth=b_params.heights["actin_bs_mechanical"],
      enter_editmode=False, location=(0,0,0)
    )
    set_active_object_shade_smooth(bpy.context.object)
    b_primitives["unavail_bs"] = bpy.context.object
    b_primitives["unavail_bs"].name = "unavail_bs"
    b_primitives["unavail_bs"].data.materials.append(materials["bs_unavailable"])

    # Available binding site.
    bpy.ops.mesh.primitive_cylinder_add(
      vertices=15,
      radius=b_params.radii["actin_bs_mechanical"],
      depth=b_params.heights["actin_bs_mechanical"],
      enter_editmode=False,
      location=(0,0,0)
    )
    set_active_object_shade_smooth(bpy.context.object)
    b_primitives["avail_bs"] = bpy.context.object
    b_primitives["avail_bs"].name = "avail_bs"
    b_primitives["avail_bs"].data.materials.append(materials["bs_available"])

  change_object_collection(b_primitives["unavail_bs"], bpy.data.collections["Primitives"])
  change_object_collection(b_primitives["avail_bs"], bpy.data.collections["Primitives"])

  return

def create_titin_primitive(b_params, args, b_primitives, materials):
  """Appends instructions to b_script to create the titin primitives."""
  if args["render_mode"] == 0:
    # Indicates we want a biological representation of titin in the visualization.
    bpy.ops.mesh.primitive_cylinder_add(
      radius=b_params.radii["titin"],
      depth=1,
      enter_editmode=False,
      location=(0,0,0),
      rotation=(0,3.14159/2,0)
    )
    b_primitives["titin_cylinder"] = bpy.context.object
    b_primitives["titin_cylinder"].name = "titin_cylinder"
    b_primitives["titin_cylinder"].data.materials.append(materials["titin"])

  else:
    # Indicates we want a mechanical representation of titin in the visualization.
    # TODO: Convert this to a spring model instead of current titin filament cylinder.
    bpy.ops.mesh.primitive_cylinder_add(
      radius=b_params.radii["titin"],
      depth=1,
      enter_editmode=False,
      location=(0,0,0),
      rotation=(0,3.14159/2,0)
    )
    b_primitives["titin_cylinder"] = bpy.context.object
    b_primitives["titin_cylinder"].name = "titin_cylinder"
    b_primitives["titin_cylinder"].data.materials.append(materials["titin"])
  set_active_object_shade_smooth(b_primitives["titin_cylinder"])
  change_object_collection(b_primitives["titin_cylinder"], bpy.data.collections["Primitives"])

  return

def create_m_node_primitive(b_params, b_primitives, materials):
  """Creates the mechanical mesh primitive for a thick filament node."""
  bpy.ops.mesh.primitive_cylinder_add(
    radius=b_params.radii["thick"],
    depth=b_params.lengths["m_node"],
    enter_editmode=False,
    location=(0,0,0),
    vertices=64,
    rotation=(0,3.14159/2,0)
  )
  set_active_object_shade_smooth(bpy.context.object)
  b_primitives["m_node"] = bpy.context.object
  b_primitives["m_node"].name = "m_node"
  b_primitives["m_node"].data.materials.append(materials["thick"])
  change_object_collection(b_primitives["m_node"], bpy.data.collections["Primitives"])

  return

def create_a_node_primitive(b_params, b_primitives, materials):
  """Create the mechanical mesh primitive for a thin filament node."""
  bpy.ops.mesh.primitive_cylinder_add(
    radius=b_params.radii["thin"],
    depth=b_params.lengths["a_node"],
    vertices=64,
    enter_editmode=False,
    location=(0,0,0),
    rotation=(0,3.14159/2,0)
  )
  set_active_object_shade_smooth(bpy.context.object)
  b_primitives["a_node"] = bpy.context.object
  b_primitives["a_node"].name = "a_node"
  b_primitives["a_node"].data.materials.append(materials["thin"])
  change_object_collection(b_primitives["a_node"], bpy.data.collections["Primitives"])

  return

def create_z_line(half_sarcomere, b_params, object_dict, materials):
  """Creates a z-line."""
  # Create a box to act as a primitive z-line. This starts at the origin of the filament.
  z_coords = ([thick_obj.m_z for thick_obj in half_sarcomere.thick]
              + [thin_obj.a_z for thin_obj in half_sarcomere.thin])
  y_coords = ([thick_obj.m_y for thick_obj in half_sarcomere.thick]
              + [thin_obj.a_y for thin_obj in half_sarcomere.thin])
  max_y = max(y_coords)
  min_y = min(y_coords)
  max_z = max(z_coords)
  min_z = min(z_coords)
  
  box_center_y = (max_y + min_y) / 2.
  box_scale_y = abs(max_y) + abs(min_y) + 2 * b_params.boxes["padding"]
  box_center_z = (max_z + min_z) / 2.
  box_scale_z = abs(max_z) + abs(min_z) + 2 * b_params.boxes["padding"]
  
  bpy.ops.mesh.primitive_cube_add(size=1,
    enter_editmode=False, 
    location=(-0.5 * b_params.boxes["thickness"], box_center_y, box_center_z)
  )
  object_dict["z_line"] = bpy.context.object
  object_dict["z_line"].scale[0] = b_params.boxes["thickness"]
  object_dict["z_line"].scale[1] = box_scale_y
  object_dict["z_line"].scale[2] = box_scale_z
  object_dict["z_line"].data.materials.append(materials["z_line"])
  object_dict["z_line"].name = "z_line"

  return

def create_m_line(half_sarcomere, b_params, object_dict, materials):
  """Creates an m-line.
  
  Note: This is essentially just creating another z-line box like in the function create_z_line()
  but as my understanding of FiberSim matures and simulations grow in scale, encompassing more 
  than one myofibril, I will have to make sure this creates an m line for just a single myofiber
  at a time.
  """
  # Create a box to act as a primitive m-line.
  z_coords = ([thick_obj.m_z for thick_obj in half_sarcomere.thick]
              + [thin_obj.a_z for thin_obj in half_sarcomere.thin])
  y_coords = ([thick_obj.m_y for thick_obj in half_sarcomere.thick]
              + [thin_obj.a_y for thin_obj in half_sarcomere.thin])
  max_y = max(y_coords)
  min_y = min(y_coords)
  max_z = max(z_coords)
  min_z = min(z_coords)
  
  box_center_y = (max_y + min_y) / 2.
  box_scale_y = abs(max_y) + abs(min_y) + 2 * b_params.boxes["padding"]
  box_center_z = (max_z + min_z) / 2.
  box_scale_z = abs(max_z) + abs(min_z) + 2 * b_params.boxes["padding"]
  box_center_x = half_sarcomere.hs_length + 0.5 * b_params.boxes["thickness"]
  
  bpy.ops.mesh.primitive_cube_add(
    size=1,
    enter_editmode=False,
    location=(box_center_x, box_center_y, box_center_z)
  )
  object_dict["m_line"] = bpy.context.object
  object_dict["m_line"].scale[0] = b_params.boxes["thickness"]
  object_dict["m_line"].scale[1] = box_scale_y
  object_dict["m_line"].scale[2] = box_scale_z
  object_dict["m_line"].data.materials.append(materials["m_line"])
  object_dict["m_line"].name = "m_line"

  return 

def create_m_line_m_node_spring_primitive(b_params, b_primitives, materials):
  """Creates the spring that connects the M line with the first myosin node."""
  spring_core_gap = 0.6 * b_params.radii["thick"]
  offset = (-(0.5 * spring_core_gap + b_params.radii["thick_spring"]),0,0)
  spring_height = 5.5 * b_params.radii["thick_spring"]

  spring_data = bpy.data.meshes.new("m_line_m_node_spring_primitive")
  line_data = bpy.data.meshes.new("line")

  verts,edges = give_circle_vertices_and_edges(
    radius=b_params.radii["thick_spring"], 
    num_vertices=b_params.num_spring_vertices, 
    offset=offset
  )

  # height vector.
  line_verts = [(0,0,0.5), (0,0,-0.5)]
  line_edge = [(0,1)]

  spring_data.from_pydata(verts, edges, [])
  spring_data.update()
  line_data.from_pydata(line_verts, line_edge, [])
  line_data.update()

  # Create the object
  obj = bpy.data.objects.new("m_line_m_node_spring_primitive", spring_data)
  line_obj = bpy.data.objects.new("line", line_data)
  bpy.data.collections[0].objects.link(obj)
  bpy.data.collections[0].objects.link(line_obj)

  # Add the screw modifier
  obj.select_set(True)

  obj.modifiers.new("Spring",'SCREW')
  obj.modifiers["Spring"].object = bpy.data.objects["line"]
  obj.modifiers["Spring"].screw_offset = spring_height
  obj.modifiers["Spring"].iterations = b_params.lengths["lambda"] / spring_height
  obj.select_set(True)
  bpy.context.view_layer.objects.active = obj
  bpy.ops.object.modifier_apply(modifier="Spring")
  obj.data.materials.append(materials["thick_spring"])
  obj.rotation_euler[1] = -3.14159 / 2
  obj.scale[0] = 1.5
  obj.scale[1] = 1.5

  bpy.ops.object.select_all(action="DESELECT")
  line_obj.select_set(True)
  bpy.ops.object.delete()

  b_primitives["m_line_m_node_spring"] = obj
  set_active_object_shade_smooth(obj)
  change_object_collection(obj, bpy.data.collections["Primitives"])
  b_params.lengths["original_m_line_m_node_spring"] = obj.dimensions[2] - b_params.radii["thick_spring"]
  return

###################################################################################################
### Functions for Materials
###################################################################################################

def create_materials(args):
  """Appends instructions to create the materials for the model to b_script."""
  materials = {}

  # Form materials for coloring objects.
  if args["render_mode"] == 0:
    # This is the biological materials. These are taken from the POV-Ray repo Ken shared with me.
    materials["thick"] = bpy.data.materials.new(name="thick_material")
    materials["thick"].diffuse_color = (0.0026586, 0, 0.318539, 1)  # Deep blue

    materials["thick_spring"] = bpy.data.materials.new(name="thick_spring_material")
    materials["thick_spring"].diffuse_color = (0.0720458, 0.140784, 0.318539, 1)

    materials["thin"] = bpy.data.materials.new(name="thin_material")
    materials["thin"].diffuse_color = (0.310813, 0.135074, 0.176204, 1)  # Light purple (thistle)

    materials["thin_spring"] = bpy.data.materials.new(name="thin_spring_material")
    materials["thin_spring"].diffuse_color = (0.510566, 0.269025, 0.799338, 1) # Chalky purple.

    materials["z_line"] = bpy.data.materials.new(name="z_line_material")
    materials["z_line"].diffuse_color = (0.102474, 0.0216006, 0, 1)  # Brown

    materials["titin"] = bpy.data.materials.new(name="titin_material")
    materials["titin"].diffuse_color = (0.990667, 1, 0.00251775, 1)  # Yellow

    materials["tropomyosin"] = bpy.data.materials.new(name="tropomyosin_material")
    materials["tropomyosin"].diffuse_color = (0.827558, 0.00127092, 0.512953, 1)  # Purple

    materials["m_line"] = bpy.data.materials.new(name="m_line_material")
    materials["m_line"].diffuse_color = (0.147314, 0.147314, 0.147314, 1)  # Gray

    materials["bs_available"] = bpy.data.materials.new(name="bs_available_material")
    materials["bs_available"].diffuse_color = (0.0194498, 0.100481, 0.0159066, 1)  # Forest Green

    materials["bs_unavailable"] = bpy.data.materials.new(name="bs_unavailable_material")
    materials["bs_unavailable"].diffuse_color = (0.259799, 0.00109354, 0.00270012, 1)  # Red

    materials["background_plane"] = bpy.data.materials.new(name="background_plane_material")
    materials["background_plane"].diffuse_color = (1, 1, 1, 1) # Pure white.

    # Cross Bridge binding state materials
    materials["cb_state_0"] = materials["thick"]
    materials["cb_state_1"] = bpy.data.materials.new(name="binding_state_1_material")
    materials["cb_state_1"].diffuse_color = (0.0179366, 0.23372, 0.295693, 1)
    materials["cb_state_2"] = bpy.data.materials.new(name="binding_state_2_material")
    materials["cb_state_2"].diffuse_color = (0.854981, 0.233305, 0.00951545, 1)
    materials["cb_state_3"] = bpy.data.materials.new(name="binding_state_3_material")
    materials["cb_state_3"].diffuse_color = (0.854981, 0.00793651, 0.00508919, 1)
    materials["cb_state_4"] = bpy.data.materials.new(name="binding_state_4_material")
    #materials["cb_state_4"].diffuse_color = (0.854981, 0.00793651, 0.00508919, 1)
    materials["cb_state_5"] = bpy.data.materials.new(name="binding_state_5_material")
    #materials["cb_state_5"].diffuse_color = (0.854981, 0.00793651, 0.00508919, 1)
  else:
    # This is the mechanical materials. These are taken from the matlab files from within the FiberSim repo.
    materials["thick"] = bpy.data.materials.new(name="thick_material")
    materials["thick"].diffuse_color = (0.154682, 0.151359, 0.661864, 1)  # "Purpley blue"

    materials["thick_spring"] = bpy.data.materials.new(name="thick_spring_material")
    materials["thick_spring"].diffuse_color = (0, 0.925826, 0.901533, 1) # I'd call it a teal?

    materials["thin"] = bpy.data.materials.new(name="thin_material")
    materials["thin"].diffuse_color = (0.137584, 0.137584, 0.137584, 1)  # Gray.

    materials["thin_spring"] = bpy.data.materials.new(name="thin_spring_material")
    materials["thin_spring"].diffuse_color = (0.447979, 0.447979, 0.447979, 1) # Very light gray.

    materials["z_line"] = bpy.data.materials.new(name="z_line_material")
    materials["z_line"].diffuse_color = (0.102474, 0.0216006, 0, 1)  # Brown

    materials["titin"] = bpy.data.materials.new(name="titin_material")
    materials["titin"].diffuse_color = (0.990667, 1, 0.00251775, 1)  # Yellow

    materials["tropomyosin"] = bpy.data.materials.new(name="tropomyosin_material")
    materials["tropomyosin"].diffuse_color = (0.827558, 0.00127092, 0.512953, 1)  # Purple

    materials["m_line"] = bpy.data.materials.new(name="m_line_material")
    materials["m_line"].diffuse_color = (0.147314, 0.147314, 0.147314, 1)  # Gray

    materials["bs_available"] = bpy.data.materials.new(name="bs_available_material")
    materials["bs_available"].diffuse_color = (0.00323825, 0.496923, 0.00741916, 1)  # Lime (?) green.

    materials["bs_unavailable"] = bpy.data.materials.new(name="bs_unavailable_material")
    materials["bs_unavailable"].diffuse_color = (0.47699, 0.00462422, 0, 1)  # Red

    materials["background_plane"] = bpy.data.materials.new(name="background_plane_material")
    materials["background_plane"].diffuse_color = (1, 1, 1, 1) # Pure white.

    # Cross Bridge binding state materials
    materials["cb_state_0"] = materials["thick"]
    materials["cb_state_1"] = bpy.data.materials.new(name="binding_state_1_material")
    materials["cb_state_1"].diffuse_color = (0.00153204, 0.280525, 0.855469, 1)
    materials["cb_state_2"] = bpy.data.materials.new(name="binding_state_2_material")
    materials["cb_state_2"].diffuse_color = (0.37575, 0.00625245, 0.1439, 1)
    materials["cb_state_3"] = bpy.data.materials.new(name="binding_state_3_material")
    materials["cb_state_3"].diffuse_color = (0.854981, 0.241046, 0.318406, 1)
    materials["cb_state_4"] = bpy.data.materials.new(name="binding_state_4_material")
    materials["cb_state_4"].diffuse_color = (0.854981, 0.00793651, 0.00508919, 1)
    materials["cb_state_5"] = bpy.data.materials.new(name="binding_state_5_material")
    materials["cb_state_5"].diffuse_color = (0.20607, 0.221692, 0.384375, 1)

  # if args.render_mode == 1:
  if False:
    # Make all of the materials metallic for fun since it's "mechanical".
    materials["thick"].metallic = 0.55
    materials["thick"].roughness = 1
    materials["thin"].metallic = 0.55
    materials["thin"].roughness = 1
    materials["z_line"].metallic = 0.55
    materials["z_line"].roughness = 1
    materials["m_line"].metallic = 0.55
    materials["m_line"].roughness = 1
    materials["titin"].metallic = 0.55
    materials["titin"].roughness = 1
    materials["tropomyosin"].metallic = 0.55
    materials["tropomyosin"].roughness = 1
    materials["bs_available"].metallic = 0.55
    materials["bs_available"].roughness = 1
    materials["bs_unavailable"].metallic = 0.55
    materials["bs_unavailable"].roughness = 1

  # Make every material the same roughness and specularity.
  for key, mat in materials.items():
    mat.specular_intensity = 0.5
    mat.roughness = 0.5

  return materials

def update_object_material(obj, new_material):
  """Checks if the object needs to update its material and updates it if so."""
  if obj.data.materials[0] != new_material:
    obj.data.materials.pop()
    obj.data.materials.append(new_material)

def set_shading_smooth():
  """Set the shading to smooth for all objects."""
  bpy.ops.object.select_all(action='SELECT')
  # Set the shading as smooth
  bpy.ops.object.shade_smooth()
  bpy.ops.object.select_all(action='DESELECT')

  return

###################################################################################################
### Functions for Creating Myofilament Structures from Primitives
###################################################################################################

def create_thick_filament(half_sarcomere, thick_index, b_params, args, materials, obj_lists, 
                          b_primitives):
  """Creates a thick filament.
  
  This is called for each thick filament to be created, creating the filament by sequentially
  copying mesh primitives created in this script.
  """
  thick_obj = half_sarcomere.thick[thick_index]
  thick_key = give_thick_key(thick_obj)

  if thick_obj.id % 5 == 0:
    print ("Drawing thick filament #{}".format(thick_obj.id))
  
  # Get list of adjacent filament keys.
  adjacent_thin_fil_indexes = find_adjacent_thin_fils_to_thick_fils(half_sarcomere, thick_index)

  # Get x of the last cross bridge.
  m_fin_x = thick_obj.cb_x[-1]
  m_first_x = thick_obj.cb_x[0]
  m_fil_length = m_first_x - m_fin_x
  
  # Draw titin.
  create_titin(half_sarcomere, thick_index, obj_lists, b_primitives, args, b_params, 
    adjacent_thin_fil_indexes)

  # Draw a spring from the M line to the first myosin node.
  create_m_line_m_node_spring(half_sarcomere, thick_index, obj_lists, b_params, b_primitives)

  cbs_per_node = 6

  if args["render_mode"] == 0:
    # Create the thick filament as a biological representation.
    raise RuntimeError("Biochemical rendering not supported currently.")
    # bpy.ops.mesh.primitive_cylinder_add(
    #   radius=b_params.radii["thick"], 
    #   depth=m_fil_length, 
    #   enter_editmode=False,
    #   location=(m_fil_length / 2 + m_fin_x, dump_dict[key]['y'], dump_dict[key]['z']), 
    #   rotation=(0, 3.14159/2, 0)
    # )
    # obj_lists[key] = bpy.context.object
    # obj_lists[key].data.materials.append(materials["thick"])
    # obj_lists[key].name = key

  elif args["render_mode"] == 1:
    # Create each node in the filament.
    cb_indexes = np.arange(0, thick_obj.m_no_of_cbs, cbs_per_node)
    name_template = thick_key + "_{}_{}"
    for i,cb_index in enumerate(cb_indexes):
      this_node = b_primitives["m_node"].copy()
      this_node.name = name_template.format("node", i+1)
      obj_lists[thick_key].append(this_node)

      this_node.location.x = thick_obj.cb_x[cb_index]
      this_node.location.y = thick_obj.m_y
      this_node.location.z = thick_obj.m_z

      # Draw springs from node-to-node.
      if cb_index != cb_indexes[-1]:
        this_spring = b_primitives["thick_spring"].copy()
        this_spring.name = name_template.format("thick_spring", i+1)
        obj_lists[thick_key].append(this_spring)

        cb_n_0 = cb_index
        cb_n_1 = cb_indexes[i+1]
        # The spring's origin is actually at the bottom of the spring instead of in the middle
        #   so we just start at the first node.
        spring_x = thick_obj.cb_x[cb_n_1]
        spring_x_length = (thick_obj.cb_x[cb_n_0] - thick_obj.cb_x[cb_n_1] 
          + 2 * b_params.radii["thick_spring"])

        this_spring.dimensions.z = spring_x_length
        this_spring.location.x = spring_x
        this_spring.location.y = thick_obj.m_y
        this_spring.location.z = thick_obj.m_z

  # Draw each myosin head.
  for cb_index in range(len(thick_obj.cb_x)):
    create_m_head(half_sarcomere, thick_index, cb_index, b_params, b_primitives, obj_lists, args,
      adjacent_thin_fil_indexes)

  return

def create_m_head(half_sarcomere, thick_index, cb_index, b_params, b_primitives, obj_lists, args,
    adjacent_thin_fil_indexes):
  """Creates a myosin head."""
  # Get a copy of all of the object primitives for this myosin head.
  thick_obj = half_sarcomere.thick[thick_index]
  thick_key = give_thick_key(thick_obj)
  name_template = thick_key + "_{}_" + str(cb_index)
  cb_dict = b_primitives["cb_state_"+str(thick_obj.cb_state[cb_index])]

  this_stalk = cb_dict["stalk"].copy()
  this_stalk.name = name_template.format("stalk")
  obj_lists[thick_key].append(this_stalk)
  
  this_cat = cb_dict["catalytic_domain"].copy()
  this_cat.name = name_template.format("cat")
  obj_lists[thick_key].append(this_cat)
  
  this_arm = cb_dict["lever_arm"].copy()
  this_arm.name = name_template.format("arm")
  obj_lists[thick_key].append(this_arm)
  
  this_spring = cb_dict["spring"].copy()
  this_spring.name = name_template.format("spring")
  obj_lists[thick_key].append(this_spring)

  update_m_head_location(half_sarcomere, thick_index, cb_index, b_params, args, this_stalk, 
                         this_arm, this_cat, this_spring, adjacent_thin_fil_indexes)

def create_m_line_m_node_spring(half_sarcomere, thick_index, obj_lists, b_params, b_primitives):
  """Creates the spring that goes from the M line to the first myosin node."""
  # Copy the m_line_m_node primitive.
  this_spring = b_primitives["m_line_m_node_spring"].copy()
  thick_key = give_thick_key(half_sarcomere.thick[thick_index])
  this_spring.name = thick_key + "_m_line_m_node_spring"
  obj_lists[thick_key].append(this_spring)

  # Update the location and stretch of the spring.
  update_m_line_m_node_spring(half_sarcomere, thick_index, b_params)
  return

def create_thin_filament(half_sarcomere, thin_index, b_params, args, materials, obj_lists, 
                         b_primitives):
  """Creates a thin filament.
  
  This is called for each thin filament to be created, creating the filament by sequentially
  copying mesh primitives created in the blender_primitives.py script.
  """
  # Print progress.
  thin_obj = half_sarcomere.thin[thin_index]
  thin_key = "a_" + str(thin_obj.id)

  if thin_obj.id % 5 == 0:
    print ("Drawing thin filament #{}".format(thin_obj.id))
  
  # Create each node in the filament.
  no_bs_per_node = 2
  bs_indexes = np.arange(0, thin_obj.a_no_of_bs, no_bs_per_node)
  name_template = thin_key + "_{}_{}"
  for i, bs_index in enumerate(bs_indexes):
    this_node = b_primitives["a_node"].copy()
    this_node.location.x = thin_obj.bs_x[bs_index]
    this_node.location.y = thin_obj.a_y
    this_node.location.z = thin_obj.a_z
    this_node.name = name_template.format("node", i)
    obj_lists[thin_key].append(this_node)

    # Draw springs from node-to-node.
    if bs_index != bs_indexes[-1]:
      bs_n_0 = bs_index
      bs_n_1 = bs_indexes[i + 1]
      # The spring's origin is actually at the bottom of the spring instead of in the middle
      #   so we just start at the first node.
      spring_x = thin_obj.bs_x[bs_n_0]
      spring_x_length = (thin_obj.bs_x[bs_n_1] - thin_obj.bs_x[bs_n_0]
                        + 3 * b_params.radii["thin_spring"])
      this_spring = b_primitives["thin_spring"].copy()
      this_spring.dimensions.z = spring_x_length
      this_spring.location.x = spring_x
      this_spring.location.y = thin_obj.a_y
      this_spring.location.z = thin_obj.a_z
      this_spring.name = name_template.format("spring", i)
      obj_lists[thin_key].append(this_spring)

  # Create each binding site
  for bs_index in range(len(thin_obj.bs_x)):
    if args["render_mode"] == 0:
      # Get the y and z projections since the binding sites jut out at each node.
      bs_proj_y = (b_params.radii["actin_bs_biological"] 
                   * np.sin(thin_obj.bs_angle[bs_index] * np.pi / 180))
      bs_proj_z = (b_params.radii["actin_bs_biological"] 
                   * np.cos(thin_obj.bs_angle[bs_index] * np.pi / 180))

      # Get the material of the binding site based on availability
      if thin_obj.bs_state[bs_index] == 1:
        # Indicates the site is unavailable for binding.
        bs_mesh = b_primitives["unavail_bs"]
      elif thin_obj.bs_state[bs_index] == 2:
        # Indicates the site is available for binding.
        bs_mesh = b_primitives["avail_bs"]

      # Make a copy of the original binding site and append to the script.
      this_bs = bs_mesh.copy()
      this_bs.location.x = thin_obj.bs_x[bs_index]
      this_bs.location.y = thin_obj.a_y + bs_proj_y
      this_bs.location.z = thin_obj.a_z + bs_proj_z
      this_bs.name = thin_key + "_bs_" + str(bs_index)
      obj_lists[thin_key].append(this_bs)

    else:
      # Get the y and z projections since the binding sites jut out at each node.
      bs_proj_y = ((b_params.radii["thin"] + 0.5 * b_params.heights["actin_bs_mechanical"]) 
                    * np.sin(thin_obj.bs_angle[bs_index] * np.pi / 180))
      bs_proj_z = ((b_params.radii["thin"] + 0.5 * b_params.heights["actin_bs_mechanical"]) 
                    * np.cos(thin_obj.bs_angle[bs_index] * np.pi / 180))
      
      # Get the material of the binding site based on availability
      if thin_obj.bs_state[bs_index] == 0:
        # Indicates the site is unavailable for binding.
        bs_mesh = b_primitives["unavail_bs"]
      elif thin_obj.bs_state[bs_index] == 1:
        # Indicates the site is available for binding.
        bs_mesh = b_primitives["avail_bs"]
      this_bs = bs_mesh.copy()
      this_bs.location.x = thin_obj.bs_x[bs_index]
      this_bs.location.y = thin_obj.a_y + bs_proj_y
      this_bs.location.z = thin_obj.a_z + bs_proj_z
      this_bs.rotation_euler.x = -(thin_obj.bs_angle[bs_index] * np.pi / 180.)
      this_bs.name = thin_key + "_bs_" + str(bs_index)
      obj_lists[thin_key].append(this_bs)

  
  # Now we create the tropomyosin
  # NOTE: Leaving this for later when it's explicitly incorporated into the model.
  #create_tropomyosin(dump_dict, b_script, b_params)

  return

def create_titin(half_sarcomere, thick_index, 
  # key, dump_dict, 
  obj_lists, b_primitives, args, b_params, adjacent_thin_fil_indexes
  ):
  """Draws titin connected from the Z-line + offset to all 8-neighborhood thick filaments."""
  print ("ASSUMING TITIN OFFSET OF ZERO HERE!")
  t_offset = 0
  # Get the end coordinates for this thick filament.
  thick_obj = half_sarcomere.thick[thick_index]
  end_vec = np.asarray([thick_obj.cb_x[-1], thick_obj.m_y, thick_obj.m_z])

  thick_key = give_thick_key(thick_obj)
  t_name = "titin_"+thick_key+"_a_{}"

  # Draw titin for each adjacent thick filament.
  for thin_fil_index in adjacent_thin_fil_indexes:
    this_t_name = t_name.format(thin_fil_index)
    this_titin = b_primitives["titin_cylinder"].copy()
    this_titin.name = this_t_name
    obj_lists[thick_key].append(this_titin)
    this_thin_obj = half_sarcomere.thin[thin_fil_index]
    base_vec = np.asarray([-0.5 * b_params.boxes["thickness"] + t_offset, this_thin_obj.a_y, 
      this_thin_obj.a_z])
    position_scale_orient_cylinder(
      cyl_obj=this_titin,
      base_vec=base_vec,
      end_vec=end_vec
    )




  # Get the end coordinates for this thick filament.
  # end_vec = np.asarray([
  #   dump_dict[key]["cb_x"][-1],
  #   dump_dict[key]['y'],
  #   dump_dict[key]['z']
  # ])

  # t_name = "titin_"+key+"_a_{}"

  # # Draw titin for each adjacent thick filament.
  # for thin_fil_ID in adjacent_thin_filament_IDs:
  #   this_t_name = t_name.format(thin_fil_ID)
  #   this_titin = b_primitives["titin_cylinder"].copy()
  #   this_titin.name = this_t_name
  #   obj_lists[key].append(this_titin)
  #   this_thin_key = "a_{}".format(thin_fil_ID)
  #   base_vec = np.asarray([
  #     -0.5 * b_params.boxes["thickness"],
  #     dump_dict[this_thin_key]['y'],
  #     dump_dict[this_thin_key]['z']
  #   ])
  #   position_scale_orient_cylinder(
  #     cyl_obj=this_titin,
  #     base_vec=base_vec,
  #     end_vec=end_vec
  #   )
  return

def create_tropomyosin(dump_dict, b_script, b_params):
  """Appends the instructions to create a tropomyosin filament to the b_script file."""
  
  return b_script

###################################################################################################
### Functions for Primitive Control
###################################################################################################

def create_blender_primitives(b_params, args, materials):
  """Creates the mesh primitives for the model."""
  b_primitives = dict()
  create_a_binding_site_primitives(b_params, args, b_primitives, materials)
  create_titin_primitive(b_params, args, b_primitives, materials)
  create_basic_myosin_head_primitives(b_params, args, b_primitives, materials)

  if args["render_mode"] == 1:
    # Indicates this is a mechanical rendering and we need the nodes and springs.
    create_m_node_primitive(b_params, b_primitives, materials)
    create_a_node_primitive(b_params, b_primitives, materials)
    create_thick_spring_primitive(b_params, b_primitives, materials)
    create_thin_spring_primitive(b_params, b_primitives, materials)
    create_m_line_m_node_spring_primitive(b_params, b_primitives, materials)

  # Shade everything as smooth if we're using the biological model to make things pretty.
  # if args["render_mode"] == 0:
  if True:
    set_shading_smooth()

  return b_primitives

def change_object_collection(obj, new_collection):
  """Links obj to new_collection and unlinks obj from current collection.
  
  Note: It's important to do this after you've finished modifying the primitive how you want since
  this invalidates any bpy.opys.object calls on the object (I think).
  """
  new_collection.objects.link(obj)
  obj.users_collection[0].objects.unlink(obj)

  # I have to deselect the object for some reason?
  obj.select_set(False)
  return

###################################################################################################
### Functions for Moving Filament Primitives
###################################################################################################

def update_thick_filament(half_sarcomere, thick_index, b_params, args):
  """Updates the location, orientation, and states of the blender primitives for thick filament."""
  # Find the keys of the adjacent filaments.
  # adjacent_thin_filament_IDs = find_adjacent_thin_fils_to_thick_fils(key, dump_dict, 'a')
  adjacent_thin_fil_indexes = find_adjacent_thin_fils_to_thick_fils(half_sarcomere, thick_index)

  thick_obj = half_sarcomere.thick[thick_index]
  thick_key = give_thick_key(thick_obj)

  # Update the spring between the M line and the first myosin node.
  update_m_line_m_node_spring(half_sarcomere, thick_index, b_params)

  # Update the titins connected to this myosin filament.
  update_titin_location(half_sarcomere, thick_index, b_params, adjacent_thin_fil_indexes)

  # Update the node locations and the spring stretch between the nodes.
  cbs_per_node = 6
  cb_indexes = np.arange(0, half_sarcomere.thick[thick_index].m_no_of_cbs, cbs_per_node)

  for i,cb_index in enumerate(cb_indexes):
    name_template = thick_key + "_{}_" + str(i+1)
    # Update the node.
    this_node = bpy.data.objects[name_template.format("node")]
    this_node.location.x = thick_obj.cb_x[cb_index]
    this_node.location.y = thick_obj.m_y
    this_node.location.z = thick_obj.m_z

    # Update the spring.
    if cb_index != cb_indexes[-1]:
      this_spring = bpy.data.objects[name_template.format("spring")]
      cb_n_0 = cb_index
      cb_n_1 = cb_indexes[i + 1]
      # The spring's origin is actually at the bottom of the spring instead of in the middle
      #   so we just start at the first node.
      spring_x = thick_obj.cb_x[cb_n_1]
      spring_x_length = (thick_obj.cb_x[cb_n_0] - thick_obj.cb_x[cb_n_1]
                         + 2 * b_params.radii["thick_spring"])
      this_spring.dimensions.z = spring_x_length
      this_spring.location.x = spring_x
      this_spring.location.y = thick_obj.m_y
      this_spring.location.z = thick_obj.m_z

  # Update the cross bridge and node/spring locations and scales.
  object_template = thick_key + "_{}_{}"
  for cb_index in range(len(thick_obj.cb_x)):
    # Get the objects to update.
    ID = str(cb_index)
    this_stalk = bpy.data.objects[object_template.format("stalk", ID)]
    this_arm = bpy.data.objects[object_template.format("arm", ID)]
    this_cat = bpy.data.objects[object_template.format("cat", ID)]
    this_spring = bpy.data.objects[object_template.format("spring", ID)]

    # Update the materials if necessary.
    this_binding_state = str(thick_obj.cb_state[cb_index])
    update_m_head_mesh(this_stalk, "stalk", this_binding_state)
    update_m_head_mesh(this_arm, "arm", this_binding_state)
    update_m_head_mesh(this_cat, "cat", this_binding_state)
    update_m_head_mesh(this_spring, "spring", this_binding_state)

    update_m_head_location(half_sarcomere, thick_index, cb_index, b_params, args, this_stalk, 
                           this_arm, this_cat, this_spring, adjacent_thin_fil_indexes)

  return

def update_titin_location(half_sarcomere, thick_index, b_params, adjacent_thin_fil_indexes):
  """Updates the location/stretch of titin at each time point."""
  print ("ASSUMING TITIN OFFSET OF ZERO HERE!")
  t_offset = 0
  # Get the end coordinates for this thick filament.
  thick_obj = half_sarcomere.thick[thick_index]
  end_vec = np.asarray([thick_obj.cb_x[-1], thick_obj.m_y, thick_obj.m_z])

  thick_key = give_thick_key(thick_obj)
  t_name = "titin_"+thick_key+"_a_{}"

  # Draw titin for each adjacent thick filament.
  for thin_fil_index in adjacent_thin_fil_indexes:
    this_titin = bpy.data.objects[t_name.format(thin_fil_index)]
    this_thin_obj = half_sarcomere.thin[thin_fil_index]
    base_vec = np.asarray([-0.5 * b_params.boxes["thickness"] + t_offset, this_thin_obj.a_y, 
      this_thin_obj.a_z])
    position_scale_orient_cylinder(
      cyl_obj=this_titin,
      base_vec=base_vec,
      end_vec=end_vec
    )
  return

def update_m_head_location(half_sarcomere, thick_index, cb_index, b_params, args, this_stalk, 
                           this_arm, this_cat, this_spring, adjacent_thin_fil_indexes):
  """Updates the location, scale, and state of the cross bridge in the model."""
  thick_obj = half_sarcomere.thick[thick_index]
  m_y = thick_obj.m_y
  m_z = thick_obj.m_z
  adjacent_thin_fil_IDs = [thin_obj.id for i, thin_obj in enumerate(half_sarcomere.thin) if i in 
    adjacent_thin_fil_indexes]

  # Always draw a stalk
  
  cb_y = ((0.5 * b_params.heights["cb_stalk"] + b_params.radii["thick"]) 
              * np.sin(thick_obj.cb_angle[cb_index] * np.pi / 180) + m_y)
  cb_z = ((0.5 * b_params.heights["cb_stalk"] + b_params.radii["thick"])
              * np.cos(thick_obj.cb_angle[cb_index] * np.pi / 180) + m_z)

  rot_angle = -thick_obj.cb_angle[cb_index] * np.pi / 180

  this_stalk.location.x = thick_obj.cb_x[cb_index]
  this_stalk.location.y = cb_y
  this_stalk.location.z = cb_z
  this_stalk.rotation_euler.x = rot_angle

  # Figure out where the catalytic domain, arm, and spring are.
  # No matter what, the arm Y and Z coords are the same (just varies in X).
  arm_base_y = cb_y
  arm_base_z = cb_z

  is_mirrored_cb = thick_obj.cb_bound_to_a_f[cb_index] not in adjacent_thin_fil_IDs
  draw_mirrored_cb = (is_mirrored_cb and args["draw_mirrored_cb_connections"]) or not is_mirrored_cb
  if (thick_obj.cb_bound_to_a_n[cb_index] >= 0 and draw_mirrored_cb):
    # This cb is bound.
    actin_filament_key = "a_" + str(thick_obj.cb_bound_to_a_f[cb_index])
    actin_filament_index = thick_obj.cb_bound_to_a_f[cb_index]
    binding_site_index = thick_obj.cb_bound_to_a_n[cb_index]

    # Since the cb is bound, it will always end coplanar (same YZ plane) with the actin binding
    # site.
    thin_obj = half_sarcomere.thin[actin_filament_index]
    arm_fin_x = thin_obj.bs_x[binding_site_index]
    cat_domain_x = arm_fin_x

    # Figure out bio/mechanical rendering coordinates
    if args["render_mode"] == 0:
      # Skip the biological rendering for now
      raise RuntimeError("Biochemical rendering not supported yet.")

    else:
      cat_domain_y = (thin_obj.a_y
                      + (b_params.heights["actin_bs_mechanical"] + b_params.radii["thin"]
                         + b_params.heights["cb_binding_site"])
                        * np.sin(thin_obj.bs_angle[binding_site_index] * np.pi / 180))
      cat_domain_z = (thin_obj.a_z
                      + (b_params.heights["actin_bs_mechanical"] + b_params.radii["thin"]
                         + b_params.heights["cb_binding_site"])
                        * np.cos(thin_obj.bs_angle[binding_site_index] * np.pi / 180))
      cat_rot_angle = (180 - thin_obj.bs_angle[binding_site_index]) * np.pi / 180

      # Find the polarity of the muscle so that the positive stretch due to cross bridge state is 
      # correct.
      if bpy.data.objects["z_line"].location.x < bpy.data.objects["m_line"].location.x:
        polarity = -1
      else:
        polarity = 1

      # Find the base X coordinate for the cb.
      arm_base_x = (thin_obj.bs_x[binding_site_index] 
                    + polarity * b_params.extensions[thick_obj.cb_state[cb_index]])

  else:
    # This cb is unbound.
    cat_domain_x = (thick_obj.cb_x[cb_index] 
                    + (np.random.rand() - 0.5) * b_params.lengths["cb_rand_x"])

    # Check if the myosin head is "super relaxed".
    # Temporarily turning off the "super relaxed" code since state 5 is not super relaxed.
    # if dump_dict[key][cb_key]["state"] != 5:
    if True:
      # Add in the random angle
      cat_domain_angle = (thick_obj.cb_angle[cb_index] 
                          + (np.random.rand() - 0.5) 
                            * b_params.angles["cb_separation"]) * np.pi / 180
      cat_domain_y = (thick_obj.m_y
                      + np.sin(cat_domain_angle)
                        * (b_params.heights["cb_stalk"]
                           + b_params.heights["cb_lever_arm"]
                           + 0.5 * b_params.heights["cb_binding_site"]
                           + b_params.radii["thick"]))
      cat_domain_z = (thick_obj.m_z 
                      + np.cos(cat_domain_angle)
                        * (b_params.heights["cb_stalk"]
                           + b_params.heights["cb_lever_arm"]
                           + 0.5 * b_params.heights["cb_binding_site"]
                           + b_params.radii["thick"]))
    else:
      cat_domain_angle = thick_obj.cb_angle[cb_index] * np.pi / 180
      cat_domain_y = (thick_obj.m_y
                      + np.sin(cat_domain_angle) 
                        * (b_params.heights["cb_stalk"]
                           + 0.5 * b_params.heights["cb_binding_site"]
                           + b_params.radii["thick"]))
      cat_domain_y = (thick_obj.m_z
                      + np.cos(cat_domain_angle) 
                        * (0.5 * b_params.heights["cb_stalk"]
                           + 0.5 * b_params.heights["cb_binding_site"]
                           + b_params.radii["thick"]))
      cat_domain_x = thick_obj.cb_x[cb_index] - 2.5 * b_params.heights["cb_lever_arm"]
    
    cat_rot_angle = rot_angle
    arm_base_x = thick_obj.cb_x[cb_index]

  arm_fin_x = cat_domain_x
  arm_fin_y = cat_domain_y
  arm_fin_z = cat_domain_z
  
  this_cat.location.x = cat_domain_x
  this_cat.location.y = cat_domain_y
  this_cat.location.z = cat_domain_z
  this_cat.rotation_euler.x = cat_rot_angle
  
  # Use convenient function to place the lever arm.
  position_scale_orient_cylinder(
    this_arm, 
    np.asarray([arm_base_x, arm_base_y, arm_base_z]),
    np.asarray([arm_fin_x, arm_fin_y, arm_fin_z])
  )
  
  # We now have to draw the spring.
  # The spring has to be placed differently since the origin is at the start of the spring
  spring_loc = [thick_obj.cb_x[cb_index], cb_y, cb_z]
  spring_extension = arm_base_x - spring_loc[0]
  spring_scale = spring_extension / b_params.lengths["original_cb_spring"]
  this_spring.location.x = spring_loc[0]
  this_spring.location.y = spring_loc[1]
  this_spring.location.z = spring_loc[2]
  this_spring.scale[2] = spring_scale

  return

def update_m_head_mesh(obj, obj_string, binding_state):
  """Checks and updates the cross bridge mesh in case it transitioned at this time step."""
  this_str = obj_string+"_cb_state_" + binding_state
  if obj.data != bpy.data.objects[this_str].data:
    obj.data = bpy.data.objects[this_str].data
  return

def update_m_line_location(dump_dict, b_params, args):
  """Updates the location of the M line based on the half sarcomere length."""
  print ("Need to make this update for all multiple half sarcomeres.")
  bpy.data.objects["m_line"].location[0] = dump_dict["header"]["hsl"]
  return

def update_m_line_m_node_spring(half_sarcomere, thick_index, b_params):
  """Updates spring location and stretch located between the M line and the first myosin node."""
  thick_key = give_thick_key(half_sarcomere.thick[thick_index])
  this_spring = bpy.data.objects[thick_key+"_m_line_m_node_spring"]
  spring_loc = [
    half_sarcomere.hs_length + 0.5 * b_params.boxes["thickness"],
    half_sarcomere.thick[thick_index].m_y,
    half_sarcomere.thick[thick_index].m_z
  ]
  spring_extension = ( spring_loc[0]
                     - half_sarcomere.thick[thick_index].cb_x[0]
                     + 1.25 * b_params.radii["thick_spring"])
  spring_scale = spring_extension / b_params.lengths["original_m_line_m_node_spring"]
  this_spring.location.x = spring_loc[0]
  this_spring.location.y = spring_loc[1]
  this_spring.location.z = spring_loc[2]
  this_spring.scale[2] = spring_scale
  return

def update_thin_filament(half_sarcomere, thin_index, b_params, args):
  """Updates the location, orientation, and states of the blender primitives for thin filament."""
  thin_obj = half_sarcomere.thin[thin_index]
  thin_key = give_thin_key(thin_obj)
  # Print progress.
  if thin_obj.id % 5 == 0:
    print ("Updating thin filament #{}".format(thin_obj.id))
  
  # Create each node in the filament.
  no_of_bs_per_node = 2
  bs_indexes = np.arange(0, thin_obj.a_no_of_bs, no_of_bs_per_node)
  name_template = thin_key + "_{}_{}"
  for i,bs_index in enumerate(bs_indexes):
    this_node = bpy.data.objects[name_template.format("node", i)]
    this_node.location.x = thin_obj.bs_x[bs_index]
    this_node.location.y = thin_obj.a_y
    this_node.location.z = thin_obj.a_z

    # Draw springs from node-to-node.
    if bs_index != bs_indexes[-1]:
      bs_n_0 = bs_index
      bs_n_1 = bs_indexes[i + 1]
      # The spring's origin is actually at the bottom of the spring instead of in the middle
      #   so we just start at the first node.
      spring_x = thin_obj.bs_x[bs_n_0]
      spring_x_length = (thin_obj.bs_x[bs_n_1] - thin_obj.bs_x[bs_n_0]
                         + 3 * b_params.radii["thin_spring"])
      this_spring = bpy.data.objects[name_template.format("spring", i)]
      this_spring.dimensions.z = spring_x_length
      this_spring.location.x = spring_x
      this_spring.location.y = thin_obj.a_y
      this_spring.location.z = thin_obj.a_z

  # Update each binding site.
  # NOTE: The nodes don't scale in the YZ plane, so we just need to move along the x axis.
  bs_name_template = thin_key + "_bs_{}"
  for bs_index in range(len(thin_obj.bs_x)):
    this_bs = bpy.data.objects[bs_name_template.format(bs_index)]
    this_bs.location.x = thin_obj.bs_x[bs_index]
    
    # Update the material.
    if thin_obj.bs_state[bs_index] == 0:
      material_string = "bs_unavailable_material"
    elif thin_obj.bs_state[bs_index] == 1:
      material_string = "bs_available_material"
    this_time_step_material = bpy.data.materials[material_string]
    update_object_material(this_bs, this_time_step_material)  
  
  # Now we update the tropomyosin
  # NOTE: Leaving this for later when it's explicitly incorporated into the model.
  #update_tropomyosin(dump_dict, b_params)
  
  return