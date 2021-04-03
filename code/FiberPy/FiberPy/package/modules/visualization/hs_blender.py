# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 22:22:51 2021

@author: ken
"""

# Blender objects
import bpy
import bmesh

import numpy as np

# Half-sarcomere
import half_sarcomere as hs

class hs_blender():
    """ Class for a half-sarcomere in blender """

    def __init__(self, hs, frame, template, blender):

        # Create local copies
        self.hs = hs
        self.frame = frame
        self.template = template
        self.blender = blender

        # Set up
        self.setup_blender_script()
        self.setup_render_engine()

        # Create obj lists
        self.obj_lists = dict()

        self.create_lists_and_collection()

        # Create a dictionary to hold blender primitives
        self.hs_b = dict()
        self.hs_b['primitives'] = dict()

        # Create the primitives
        self.create_primitive_m_crown()

        # Create thick filaments
        self.b_obj = dict()
        for i, thick_fil in enumerate(self.hs.thick_fil):
            self.b_obj[('m_%i' % i)] = []
            self.create_thick_filament(i)

        # Create a camera
        self.create_camera()

        # Create lights
        self.create_lights()

        # Link everything
        self.link_all_objects()

        # Render
        self.render_screenshot()

    def setup_blender_script(self):
        """Adds boilerplate blender scripting instructions."""
        # Set mode to object mode since we have to do most of the operations there.

        if bpy.context.active_object != None and \
                bpy.context.active_object.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')

        # Delete objects if they currently exist.
        if bpy.data.objects != []:
            bpy.ops.object.select_all(action='SELECT')
            bpy.ops.object.delete()
  
        # Set the viewport background color to white.
        bpy.context.preferences.themes['Default']. \
            view_3d.space.gradients.high_gradient.hsv = (0.0, 0.0, 1.0)
        bpy.context.scene.world.use_nodes = False
        bpy.context.scene.world.color = (1, 1, 1)

    def setup_render_engine(self):
        """Sets up the render engine that we want to use."""
        # Setting the render samples to 8 instead of 64 to cut render time.
        # Note: I don't think this affects much since the objects aren't very reflective.
        # bpy.context.scene.eevee.taa_render_samples = 8

        if self.blender["render_quality"] == "high":
            bpy.context.scene.render.engine = 'CYCLES'
        elif self.blender["render_quality"] == "medium":
            # Opting for Blender's workbench render engine because it does a good job and it's fast.
            bpy.context.scene.render.engine = 'BLENDER_WORKBENCH'
            bpy.context.scene.display.shading.light = 'STUDIO'
            bpy.context.scene.view_settings.view_transform = 'Filmic'
            bpy.context.scene.display.shading.show_object_outline = True
            bpy.context.scene.view_settings.exposure = 1.5
            bpy.context.scene.view_settings.gamma = 1.1
            bpy.context.scene.display.shading.show_specular_highlight = True
            bpy.context.scene.display.shading.shadow_intensity = 0.1

    def create_camera(self):
        """ Creates a camera """

        loc = (self.frame['camera']['location']['x'],
               self.frame['camera']['location']['y'],
               self.frame['camera']['location']['z'])

        rot = (self.frame['camera']['rotation']['x'] * np.pi/180,
               self.frame['camera']['rotation']['y'] * np.pi/180,
               self.frame['camera']['rotation']['z'] * np.pi/180)

        bpy.ops.object.camera_add(
            enter_editmode=False,
            align='VIEW',
            location=loc,
            rotation=rot)
        bpy.context.object.data.lens = 35

        bpy.context.scene.camera = bpy.context.object

    def create_lights(self):
        """ Create lights around the half-sarcomere """
        print("lights")

    def create_lists_and_collection(self):
        """ Lists everything together """

        thick = bpy.data.collections.new(name='Thick collection')
        bpy.context.scene.collection.children.link(thick)

        for i, thick_fil in enumerate(self.hs.thick_fil):
            thick_id = ('m_%i' % i)
            self.obj_lists[thick_id] = []
            col = bpy.data.collections.new(thick_id)
            thick.children.link(col)
    
    def link_all_objects(self):
        """Links all objects in `obj_lists` to the current collection
  
        Note: Collection is a Blender term and is how they manage large collections of objects 
            "efficiently."
                """
        for key, o_list in self.obj_lists.items():
            link = bpy.data.collections[key].objects.link
            for o in o_list:
                link(o)

    def render_screenshot(self):
        """Renders a screen capture of the current geometry."""

        bpy.context.scene.render.filepath = self.frame['image_file']
        bpy.context.scene.frame_set(1)
        bpy.ops.render.render(write_still=True)
        return

    def create_thick_filament(self, thick_id):
        """ Creates thick_fil[id] """

        # Find the thick fil
        thick_f = self.hs.thick_fil[thick_id]

        # Blender objects

        # Loop through the nodes
        crown_indices = np.arange(0, thick_f.m_no_of_cbs,
                                  thick_f.m_cbs_per_node)
        for i, ind in enumerate(crown_indices):
            crown = self.hs_b['primitives']['m_crown'].copy()
            self.b_obj['m_%i' % thick_id].append(crown)

            crown.location.x = thick_f.cb_x[i]
            crown.location.y = thick_f.m_y
            crown.location.z = thick_f.m_z
            
            print('i: %i' % i)

            crown.name = ('m_crown_%i_%i' % (thick_id,i))

    def create_primitive_m_crown(self):
        """ Creates an m_node crown"""

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thick_filament']['crown']['radius'],
            depth=self.template['thick_filament']['crown']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['thick_filament']['crown']['vertices'],
            rotation=(0, np.pi/2, 0))

        # set_active_object_shade_smooth(bpy.context.object)
        self.hs_b['primitives']['m_crown'] = bpy.context.object

    def return_all_filament_y_coords(self):
        """Returns a list of all filament y coordinates"""

        m_ys = [thick_obj.m_y for thick_obj in self.hs.thick_fil]
        a_ys = [thin_obj.a_y for thin_obj in self.hs.fhin_fil]

        return m_ys + a_ys

    def return_all_filament_z_coords(self):
        """Returns a list of all filament z coordinates"""

        m_zs = [thick_obj.m_z for thick_obj in self.hs.thick_fil]
        a_zs = [thin_obj.a_z for thin_obj in self.hs.thin.fil]

        return m_zs + a_zs

    def find_bounding_box_yz_positions(self):
        """Returns min/max yz positions of the bounding box of the geometry."""

        filament_y_positions = self.return_all_filament_y_coords()
        filament_z_positions = self.return_all_filament_z_coords()
        min_y = np.min(filament_y_positions)
        max_y = np.max(filament_y_positions)
        min_z = np.min(filament_z_positions)
        max_z = np.max(filament_z_positions)

        return ((min_y, max_y), (min_z, max_z))
