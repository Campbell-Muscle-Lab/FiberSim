# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 22:22:51 2021

@author: ken
"""

# Blender objects
import bpy
import bmesh

import math

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

        # Set up yz_scaling
        self.yz_scaling = template['lattice']['inter_thick_nm']

        # Set up
        self.setup_blender_script()
        self.setup_render_engine()

        # Create obj dictionary
        self.b_obj = dict()

        # Create materials
        self.materials = dict()
        self.create_materials()

        # self.create_lists_and_collection()

        # Create a dictionary to hold blender primitives
        self.hs_b = dict()
        self.hs_b['primitives'] = dict()

        # Create the primitives
        self.create_primitive_m_crown()
        self.create_primitive_m_stub()
        self.create_primitive_m_cat()
        self.create_primitive_m_link()
        self.create_primitive_a_node()
        self.create_primitive_a_bs()

        # Create thick filaments
        for i, thick_fil in enumerate(self.hs.thick_fil):
            self.create_thick_filament(i)

        # Create thin filaments
        for i, thin_fil in enumerate(self.hs.thin_fil):
            self.create_thin_filament(i)

        # Set nearest thin filaments
        # by adding an array into the hs.thick_fil structure
        self.set_nearest_thin_filaments()
        
        # # Create cross-bridges
        self.create_cross_bridges()

        #self.test()

        # Create a camera
        self.create_camera()

        # Create lights
        self.create_lights()

        # # Link everything
        # self.link_all_objects()
        
        # print(list(bpy.data.objects))

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

    def render_screenshot(self):
        """Renders a screen capture of the current geometry."""

        bpy.context.scene.render.filepath = self.frame['image_file']
        bpy.context.scene.frame_set(1)
        bpy.ops.render.render(write_still=True)
        return

    def create_materials(self):
        """ Creates materials """

        mat = bpy.data.materials.new(name = "m_crown")
        mat.diffuse_color = self.template['thick_filament']['crown']['color']
        self.materials['m_crown'] = mat

        mat = bpy.data.materials.new(name = "a_node")
        mat.diffuse_color = self.template['thin_filament']['node']['color']
        self.materials['a_node'] = mat

    def create_thick_filament(self, thick_id):
        """ Creates thick_fil[id] """

        print('Creating thick_fil[%i]' % thick_id)

        # Find the thick fil
        thick_f = self.hs.thick_fil[thick_id]

        # Loop through the nodes
        crown_indices = np.arange(0, thick_f.m_no_of_cbs,
                                  thick_f.m_cbs_per_node)
        for i, ind in enumerate(crown_indices):
            # crown = self.create_primitive_m_crown()
            crown = self.hs_b['primitives']['m_crown'].copy()
            crown.name = ('m_crown_%i_%i' % (thick_id, i))
            base_x = thick_f.cb_x[ind] - \
                self.template['thick_filament']['crown']['depth']/2.0
            base_y = self.yz_scaling * thick_f.m_y
            base_z = self.yz_scaling * thick_f.m_z
            distal_x = base_x + self.template['thick_filament']['crown']['depth']
            distal_y = base_y
            distal_z = base_z
            
            self.replot_primitive_cylinder(
                crown,
                np.asarray([base_x, base_y, base_z]),
                np.asarray([distal_x, distal_y, distal_z]),
                self.template['thick_filament']['crown']['depth'])
            bpy.context.collection.objects.link(crown)

    def create_thin_filament(self, thin_id):
        """ Creates thin_fil[id] """

        print('Creating thin_fil[%i]' % thin_id)

        # Set the thin file
        thin_f = self.hs.thin_fil[thin_id]

        # Loop throught the nodes
        node_indices = np.arange(0, thin_f.a_no_of_bs,
                                 thin_f.a_bs_per_node)
        
        for ind in node_indices:
            # Make the node
            node = self.hs_b['primitives']['a_node'].copy()
            node.name = ('a_node_%i_%i' % (thin_id, ind))
            base_x = thin_f.bs_x[ind] - \
                self.template['thin_filament']['node']['depth'] / 2.0
            base_y = self.yz_scaling * thin_f.a_y
            base_z = self.yz_scaling * thin_f.a_z
            distal_x = base_x + self.template['thin_filament']['node']['depth']
            distal_y = base_y
            distal_z = base_z

            self.replot_primitive_cylinder(
                node,
                np.asarray([base_x, base_y, base_z]),
                np.asarray([distal_x, distal_y, distal_z]),
                self.template['thin_filament']['node']['depth'])
            bpy.context.collection.objects.link(node)

        # Add in the binding sites
        # Loop through the bs
        bs_indices = np.arange(0, thin_f.a_no_of_bs)
        for bs_i in bs_indices:
            bs = self.hs_b['primitives']['a_bs'].copy()
            bs.name = ('a_bs_%i_%i' % (thin_id, bs_i))
            bs_angle = thin_f.bs_angle[bs_i]
            bot_x = thin_f.bs_x[bs_i]
            bot_y = self.yz_scaling * thin_f.a_y + \
                self.template['thin_filament']['node']['radius'] * \
                    np.sin(np.pi * bs_angle / 180.0)
            bot_z = self.yz_scaling * thin_f.a_z + \
                self.template['thin_filament']['node']['radius'] * \
                    np.cos(np.pi * bs_angle / 180.0)
            top_x = bot_x
            top_y = self.yz_scaling * thin_f.a_y + \
                (self.template['thin_filament']['node']['radius'] +
                  self.template['thin_filament']['bs']['depth']) * \
                    np.sin(np.pi * bs_angle / 180.0)
            top_z = self.yz_scaling * thin_f.a_z + \
                (self.template['thin_filament']['node']['radius'] +
                  self.template['thin_filament']['bs']['depth']) * \
                    np.cos(np.pi * bs_angle / 180.0)

            self.replot_primitive_cylinder(
                bs,
                np.asarray([bot_x, bot_y, bot_z]),
                np.asarray([top_x, top_y, top_z]),
                self.template['thin_filament']['bs']['depth'])
            bpy.context.collection.objects.link(bs)

    def set_nearest_thin_filaments(self):
        """ Fills an array for each hs.thick with the nearest
            thin filaments """

        # Loop through thick filaments, finding the distance to
        # each thin filament
        for thick_i, thick_f in enumerate(self.hs.thick_fil):
            hypot = np.zeros(len(self.hs.thin_fil))
            for thin_i, thin_f in enumerate(self.hs.thin_fil):
                hypot[thin_i] = np.hypot((thick_f.m_y - thin_f.a_y),
                                         (thick_f.m_z - thin_f.a_z))
            min_hypot = min(hypot)
            # Find the min distance
            self.hs.thick_fil[thick_i].nearest_a_f = -1*np.ones(6)
            counter = 0
            for i, h in enumerate(hypot):
                if (h < (1.1 * min_hypot)):
                    self.hs.thick_fil[thick_i].nearest_a_f[counter] = i
                    counter = counter + 1
        # Display
        print('thick_f: nearest_a_n values')
        for thick_f in self.hs.thick_fil:
            print(self.hs.thick_fil[thick_i].nearest_a_f)

    def create_cross_bridges(self):
        """ Draws cross-bridges """

        # Loop through myosin heads
        for thick_i, thick_f in enumerate(self.hs.thick_fil):
            for cb_i in np.arange(0, thick_f.m_no_of_cbs):
                print('Creating cross-bridges: thick_fil[%i] cb[%i]' % (thick_i, cb_i))
                stub = self.hs_b['primitives']['m_stub'].copy()
                stub.name = ('m_stub_%i_%i' % (thick_i, cb_i))
                bot_x = thick_f.cb_x[cb_i]
                bot_y = self.yz_scaling * thick_f.m_y + \
                    self.template['thick_filament']['crown']['radius'] * \
                    np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                bot_z = self.yz_scaling * thick_f.m_z + \
                    self.template['thick_filament']['crown']['radius'] * \
                    np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                top_x = bot_x
                top_y = self.yz_scaling * thick_f.m_y + \
                    (self.template['thick_filament']['crown']['radius'] +
                     self.template['thick_filament']['myosin']['stub_height']) * \
                    np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                top_z = self.yz_scaling * thick_f.m_z + \
                    (self.template['thick_filament']['crown']['radius'] +
                     self.template['thick_filament']['myosin']['stub_height']) * \
                    np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                self.replot_primitive_cylinder(
                    stub,
                    np.asarray([bot_x, bot_y, bot_z]),
                    np.asarray([top_x, top_y, top_z]),
                    self.template['thick_filament']['myosin']['stub_height'])
                bpy.context.collection.objects.link(stub)

                if (thick_f.cb_bound_to_a_f[cb_i] >= 0):
                    # Head is bound
                    # Find the cb end of the link
                    # Find the bs end of the link
                    thin_a_n = thick_f.cb_bound_to_a_f[cb_i]
                    # Check filament is a neighbor
                    if (not(any(thick_f.nearest_a_f == thin_a_n))):
                        # Continue out if not a neighbor
                        continue
                    thin_f = self.hs.thin_fil[thin_a_n]
                    thin_bs = thick_f.cb_bound_to_a_n[cb_i]
                    distal_x = thin_f.bs_x[thin_bs]
                    distal_y = self.yz_scaling * thin_f.a_y + \
                        (self.template['thin_filament']['node']['radius'] + \
                         self.template['thin_filament']['bs']['depth']) * \
                            np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    distal_z = self.yz_scaling * thin_f.a_z + \
                        (self.template['thin_filament']['node']['radius'] + \
                         self.template['thin_filament']['bs']['depth']) * \
                            np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                            
                elif (any(self.template['srx_states']==thick_f.cb_state[cb_i])):
                    # It's SRX
                    distal_x = top_x + self.template['thick_filament']['myosin']['link_height']
                    distal_y = top_y
                    distal_z = top_z
                else:
                    # It's DRX
                    distal_x = top_x
                    distal_y = top_y + \
                        self.template['thick_filament']['myosin']['link_height'] * \
                            np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                    distal_z = top_z + \
                        self.template['thick_filament']['myosin']['link_height'] * \
                            np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                # Draw link
                cb_link = self.hs_b['primitives']['m_link'].copy()
                cb_link.name = ('cb_link_%i_%i' % (thick_i, cb_i))
                self.replot_primitive_cylinder(cb_link,
                                   np.asarray([top_x, top_y, top_z]),
                                   np.asarray([distal_x, distal_y, distal_z]),
                                   self.template['thick_filament']['myosin']['link_height'])
                bpy.context.collection.objects.link(cb_link)

    def test(self):
        bot_x = 10
        bot_y = 0
        bot_z = 0
        top_x = 5
        top_y = 2
        top_z = 2
        cb_link = self.hs_b['primitives']['m_link'].copy()
        print(cb_link)
        self.replot_primitive_cylinder(cb_link,
                           np.asarray([bot_x, bot_y, bot_z]),
                           np.asarray([top_x, top_y, top_z]),
                           self.template['thick_filament']['myosin']['link_height'])
        bpy.context.collection.objects.link(cb_link)

    def replot_primitive_cylinder(self, cyl_obj, base_vec, end_vec, orig_length):
        """ Positions, scales, and orients a cylinder primitive such that
            it spans base_vec to end_vec."""

        # Calculate the new length
        new_length = np.sqrt(np.sum(np.power(base_vec - end_vec, 2)))
        length_scale = new_length / orig_length
        
        phi = np.arctan2((end_vec[1]-base_vec[1]),(end_vec[0]-base_vec[0]))
        theta = np.arccos((end_vec[2]-base_vec[2]) / new_length)

        cyl_obj.scale[2] = length_scale
        cyl_obj.rotation_euler[1] = theta
        cyl_obj.rotation_euler[2] = phi
        cyl_obj.location = 0.5*(base_vec + end_vec)

    def create_primitive_m_crown(self):
        """ Creates an m_crown """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thick_filament']['crown']['radius'],
            depth=self.template['thick_filament']['crown']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['thick_filament']['crown']['vertices'],
            rotation=(0, 0, 0))

        m_crown = bpy.context.object
        mesh = m_crown.data
        mesh.materials.append(self.materials['m_crown'])

        self.hs_b['primitives']['m_crown'] = m_crown

    def create_primitive_m_stub(self):
        """ Creates an m_stub """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thick_filament']['myosin']['stub_radius'],
            depth=self.template['thick_filament']['myosin']['stub_height'],
            enter_editmode=False,
            location=(0,0,0),
            vertices=self.template['thick_filament']['myosin']['stub_vertices'],
            rotation=(0,0,0))

        m_stub = bpy.context.object

        self.hs_b['primitives']['m_stub'] = m_stub

    def create_primitive_m_cat(self):
        """ Creates an m_cat, where myosin binds to actin """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thick_filament']['myosin']['cat_radius'],
            depth=self.template['thick_filament']['myosin']['cat_height'],
            enter_editmode=False,
            location=(0,0,0),
            vertices=self.template['thick_filament']['myosin']['cat_vertices'],
            rotation=(0,0,0))

        m_cat = bpy.context.object

        self.hs_b['primitives']['m_cat'] = m_cat

    def create_primitive_m_link(self):
        """ Creates an m_link, where myosin binds to actin """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thick_filament']['myosin']['link_radius'],
            depth=self.template['thick_filament']['myosin']['link_height'],
            enter_editmode=False,
            location=(0,0,0),
            vertices=self.template['thick_filament']['myosin']['link_vertices'],
            rotation=(0,0,0))

        m_link = bpy.context.object

        self.hs_b['primitives']['m_link'] = m_link

    def create_primitive_a_node(self):
        """ Creates an a_node """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thin_filament']['node']['radius'],
            depth=self.template['thin_filament']['node']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['thin_filament']['node']['vertices'],
            rotation=(0, 0, 0))

        a_node = bpy.context.object
        mesh = a_node.data
        mesh.materials.append(self.materials['a_node'])

        self.hs_b['primitives']['a_node'] = a_node

    def create_primitive_a_bs(self):
        """ Creates an a_bs"""

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thin_filament']['bs']['radius'],
            depth=self.template['thin_filament']['bs']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['thin_filament']['bs']['vertices'],
            rotation=(0, 0, 0))

        a_bs = bpy.context.object
        # mesh = a_bs.data
        # mesh.materials.append(self.materials['a_node'])

        self.hs_b['primitives']['a_bs'] = a_bs

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
