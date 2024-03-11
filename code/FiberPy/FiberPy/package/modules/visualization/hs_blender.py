# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 22:22:51 2021

@author: ken
"""

# Blender objects
import bpy
import bmesh

import os
import math

import numpy as np

# Half-sarcomere
import half_sarcomere as hs

class hs_blender():
    """ Class for a half-sarcomere in blender """

    def __init__(self, hs, frame, template, options, output_image_file):

        # Create local copies
        self.hs = hs
        self.frame = frame
        self.template = template
        self.options = options
        self.output_image_file = output_image_file

        # Set up yz_scaling
        self.yz_scaling = template['lattice']['inter_thick_nm']

        # Set up
        self.setup_blender_script()
        self.setup_render_engine()

        # Create obj dictionary
        self.b_obj = dict()

        self.m_colors = dict()

        # Create materials
        self.materials = dict()
        self.create_materials()

        # Create primitives and store in a dictionary
        self.hs_b = dict()
        self.hs_b['primitives'] = dict()
        self.create_primitives()

        # Create end disk
        self.create_end_disks()

        # Create thick filaments
        for i, thick_fil in enumerate(self.hs.thick_fil):
            if (thick_fil.m_y < 3):
                self.create_thick_filament(i)

        # Create thin filaments
        for i, thin_fil in enumerate(self.hs.thin_fil):
            if (thin_fil.a_y < 2):
                self.create_thin_filament(i)

        # # Set nearest thin filaments
        # # by adding an array into the hs.thick_fil structure
        # self.set_nearest_thin_filaments()

        # # # # Draws titin filaments
        # self.create_titin_filaments()

        # # Create cross-bridges
        # self.create_cross_bridges()

        # # Create mybpc
        # self.create_mybpcs()

        # Hide the primitives
        for p in self.hs_b['primitives']:
            self.hs_b['primitives'][p].hide_set(True)

        # Create a camera
        self.create_camera()

        # Create lights
        self.create_lights()

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
        bpy.context.scene.world.color = self.template['background']['color']

    def setup_render_engine(self):
        """Sets up the render engine that we want to use."""
        # Setting the render samples to 8 instead of 64 to cut render time.
        # Note: I don't think this affects much since the objects aren't very reflective.
        # bpy.context.scene.eevee.taa_render_samples = 8

        if (False):
            bpy.context.scene.render.engine = 'CYCLES'
        else:
            # Opting for Blender's workbench render engine because it does a good job and it's fast.
            bpy.context.scene.render.engine = 'BLENDER_WORKBENCH'
            bpy.context.scene.display.shading.light = 'STUDIO'
            bpy.context.scene.view_settings.view_transform = 'Filmic'
            bpy.context.scene.display.shading.show_object_outline = True
            bpy.context.scene.view_settings.exposure = 1.5
            bpy.context.scene.view_settings.gamma = 0.8
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
        bpy.context.object.data.lens = 20
        
        if ('orthographic' in self.frame['camera']):
            bpy.context.object.data.type = 'ORTHO'
            bpy.context.object.data.ortho_scale = self.frame['camera']['orthographic']['ortho_scale']
            bpy.context.object.data.clip_end = self.options['max_render_distance']

        bpy.context.scene.camera = bpy.context.object

    def create_lights(self):
        """ Create lights around the half-sarcomere """
        print("lights")

    def render_screenshot(self):
        """Renders a screen capture of the current geometry."""

        bpy.context.scene.render.filepath = self.output_image_file
        
        print(self.options)
        
        bpy.context.scene.render.resolution_x = self.options['image_pixels'][0]
        bpy.context.scene.render.resolution_y = self.options['image_pixels'][1]
                
        bpy.context.scene.frame_set(1)
        bpy.ops.render.render(write_still=True)
        return

    def create_materials(self):
        """ Creates materials """

        # Creates materials that don't depend on isotype or state
        mat = bpy.data.materials.new(name = "m_backbone")
        mat.diffuse_color = self.template['thick_filament']['backbone']['color']
        self.materials['m_backbone'] = mat

        mat = bpy.data.materials.new(name = "m_crown")
        mat.diffuse_color = self.template['thick_filament']['crown']['color']
        self.materials['m_crown'] = mat

        mat = bpy.data.materials.new(name = "a_backbone")
        mat.diffuse_color = self.template['thin_filament']['backbone']['color']
        self.materials['a_backbone'] = mat

        mat = bpy.data.materials.new(name = "a_node")
        mat.diffuse_color = self.template['thin_filament']['node']['color']
        self.materials['a_node'] = mat

        mat = bpy.data.materials.new(name = "titin")
        mat.diffuse_color = self.template['titin']['color']
        self.materials['titin'] = mat

        mat = bpy.data.materials.new(name = "z_disk")
        mat.diffuse_color = self.template['z_disk']['color']
        self.materials['z_disk'] = mat

        mat = bpy.data.materials.new(name = "m_disk")
        mat.diffuse_color = self.template['m_disk']['color']
        self.materials['m_disk'] = mat

        # Now create meshes for cb isotypes and states
        isotypes = self.template['thick_filament']['myosin']['isotypes']
        for i, iso in enumerate(isotypes):
            states = iso['states']
            for j, s in enumerate(states):
                mat_name = 'm_cb_iso_state_%i_%i' % (i,j+1)
                mat = bpy.data.materials.new(name = mat_name)
                mat.diffuse_color = s['color']
                self.materials[mat_name] = mat
                self.m_colors[mat_name] = s['color']

    def create_thick_filament(self, thick_id):
        """ Creates thick_fil[id] """

        print('Creating thick_fil[%i]' % thick_id)

        # Find the thick fil
        thick_f = self.hs.thick_fil[thick_id]

        # Draw the backbone
        backbone = self.hs_b['primitives']['m_backbone'].copy()
        backbone.name = ('m_backbone_%i' % thick_id)
        proximal_x = self.hs.hs_length
        proximal_y = self.yz_scaling * thick_f.m_y
        proximal_z = self.yz_scaling * thick_f.m_z
        distal_x = thick_f.cb_x[-1]
        distal_y = proximal_y
        distal_z = proximal_z

        self.replot_primitive_cylinder(
                backbone,
                np.asarray([proximal_x, proximal_y, proximal_z]),
                np.asarray([distal_x, distal_y, distal_z]),
                self.template['thick_filament']['backbone']['depth'])

        bpy.context.collection.objects.link(backbone)

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
        
        # Store camera location
        cam_loc = np.asarray([self.frame['camera']['location']['x'],
                              self.frame['camera']['location']['y'],
                              self.frame['camera']['location']['z']])
        
        # Branch on mode
        if (self.frame['style'] == 'x_ray'):
            
            thin_f = self.hs.thin_fil[thin_id]
            bs_indices = np.arange(0, thin_f.a_no_of_bs)
            
            for bs_i in bs_indices:
                
                bs_angle = thin_f.bs_angle[bs_i]
                
                sph = self.hs_b['primitives']['bs_sphere'].copy()
                sph.name = ('a_bs_sp_%i_%i' % (thin_id, bs_i))
                
                bot_x = thin_f.bs_x[bs_i]
                bot_y = self.yz_scaling * thin_f.a_y + \
                    self.template['thin_filament']['node']['radius'] * \
                        np.sin(np.pi * bs_angle / 180.0)
                bot_z = self.yz_scaling * thin_f.a_z + \
                    self.template['thin_filament']['node']['radius'] * \
                        np.cos(np.pi * bs_angle / 180.0)

                # Check distance
                h = np.linalg.norm(cam_loc - np.asarray([bot_x, bot_y, bot_z]))
                if (h > self.options['max_render_distance']):
                    continue
                
                sph.location=(bot_x, bot_y, bot_z)
                
                bpy.context.collection.objects.link(sph)
            
            return

        # Set the thin file
        thin_f = self.hs.thin_fil[thin_id]

        # Draw the backbone
        backbone = self.hs_b['primitives']['a_backbone'].copy()
        backbone.name = ('a_backbone_%i' % thin_id)
        proximal_x = 0
        proximal_y = self.yz_scaling * thin_f.a_y
        proximal_z = self.yz_scaling * thin_f.a_z
        distal_x = thin_f.bs_x[-1]
        distal_y = proximal_y
        distal_z = proximal_z

        self.replot_primitive_cylinder(
                backbone,
                np.asarray([proximal_x, proximal_y, proximal_z]),
                np.asarray([distal_x, distal_y, distal_z]),
                self.template['thin_filament']['backbone']['depth'])

        bpy.context.collection.objects.link(backbone)

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
            s = thin_f.bs_state[bs_i]
            bs_mesh = 'a_state_%i' % (s)
            bs = self.hs_b['primitives'][bs_mesh].copy()

            bs.name = ('a_bs_%i_%i' % (thin_id, bs_i))
            bs_angle = thin_f.bs_angle[bs_i]
            bot_x = thin_f.bs_x[bs_i]
            bot_y = self.yz_scaling * thin_f.a_y + \
                self.template['thin_filament']['node']['radius'] * \
                    np.sin(np.pi * bs_angle / 180.0)
            bot_z = self.yz_scaling * thin_f.a_z + \
                self.template['thin_filament']['node']['radius'] * \
                    np.cos(np.pi * bs_angle / 180.0)

            # Check distance
            h = np.linalg.norm(cam_loc - np.asarray([bot_x, bot_y, bot_z]))
            if (h > self.options['max_render_distance']):
                continue

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

    def create_titin_filaments(self):
        """ Draws titin filaments """

        # Loop through myosin filaments
        for thick_i, thick_f in enumerate(self.hs.thick_fil):
            for a_f in thick_f.nearest_a_f:
                # Convert the a_f to an integer
                a_f = int(a_f)
                # Check filament is a neighbor
                if (a_f == -1):
                    # Continue out if not a neighbor
                    continue

                # Draw titin
                # Copy the titin primitive and set its name
                tit = self.hs_b['primitives']['titin'].copy()
                tit.name = ('titin_%i_%i' % (thick_i, a_f))

                # Work out the start and points
                cb_ind = (6 * self.hs.t_attach_m_node) - 1
                thick_x = thick_f.cb_x[cb_ind]
                thick_y = self.yz_scaling * thick_f.m_y
                thick_z = self.yz_scaling * thick_f.m_z
                
                bs_ind = (2 * self.hs.t_attach_a_node) - 1
                thin_x = self.hs.thin_fil[a_f].bs_x[bs_ind]
                thin_y = self.yz_scaling * self.hs.thin_fil[a_f].a_y
                thin_z = self.yz_scaling * self.hs.thin_fil[a_f].a_z

                # Draw and link
                self.replot_primitive_cylinder(
                    tit,
                    np.asarray([thick_x, thick_y, thick_z]),
                    np.asarray([thin_x, thin_y, thin_z]),
                    self.template['titin']['depth'])

                bpy.context.collection.objects.link(tit)

    def create_cross_bridges(self):
        """ Draws cross-bridges """

        # Store camera location
        cam_loc = np.asarray([self.frame['camera']['location']['x'],
                              self.frame['camera']['location']['y'],
                              self.frame['camera']['location']['z']])

        # Loop through myosin filaments
        for thick_i, thick_f in enumerate(self.hs.thick_fil):

            # Loop through myosin heads
            for cb_i in np.arange(0, thick_f.m_no_of_cbs):
                if ((cb_i % 50)==0):
                    print('Creating cross-bridges: thick_fil[%i] cb[%i]' % (thick_i, cb_i))

                # Get the isotype and state
                cb_isotype = thick_f.cb_iso[cb_i] 
                cb_state = thick_f.cb_state[cb_i] 
                
                isotype_data = self.template['thick_filament']['myosin']['isotypes'][cb_isotype - 1]
                state_data = isotype_data['states'][cb_state -1]
                state_type = state_data['type']

                # Generate names - these change depending on the isotype and stub
                prim_names = dict()
                for n in ['stub', 'link', 'extension', 'lever', 'cat']:
                    prim_names[n] = 'm_%s_iso_%i_state_%i' % (n, cb_isotype, cb_state)

                # Prepare to draw the stub
                bot_x = thick_f.cb_x[cb_i]
                bot_y = self.yz_scaling * thick_f.m_y + \
                    self.template['thick_filament']['crown']['radius'] * \
                    np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                bot_z = self.yz_scaling * thick_f.m_z + \
                    self.template['thick_filament']['crown']['radius'] * \
                    np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                # Check distance
                draw_sphere = False
                h = np.linalg.norm(cam_loc - np.asarray([bot_x, bot_y, bot_z]))
                if (h > self.options['max_render_distance']):
                    continue
                if (h < self.options['max_smooth_distance']):
                    draw_sphere = True

                top_x = bot_x
                top_y = self.yz_scaling * thick_f.m_y + \
                    (self.template['thick_filament']['crown']['radius'] +
                     self.template['thick_filament']['myosin']['stub_height']) * \
                    np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                top_z = self.yz_scaling * thick_f.m_z + \
                    (self.template['thick_filament']['crown']['radius'] +
                     self.template['thick_filament']['myosin']['stub_height']) * \
                    np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                # Copy the stub primitive and set its name
                stub = self.hs_b['primitives'][prim_names['stub']].copy()
                stub.name = ('m_stub_%i_%i' % (thick_i, cb_i))

                # Draw and link
                self.replot_primitive_cylinder(
                    stub,
                    np.asarray([bot_x, bot_y, bot_z]),
                    np.asarray([top_x, top_y, top_z]),
                    self.template['thick_filament']['myosin']['stub_height'])

                bpy.context.collection.objects.link(stub)

                if (draw_sphere):
                    # Draw sphere and link
                    sp_name = 'm_sphere_stub_dist_iso_%i_state_%i' % \
                        (cb_isotype, cb_state)
                    stub_top = self.hs_b['primitives'][sp_name].copy()
                    stub_top.name = ('m_stub_top_%i_%i' % (thick_i, cb_i))
                    stub_top.location=(top_x, top_y, top_z)

                    bpy.context.collection.objects.link(stub_top)

                    sp_name = 'm_sphere_stub_prox_iso_%i_state_%i' % \
                        (cb_isotype, cb_state)
                    stub_bot = self.hs_b['primitives'][sp_name].copy()
                    stub_bot.name = ('m_stub_bot_%i_%i' % (thick_i, cb_i))
                    stub_bot.location=(bot_x, bot_y, bot_z)

                    bpy.context.collection.objects.link(stub_bot)

                # Now move on to the remainder
                if (thick_f.cb_bound_to_a_f[cb_i] >= 0):
                    # Head is bound
                    thin_a_n = thick_f.cb_bound_to_a_f[cb_i]

                    # Check filament is a neighbor
                    if (not(any(thick_f.nearest_a_f == thin_a_n))):
                        # Continue out if not a neighbor
                        continue

                    # Deduce the thin filament
                    thin_f = self.hs.thin_fil[thin_a_n]
                    thin_bs = thick_f.cb_bound_to_a_n[cb_i]

                    # Deduce the cat domain coordinates
                    cat_distal_x = thin_f.bs_x[thin_bs]
                    cat_distal_y = self.yz_scaling * thin_f.a_y + \
                        (self.template['thin_filament']['node']['radius'] +
                         self.template['thin_filament']['bs']['depth']) * \
                             np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    cat_distal_z = self.yz_scaling * thin_f.a_z + \
                        (self.template['thin_filament']['node']['radius'] +
                         self.template['thin_filament']['bs']['depth']) * \
                             np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    cat_proximal_x = cat_distal_x
                    cat_proximal_y = cat_distal_y + \
                         self.template['thick_filament']['myosin']['cat_height'] * \
                             np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    cat_proximal_z = cat_distal_z + \
                        self.template['thick_filament']['myosin']['cat_height'] * \
                             np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)

                    # Copy the cat primitive and set its name
                    cat = self.hs_b['primitives'][prim_names['cat']].copy()
                    cat.name = ('m_cat_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        cat,
                        np.asarray([cat_distal_x, cat_distal_y, cat_distal_z]),
                        np.asarray([cat_proximal_x, cat_proximal_y, cat_proximal_z]),
                        self.template['thick_filament']['myosin']['cat_height'])

                    bpy.context.collection.objects.link(cat)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_cat_prox_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        cat_prox = self.hs_b['primitives'][sp_name].copy()
                        cat_prox.name = ('m_stub_top_%i_%i' % (thick_i, cb_i))
                        cat_prox.location=(cat_proximal_x, cat_proximal_y, cat_proximal_z)
    
                        bpy.context.collection.objects.link(cat_prox)

                    # Lever arm coordinates
                    lever_proximal_x = cat_distal_x - state_data['extension']
                    lever_proximal_y = cat_distal_y + \
                         self.template['thick_filament']['myosin']['lever_height'] * \
                             np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    lever_proximal_z = cat_distal_z + \
                        self.template['thick_filament']['myosin']['lever_height'] * \
                             np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)

                    # Copy the lever primitive and set its name
                    lever = self.hs_b['primitives'][prim_names['lever']].copy()
                    lever.name = ('m_lever_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        lever,
                        np.asarray([cat_proximal_x, cat_proximal_y, cat_proximal_z]),
                        np.asarray([lever_proximal_x, lever_proximal_y, lever_proximal_z]),
                        self.template['thick_filament']['myosin']['lever_height'])

                    bpy.context.collection.objects.link(lever)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_lever_prox_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        lever_prox = self.hs_b['primitives'][sp_name].copy()
                        lever_prox.name = ('m_lever_prox_%i_%i' % (thick_i, cb_i))
                        lever_prox.location=(lever_proximal_x, lever_proximal_y, lever_proximal_z)
    
                        bpy.context.collection.objects.link(lever_prox)

                    # Link distal coordinates
                    link_distal_x = top_x
                    link_distal_y = lever_proximal_y
                    link_distal_z = lever_proximal_z

                    # Copy the extension primitive and set its name
                    ext = self.hs_b['primitives'][prim_names['extension']].copy()
                    ext.name = ('m_extension_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        ext,
                        np.asarray([lever_proximal_x, lever_proximal_y, lever_proximal_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['myosin']['extension_height'])

                    bpy.context.collection.objects.link(ext)

                    # Copy the link primitive and set its name
                    link = self.hs_b['primitives'][prim_names['link']].copy()
                    link.name = ('m_link_%i_%i' % (thick_i, cb_i))

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_link_dist_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        link_dist = self.hs_b['primitives'][sp_name].copy()
                        link_dist.name = ('m_link_dist_%i_%i' % (thick_i, cb_i))
                        link_dist.location=(link_distal_x, link_distal_y, link_distal_z)
    
                        bpy.context.collection.objects.link(link_dist)

                    # Draw and link
                    self.replot_primitive_cylinder(
                        link,
                        np.asarray([top_x, top_y, top_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['myosin']['link_height'])

                    bpy.context.collection.objects.link(link)

                elif (state_type == 'S'):
                    # It's SRX
                    
                    link_distal_x = top_x + self.template['thick_filament']['myosin']['link_height']
                    link_distal_y = top_y
                    link_distal_z = top_z

                    # Copy the link primitive and set its name
                    link = self.hs_b['primitives'][prim_names['link']].copy()
                    link.name = ('m_link_%i_%i' % (thick_i, cb_i))

                    bpy.context.collection.objects.link(link)

                    # Draw and link
                    self.replot_primitive_cylinder(
                        link,
                        np.asarray([top_x, top_y, top_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['myosin']['link_height'])

                    # Lever arm coordinates
                    lever_distal_x = link_distal_x + self.template['thick_filament']['myosin']['lever_height']
                    lever_distal_y = link_distal_y
                    lever_distal_z = link_distal_z

                    # Copy the lever primitive and set its name
                    lever = self.hs_b['primitives'][prim_names['lever']].copy()
                    lever.name = ('m_lever_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        lever,
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        np.asarray([lever_distal_x, lever_distal_y, lever_distal_z]),
                        self.template['thick_filament']['myosin']['lever_height'])

                    bpy.context.collection.objects.link(lever)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_link_dist_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        link_dist = self.hs_b['primitives'][sp_name].copy()
                        link_dist.name = ('m_link_dist_%i_%i' % (thick_i, cb_i))
                        link_dist.location=(link_distal_x, link_distal_y, link_distal_z)
    
                        bpy.context.collection.objects.link(link_dist)

                    # Lever arm coordinates
                    cat_distal_x = lever_distal_x + self.template['thick_filament']['myosin']['cat_height']
                    cat_distal_y = link_distal_y
                    cat_distal_z = link_distal_z

                     # Copy the cat primitive and set its name
                    cat = self.hs_b['primitives'][prim_names['cat']].copy()
                    cat.name = ('m_cat_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        cat,
                        np.asarray([lever_distal_x, lever_distal_y, lever_distal_z]),
                        np.asarray([cat_distal_x, cat_distal_y, cat_distal_z]),
                        self.template['thick_filament']['myosin']['cat_height'])

                    bpy.context.collection.objects.link(cat)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_cat_prox_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        cat_prox = self.hs_b['primitives'][sp_name].copy()
                        cat_prox.name = ('m_stub_top_%i_%i' % (thick_i, cb_i))
                        cat_prox.location=(lever_distal_x, lever_distal_y, lever_distal_z)
    
                        bpy.context.collection.objects.link(cat_prox)

                else:
                    # It's DRX
                    
                    link_distal_x = top_x
                    link_distal_y = top_y + \
                        self.template['thick_filament']['myosin']['link_height'] * \
                            np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                    link_distal_z = top_z + \
                        self.template['thick_filament']['myosin']['link_height'] * \
                            np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                    # Copy the link primitive and set its name
                    link = self.hs_b['primitives'][prim_names['link']].copy()
                    link.name = ('m_link_%i_%i' % (thick_i, cb_i))

                    bpy.context.collection.objects.link(link)

                    # Draw and link
                    self.replot_primitive_cylinder(
                        link,
                        np.asarray([top_x, top_y, top_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['myosin']['link_height'])

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_link_dist_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        link_dist = self.hs_b['primitives'][sp_name].copy()
                        link_dist.name = ('m_link_dist_%i_%i' % (thick_i, cb_i))
                        link_dist.location=(link_distal_x, link_distal_y, link_distal_z)
    
                        bpy.context.collection.objects.link(link_dist)

                    # Lever arm coordinates
                    lever_distal_x = link_distal_x
                    lever_distal_y = link_distal_y + \
                        self.template['thick_filament']['myosin']['lever_height'] * \
                            np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                    lever_distal_z = link_distal_z  + \
                        self.template['thick_filament']['myosin']['lever_height'] * \
                            np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                    # Copy the lever primitive and set its name
                    lever = self.hs_b['primitives'][prim_names['lever']].copy()
                    lever.name = ('m_lever_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        lever,
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        np.asarray([lever_distal_x, lever_distal_y, lever_distal_z]),
                        self.template['thick_filament']['myosin']['lever_height'])

                    bpy.context.collection.objects.link(lever)

                    # Cat arm coordinates
                    cat_distal_x = lever_distal_x
                    cat_distal_y = lever_distal_y + \
                        self.template['thick_filament']['myosin']['cat_height'] * \
                            np.sin(np.pi * thick_f.cb_angle[cb_i] / 180.0)
                    cat_distal_z = lever_distal_z + \
                        self.template['thick_filament']['myosin']['cat_height'] * \
                            np.cos(np.pi * thick_f.cb_angle[cb_i] / 180.0)

                     # Copy the cat primitive and set its name
                    cat = self.hs_b['primitives'][prim_names['cat']].copy()
                    cat.name = ('m_cat_%i_%i' % (thick_i, cb_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        cat,
                        np.asarray([lever_distal_x, lever_distal_y, lever_distal_z]),
                        np.asarray([cat_distal_x, cat_distal_y, cat_distal_z]),
                        self.template['thick_filament']['myosin']['cat_height'])

                    bpy.context.collection.objects.link(cat)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'm_sphere_cat_prox_iso_%i_state_%i' % \
                            (cb_isotype, cb_state)
                        cat_prox = self.hs_b['primitives'][sp_name].copy()
                        cat_prox.name = ('m_stub_top_%i_%i' % (thick_i, cb_i))
                        cat_prox.location=(lever_distal_x, lever_distal_y, lever_distal_z)
    
                        bpy.context.collection.objects.link(cat_prox)


    def create_mybpcs(self):
        """ Draws myosin binding protein c molecules """

        # Store camera location
        cam_loc = np.asarray([self.frame['camera']['location']['x'],
                              self.frame['camera']['location']['y'],
                              self.frame['camera']['location']['z']])

        # Loop through myosin filaments
        for thick_i, thick_f in enumerate(self.hs.thick_fil):
            print('Creating mybpcs: thick_fil[%i]' % thick_i)

            # Loop through mybpc
            for pc_i in np.arange(0, thick_f.c_no_of_pcs):

                # Get the isotype and state
                pc_isotype = thick_f.pc_iso[pc_i]
                pc_state = thick_f.pc_state[pc_i]
                
                isotype_data = self.template['thick_filament']['mybpc']['isotypes'][pc_isotype-1]
                state_data = isotype_data['states'][pc_state-1]
                state_type = state_data['type']

                # Generate names - these change depending on the isotype and stub
                prim_names = dict()
                for n in ['stub', 'link', 'extension', 'lever', 'cat']:
                    prim_names[n] = 'c_%s_iso_%i_state_%i' % (n, pc_isotype, pc_state)

                # Always draw the stub
                bot_x = thick_f.pc_x[pc_i]
                bot_y = self.yz_scaling * thick_f.m_y + \
                    self.template['thick_filament']['crown']['radius'] * \
                    np.sin(np.pi * thick_f.pc_angle[pc_i] / 180.0)
                bot_z = self.yz_scaling * thick_f.m_z + \
                    self.template['thick_filament']['crown']['radius'] * \
                    np.cos(np.pi * thick_f.pc_angle[pc_i] / 180.0)

                # Check distance
                draw_sphere = False
                h = np.linalg.norm(cam_loc - np.asarray([bot_x, bot_y, bot_z]))
                if (h > self.options['max_render_distance']):
                    continue
                if (h < self.options['max_smooth_distance']):
                    draw_sphere = True

                top_x = bot_x
                top_y = self.yz_scaling * thick_f.m_y + \
                    (self.template['thick_filament']['crown']['radius'] +
                     self.template['thick_filament']['mybpc']['stub_height']) * \
                    np.sin(np.pi * thick_f.pc_angle[pc_i] / 180.0)
                top_z = self.yz_scaling * thick_f.m_z + \
                    (self.template['thick_filament']['crown']['radius'] +
                     self.template['thick_filament']['mybpc']['stub_height']) * \
                    np.cos(np.pi * thick_f.pc_angle[pc_i] / 180.0)

                # Copy the stub primitive and set its name
                stub = self.hs_b['primitives'][prim_names['stub']].copy()
                stub.name = ('c_stub_%i_%i' % (thick_i, pc_i))

                # Draw and link
                self.replot_primitive_cylinder(
                    stub,
                    np.asarray([bot_x, bot_y, bot_z]),
                    np.asarray([top_x, top_y, top_z]),
                    self.template['thick_filament']['mybpc']['stub_height'])

                bpy.context.collection.objects.link(stub)

                if (draw_sphere):
                    # Draw sphere and link
                    sp_name = 'c_sphere_stub_dist_iso_%i_state_%i' % \
                        (pc_isotype, pc_state)
                    stub_top = self.hs_b['primitives'][sp_name].copy()
                    stub_top.name = ('c_stub_top_%i_%i' % (thick_i, pc_i))
                    stub_top.location=(top_x, top_y, top_z)
    
                    bpy.context.collection.objects.link(stub_top)
    
                    sp_name = 'c_sphere_stub_prox_iso_%i_state_%i' % \
                        (pc_isotype, pc_state)
                    stub_bot = self.hs_b['primitives'][sp_name].copy()
                    stub_bot.name = ('c_stub_bot_%i_%i' % (thick_i, pc_i))
                    stub_bot.location=(bot_x, bot_y, bot_z)
    
                    bpy.context.collection.objects.link(stub_bot)


                # Now move on to the remainder
                if (thick_f.pc_bound_to_a_f[pc_i] >= 0):
                    # MyBPC is bound
                    thin_a_n = thick_f.pc_bound_to_a_f[pc_i]

                    # Check filament is a neighbor
                    if (not(any(thick_f.nearest_a_f == thin_a_n))):
                        # Continue out if not a neighbor
                        continue

                    # Deduce the thin filament
                    thin_f = self.hs.thin_fil[thin_a_n]
                    thin_bs = thick_f.pc_bound_to_a_n[pc_i]

                    # Deduce the cat domain coordinates
                    cat_distal_x = thin_f.bs_x[thin_bs]
                    cat_distal_y = self.yz_scaling * thin_f.a_y + \
                        (self.template['thin_filament']['node']['radius'] +
                         self.template['thin_filament']['bs']['depth']) * \
                             np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    cat_distal_z = self.yz_scaling * thin_f.a_z + \
                        (self.template['thin_filament']['node']['radius'] +
                         self.template['thin_filament']['bs']['depth']) * \
                             np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    cat_proximal_x = cat_distal_x
                    cat_proximal_y = cat_distal_y + \
                         self.template['thick_filament']['mybpc']['cat_height'] * \
                             np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    cat_proximal_z = cat_distal_z + \
                        self.template['thick_filament']['mybpc']['cat_height'] * \
                             np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)

                    # Copy the cat primitive and set its name
                    cat = self.hs_b['primitives'][prim_names['cat']].copy()
                    cat.name = ('c_cat_%i_%i' % (thick_i, pc_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        cat,
                        np.asarray([cat_distal_x, cat_distal_y, cat_distal_z]),
                        np.asarray([cat_proximal_x, cat_proximal_y, cat_proximal_z]),
                        self.template['thick_filament']['mybpc']['cat_height'])

                    bpy.context.collection.objects.link(cat)

                    # Lever arm coordinates
                    lever_proximal_x = cat_distal_x
                    lever_proximal_y = cat_distal_y + \
                         self.template['thick_filament']['mybpc']['lever_height'] * \
                             np.sin(np.pi * thin_f.bs_angle[thin_bs] / 180.0)
                    lever_proximal_z = cat_distal_z + \
                        self.template['thick_filament']['mybpc']['lever_height'] * \
                             np.cos(np.pi * thin_f.bs_angle[thin_bs] / 180.0)

                    # Copy the lever primitive and set its name
                    lever = self.hs_b['primitives'][prim_names['lever']].copy()
                    lever.name = ('c_lever_%i_%i' % (thick_i, pc_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        lever,
                        np.asarray([cat_proximal_x, cat_proximal_y, cat_proximal_z]),
                        np.asarray([lever_proximal_x, lever_proximal_y, lever_proximal_z]),
                        self.template['thick_filament']['mybpc']['lever_height'])

                    bpy.context.collection.objects.link(lever)

                    # Link distal coordinates
                    link_distal_x = top_x
                    link_distal_y = lever_proximal_y
                    link_distal_z = lever_proximal_z

                    # Copy the extension primitive and set its name
                    ext = self.hs_b['primitives'][prim_names['extension']].copy()
                    ext.name = ('c_extension_%i_%i' % (thick_i, pc_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        ext,
                        np.asarray([lever_proximal_x, lever_proximal_y, lever_proximal_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['mybpc']['extension_height'])

                    bpy.context.collection.objects.link(ext)

                    # Copy the link primitive and set its name
                    link = self.hs_b['primitives'][prim_names['link']].copy()
                    link.name = ('c_link_%i_%i' % (thick_i, pc_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        link,
                        np.asarray([top_x, top_y, top_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['mybpc']['link_height'])

                    bpy.context.collection.objects.link(link)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'c_sphere_link_dist_iso_%i_state_%i' % \
                            (pc_isotype, pc_state)
                        link_dist = self.hs_b['primitives'][sp_name].copy()
                        link_dist.name = ('c_link_dist_%i_%i' % (thick_i, pc_i))
                        link_dist.location=(link_distal_x, link_distal_y, link_distal_z)
    
                        bpy.context.collection.objects.link(link_dist)

                else:
                    # It's not attached
                    link_distal_x = top_x
                    link_distal_y = top_y + \
                        self.template['thick_filament']['mybpc']['link_height'] * \
                            np.sin(np.pi * thick_f.pc_angle[pc_i] / 180.0)
                    link_distal_z = top_z + \
                        self.template['thick_filament']['mybpc']['link_height'] * \
                            np.cos(np.pi * thick_f.pc_angle[pc_i] / 180.0)

                    # Copy the link primitive and set its name
                    link = self.hs_b['primitives'][prim_names['link']].copy()
                    link.name = ('c_link_%i_%i' % (thick_i, pc_i))

                    bpy.context.collection.objects.link(link)

                    # Draw and link
                    self.replot_primitive_cylinder(
                        link,
                        np.asarray([top_x, top_y, top_z]),
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        self.template['thick_filament']['mybpc']['link_height'])

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'c_sphere_link_dist_iso_%i_state_%i' % \
                            (pc_isotype, pc_state)
                        link_dist = self.hs_b['primitives'][sp_name].copy()
                        link_dist.name = ('c_link_dist_%i_%i' % (thick_i, pc_i))
                        link_dist.location=(link_distal_x, link_distal_y, link_distal_z)
    
                        bpy.context.collection.objects.link(link_dist)

                    # Lever arm coordinates
                    lever_distal_x = link_distal_x
                    lever_distal_y = link_distal_y + \
                        self.template['thick_filament']['mybpc']['lever_height'] * \
                            np.sin(np.pi * thick_f.pc_angle[pc_i] / 180.0)
                    lever_distal_z = link_distal_z  + \
                        self.template['thick_filament']['mybpc']['lever_height'] * \
                            np.cos(np.pi * thick_f.pc_angle[pc_i] / 180.0)

                    # Copy the lever primitive and set its name
                    lever = self.hs_b['primitives'][prim_names['lever']].copy()
                    lever.name = ('c_lever_%i_%i' % (thick_i, pc_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        lever,
                        np.asarray([link_distal_x, link_distal_y, link_distal_z]),
                        np.asarray([lever_distal_x, lever_distal_y, lever_distal_z]),
                        self.template['thick_filament']['mybpc']['lever_height'])

                    bpy.context.collection.objects.link(lever)

                    # Cat arm coordinates
                    cat_distal_x = lever_distal_x
                    cat_distal_y = lever_distal_y + \
                        self.template['thick_filament']['mybpc']['cat_height'] * \
                            np.sin(np.pi * thick_f.pc_angle[pc_i] / 180.0)
                    cat_distal_z = lever_distal_z + \
                        self.template['thick_filament']['mybpc']['cat_height'] * \
                            np.cos(np.pi * thick_f.pc_angle[pc_i] / 180.0)

                     # Copy the cat primitive and set its name
                    cat = self.hs_b['primitives'][prim_names['cat']].copy()
                    cat.name = ('c_cat_%i_%i' % (thick_i, pc_i))

                    # Draw and link
                    self.replot_primitive_cylinder(
                        cat,
                        np.asarray([lever_distal_x, lever_distal_y, lever_distal_z]),
                        np.asarray([cat_distal_x, cat_distal_y, cat_distal_z]),
                        self.template['thick_filament']['myosin']['cat_height'])

                    bpy.context.collection.objects.link(cat)

                    if (draw_sphere):
                        # Draw sphere and link
                        sp_name = 'c_sphere_cat_prox_iso_%i_state_%i' % \
                            (pc_isotype, pc_state)
                        cat_prox = self.hs_b['primitives'][sp_name].copy()
                        cat_prox.name = ('m_stub_top_%i_%i' % (thick_i, pc_i))
                        cat_prox.location=(lever_distal_x, lever_distal_y, lever_distal_z)
    
                        bpy.context.collection.objects.link(cat_prox)


    def create_end_disks(self):
        """ Creates the z and m disks """

        # Get some dimensions
        y,z = self.find_bounding_box_yz_positions()
        y_min = y[0]
        y_max = y[1]
        z_min = z[0]
        z_max = z[1]

        center_y = self.yz_scaling * (y_min + y_max) / 2.0
        center_z = self.yz_scaling * (z_min + z_max) / 2.0

        # Make the z-disk
        z = self.hs_b['primitives']['z_disk'].copy()
        z.scale[1] = self.yz_scaling * (y_max - y_min) * 1.1
        z.scale[2] = self.yz_scaling * (z_max - z_min) * 1.1
        z.location = (-0.5, center_y, center_z)

        bpy.context.collection.objects.link(z)

        # Make the m-disk
        m = self.hs_b['primitives']['m_disk'].copy()
        m.scale[1] = self.yz_scaling * (y_max - y_min) * 1.1
        m.scale[2] = self.yz_scaling * (z_max - z_min) * 1.1
        m.location = (self.hs.hs_length + 0.5, center_y, center_z)

        bpy.context.collection.objects.link(m)

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

    def create_primitives(self):
        """ Creates all primitives """

        self.create_primitive_m_backbone()
        self.create_primitive_m_crown()
        self.create_cb_primitives()
        self.create_pc_primitives()
        self.create_primitive_a_backbone()
        self.create_primitive_a_node()
        self.create_bs_primitives()
        self.create_primitive_titin()
        self.create_primitive_end_disks()

    def create_primitive_m_backbone(self):
        """ Creates an m_backbone """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thick_filament']['backbone']['radius'],
            depth=self.template['thick_filament']['backbone']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['thick_filament']['backbone']['vertices'],
            rotation=(0, 0, 0))

        m_backbone = bpy.context.object
        mesh = m_backbone.data
        mesh.materials.append(self.materials['m_backbone'])

        self.hs_b['primitives']['m_backbone'] = m_backbone

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
    
    def create_cb_primitives(self):
        """ Cycle though the cb isotypes and states, making a mesh for the
            components of each one. This seems to be the easiest way of
            assigning separate colors to each state """
        
        # Loop through isotypes
        isotypes = self.template['thick_filament']['myosin']['isotypes']
        for i, iso in enumerate(isotypes):
            states = iso['states']
            for j, s in enumerate(states):
                for f in ['stub','link','extension','lever','cat']:
                    bpy.ops.mesh.primitive_cylinder_add(
                        radius=self.template['thick_filament']['myosin']\
                            ['%s_radius' % f],
                        depth=self.template['thick_filament']['myosin']\
                            ['%s_height' % f],
                        enter_editmode=False,
                        location=(0,0,0),
                        vertices=self.template['thick_filament']['myosin']\
                            ['%s_vertices' % f],
                        rotation=(0,0,0))

                    m = bpy.context.object
                    mesh = m.data
                    m_name = 'm_%s_iso_%i_state_%i' % (f, i+1,j+1)
                    mat = bpy.data.materials.new(name=m_name)
                    mat.diffuse_color = s['color']
                    mesh.materials.append(mat)

                    self.hs_b['primitives'][m_name] = m
                    
            # Add in spheres to 'round' corners
            for j, s in enumerate(states):
                for f in ['stub_prox','stub_dist','link_dist','lever_prox','cat_prox']:
                    if ('stub' in f):
                        r = self.template['thick_filament']['myosin']['stub_radius']
                        v = self.template['thick_filament']['myosin']['stub_vertices']
                    elif ('link' in f):
                        r = self.template['thick_filament']['myosin']['link_radius']
                        v = self.template['thick_filament']['myosin']['link_vertices']
                    elif ('lever' in f):
                        r = self.template['thick_filament']['myosin']['lever_radius']
                        v = self.template['thick_filament']['myosin']['lever_vertices']
                    else:
                        r = self.template['thick_filament']['myosin']['cat_radius']
                        v = self.template['thick_filament']['myosin']['cat_vertices']

                    bpy.ops.mesh.primitive_ico_sphere_add(
                        radius=r,
                        enter_editmode=False,
                        subdivisions=4,
                        location=(0,0,0))
                    
                    m = bpy.context.object
                    mesh = m.data
                    m_name = 'm_sphere_%s_iso_%i_state_%i' % (f, i+1, j+1)
                    mat = bpy.data.materials.new(name=m_name)
                    mat.diffuse_color = s['color']
                    mesh.materials.append(mat)
                    
                    self.hs_b['primitives'][m_name] = m

    def create_pc_primitives(self):
        """ Cycle though the mybpc isotypes and states, making a mesh for the
            components of each one. This seems to be the easiest way of
            assigning separate colors to each state """
        
        # Loop through isotypes
        isotypes = self.template['thick_filament']['mybpc']['isotypes']
        for i, iso in enumerate(isotypes):
            states = iso['states']
            for j, s in enumerate(states):
                for f in ['stub','link','extension','lever', 'cat']:
                    bpy.ops.mesh.primitive_cylinder_add(
                        radius=self.template['thick_filament']['mybpc']\
                            ['%s_radius' % f],
                        depth=self.template['thick_filament']['mybpc']\
                            ['%s_height' % f],
                        enter_editmode=False,
                        location=(0,0,0),
                        vertices=self.template['thick_filament']['mybpc']\
                            ['%s_vertices' % f],
                        rotation=(0,0,0))

                    c = bpy.context.object
                    mesh = c.data
                    c_name = 'c_%s_iso_%i_state_%i' % (f, i+1,j+1)
                    mat = bpy.data.materials.new(name=c_name)
                    mat.diffuse_color = s['color']
                    mesh.materials.append(mat)

                    self.hs_b['primitives'][c_name] = c

            # Add in spheres to 'round' corners
            for j, s in enumerate(states):
                for f in ['stub_prox','stub_dist','link_dist','lever_prox','cat_prox']:
                    if ('stub' in f):
                        r = self.template['thick_filament']['mybpc']['stub_radius']
                        v = self.template['thick_filament']['mybpc']['stub_vertices']
                    elif ('link' in f):
                        r = self.template['thick_filament']['mybpc']['link_radius']
                        v = self.template['thick_filament']['mybpc']['link_vertices']
                    elif ('lever' in f):
                        r = self.template['thick_filament']['mybpc']['lever_radius']
                        v = self.template['thick_filament']['mybpc']['lever_vertices']
                    else:
                        r = self.template['thick_filament']['mybpc']['cat_radius']
                        v = self.template['thick_filament']['mybpc']['cat_vertices']

                    bpy.ops.mesh.primitive_ico_sphere_add(
                        radius=r,
                        enter_editmode=False,
                        subdivisions=4,
                        location=(0,0,0))
                    
                    m = bpy.context.object
                    mesh = m.data
                    m_name = 'c_sphere_%s_iso_%i_state_%i' % (f, i+1, j+1)
                    mat = bpy.data.materials.new(name=m_name)
                    mat.diffuse_color = s['color']
                    mesh.materials.append(mat)
                    
                    self.hs_b['primitives'][m_name] = m


    def create_primitive_a_backbone(self):
        """ Creates an a_backbone """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['thin_filament']['backbone']['radius'],
            depth=self.template['thin_filament']['backbone']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['thin_filament']['backbone']['vertices'],
            rotation=(0, 0, 0))

        a_backbone = bpy.context.object
        mesh = a_backbone.data
        mesh.materials.append(self.materials['a_backbone'])

        self.hs_b['primitives']['a_backbone'] = a_backbone

    def create_primitive_titin(self):
        """ Create a titin primitive """

        bpy.ops.mesh.primitive_cylinder_add(
            radius=self.template['titin']['radius'],
            depth=self.template['titin']['depth'],
            enter_editmode=False,
            location=(0, 0, 0),
            vertices=self.template['titin']['vertices'],
            rotation=(0, 0, 0))

        titin = bpy.context.object
        mesh = titin.data
        mesh.materials.append(self.materials['titin'])

        self.hs_b['primitives']['titin'] = titin

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

    def create_primitive_end_disks(self):
        """ Creates M and Z disks """

        for x in ['m', 'z']:
            bpy.ops.mesh.primitive_cube_add(
                size=1,
                enter_editmode=False,
                location=(0,0,0))
            
            b = bpy.context.object
            mesh = b.data
            m_name = '%s_disk' % x
            mesh.materials.append(self.materials[m_name])

            self.hs_b['primitives'][m_name] = b

    def create_bs_primitives(self):
        """ Cycles through bs states making a mesh for each state """

        # Loop through states
        for i in range(0,2):
            s = self.template['thin_filament']['states'][i]

            # Make the mesh
            bpy.ops.mesh.primitive_cylinder_add(
                radius=self.template['thin_filament']['bs']['radius'],
                depth=self.template['thin_filament']['bs']['depth'],
                enter_editmode=False,
                location=(0, 0, 0),
                vertices=self.template['thin_filament']['bs']['vertices'],
                rotation=(0, 0, 0))

            # Add the material to get the color
            m = bpy.context.object
            mesh = m.data
            m_name = 'a_state_%i' % (i+1)
            mat = bpy.data.materials.new(name=m_name)
            mat.diffuse_color = s['color']
            mesh.materials.append(mat)

            self.hs_b['primitives'][m_name] = m
            
        # Create bs_sphere
        bpy.ops.mesh.primitive_ico_sphere_add(
            radius = self.template['thin_filament']['sphere']['radius'],
            enter_editmode=False,
            subdivisions=4,
            location=(0,0,0))
        
        m = bpy.context.object
        mesh = m.data
        m_name = 'bs_sphere'
        mat = bpy.data.materials.new(name=m_name)
        mat.diffuse_color = self.template['thin_filament']['sphere']['color']
        mesh.materials.append(mat)
        
        self.hs_b['primitives']['bs_sphere'] = m

    def return_all_filament_y_coords(self):
        """Returns a list of all filament y coordinates"""

        m_ys = [thick_obj.m_y for thick_obj in self.hs.thick_fil]
        a_ys = [thin_obj.a_y for thin_obj in self.hs.thin_fil]

        return m_ys + a_ys

    def return_all_filament_z_coords(self):
        """Returns a list of all filament z coordinates"""

        m_zs = [thick_obj.m_z for thick_obj in self.hs.thick_fil]
        a_zs = [thin_obj.a_z for thin_obj in self.hs.thin_fil]

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
