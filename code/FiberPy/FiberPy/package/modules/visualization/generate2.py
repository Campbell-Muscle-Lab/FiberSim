import os
import time
import json
import sys

import threading
import multiprocessing

# Blender packages that are located in Blender's Python Path.
import bpy
import bmesh
import mathutils

from pathlib import Path

# Generate the path that allows us to add half-sarcomere module
parent_dir = os.path.realpath(os.path.dirname(__file__))
sys.path.append(parent_dir)
print(sys.path)

import half_sarcomere as hs
import hs_blender as hs_b

def run_render_batch():
    """ Entry point for the generate_model
        last argument is render_job_file """

    # Unpack the arguments
    render_job_file = sys.argv[-1]

    # Open the render file and pull off the job
    with open(render_job_file, 'r') as f:
        render_job = json.load(f)
    
    render_job = render_job['render_job']

    # Pull out the files
    if (render_job['relative_to'] == 'this_file'):
        parent_path = Path(render_job_file).parent
        template_file = os.path.join(parent_path,
                                      render_job['template_file'])
        frame_file = os.path.join(parent_path,
                                    render_job['frame_file'])
        blender_file = os.path.join(parent_path,
                                    render_job['blender_file'])

    # Load json files
    with open(template_file, 'r') as f:
        template = json.load(f)
    template_data = template['render_template']

    with open(frame_file, 'r') as f:
        frame_data = json.load(f)
    frame_data = frame_data['frame_data']
    frame_data['parent_file'] = os.path.abspath(frame_file)

    with open(blender_file, 'r') as f:
        blender = json.load(f)
    blender_data = blender['blender_data']
    
    print(frame_data)
    print(template_data)
    print(blender_data)
    
    # Run the job
    run_render_job(frame_data, template_data, blender_data)

def run_render_job(frame_data, template_data, blender_data):
    
    parent_path = Path(frame_data['parent_file']).parent
    status_file = os.path.join(parent_path, frame_data['status_file'])
    output_image_file = os.path.join(parent_path, frame_data['image_file'])
    
        
    h = hs.half_sarcomere(status_file)

    hs_blend = hs_b.hs_blender(h, frame_data, template_data, blender_data,
                               output_image_file)

    #bpy.ops.wm.quit_blender()

    return

if __name__ == "__main__":
    run_render_batch()