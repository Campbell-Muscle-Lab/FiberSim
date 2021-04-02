import os
import time
import json
import sys


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

def run_render_batch():
    """ Entry point for the generate_model
        second_last argument is render_file
        last argument is job index """

    # Unpack the arguments
    render_file = sys.argv[-2]
    job_index = int(sys.argv[-1])

    # Open the render file and pull off the correct job
    with open(render_file, 'r') as f:
        render_batch = json.load(f)
    
    render_job = render_batch['render_batch']['render_jobs'][job_index]

    # Pull out the files
    if (render_job['relative_to'] == 'this_file'):
        parent_path = Path(render_file).parent
        template_file = os.path.join(parent_path,
                                     render_job['template_file'])
        frames_file = os.path.join(parent_path,
                                   render_job['frames_file'])
        blender_file = os.path.join(parent_path,
                                    render_job['blender_file'])
    
    # Load json files
    with open(template_file, 'r') as f:
        template = json.load(f)
    template_data = template['render_template']
    
    with open(frames_file, 'r') as f:
        frames = json.load(f)
    frames_data = frames['frames']
    
    with open(blender_file, 'r') as f:
        blender = json.load(f)
    blender_data = blender['blender_data']
    
    print(frames_data)
    print(template_data)
    print(blender_data)
    
    for fr in frames_data:
        print(fr)
        fr['parent_file'] = render_file
        run_render_job(fr, template_data, blender_data)
    

def run_render_job(frame_data, template_data, blender_data):
    
    parent_path = Path(frame_data['parent_file']).parent
    status_file = os.path.join(parent_path, frame_data['status_file'])
        
    h = hs.half_sarcomere(status_file)
    print(h)
    

    return

if __name__ == "__main__":
    run_render_batch()