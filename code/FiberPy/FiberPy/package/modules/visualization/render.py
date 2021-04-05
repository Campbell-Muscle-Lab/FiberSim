"""
Interface for blender rendering of images
@author: Ken Campbell
"""

import os
import json
import subprocess

from pathlib import Path

def generate_images(render_file):
    """ Function generates images defined by image_dict
        Opens a render batch file, finds the information about blender,
        then launches blender, which calls a python file """

    # Load the render batch file as a dict
    with open(render_file, 'r') as f:
        render_batch = json.load(f)
    
    render_jobs = render_batch['render_batch']['render_jobs']

    # Loop through batches
    for job_number, rj in enumerate(render_jobs):

        if (rj['relative_to'] == 'this_file'):
            parent_path = Path(render_file).parent
            blender_file = os.path.join(parent_path, rj['blender_file'])

        with open(blender_file, 'r') as f:
            blender_data = json.load(f)
    
        # Find the blender executable
        blender_exe_path = blender_data['blender_data']['blender_exe_path']
    
        # Generate the command string
        command_string = ('cd "%s"\n ' % blender_exe_path)
    
        # Deduce the path to this folder
        parent_path = Path(__file__).parent
        
        # Set path to generate.py
        generate_path = os.path.join(parent_path, 'generate2.py')
    
        # Complete the commmand line
        command_string = command_string + \
              ('blender -noaudio --python "%s" -- -j %s %i' %
                  (generate_path,
                   os.path.abspath(render_file),
                   job_number))
        
        # Write command to temp.bat
        with open('run.bat', 'w') as f:
            f.write('%s' % command_string)
        
        print('run.bat')
        subprocess.call('run.bat')

    return

