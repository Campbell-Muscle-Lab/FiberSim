"""
Interface for blender rendering of images
@author: Ken Campbell
"""

import os
import json
import subprocess

from pathlib import Path

def generate_images(render_file_string):
    """ function generates images defined by image_dict """

    # Load the render file as a dict
    with open(render_file_string, 'r') as f:
        data = json.load(f)
    
    print(data)
    
    # Generate the command string
    command_string = ('cd "%s"\n ' % data['blender_exe_path'])

    # Deduce the path to this folder
    parent_path = Path(__file__).parent
    
    # Set path to generate.py
    generate_path = os.path.join(parent_path, 'generate.py')

    # Complete the commmand line
    command_string = command_string + \
         ('blender --python "%s" -- -j %s' %
              (generate_path, render_file_string))
    
    # Write command to temp.bat
    with open('run.bat', 'w') as f:
        f.write('%s' % command_string)
    
    print('run.bat')
    subprocess.call('run.bat')

    return

