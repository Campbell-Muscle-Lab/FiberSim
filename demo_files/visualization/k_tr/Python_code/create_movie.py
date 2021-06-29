# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 16:30:13 2021

@author: kscamp3
"""

import os
import cv2
import natsort

def create_movie(im_folder, movie_file):
    """ Create movie from folder """
    
    # Get files in folder
    ken=1
    im_files = [os.path.join(im_folder, f) for f in os.listdir(im_folder)
                 if os.path.isfile(os.path.join(im_folder,f))]
    im_files = natsort.natsorted(im_files)
    print(im_files)
    
    im_array=[]
    for f in im_files:
        print(f)
        im = cv2.imread(f)
        height, width, layers = im.shape
        size = (width, height)
        im_array.append(im)
    
    out = cv2.VideoWriter(movie_file,
                          cv2.VideoWriter_fourcc(*'MP4V'),
                          15,
                          size)
    
    for i in range(len(im_array)):
        print(i)
        out.write(im_array[i])
    out.release()
    print('Move written to: %s' % movie_file)
    
    
    
if __name__ == "__main__":
    im_folder = '../renders'
    movie_file = '../FiberSim_movie.mp4'
    
    create_movie(im_folder, movie_file)