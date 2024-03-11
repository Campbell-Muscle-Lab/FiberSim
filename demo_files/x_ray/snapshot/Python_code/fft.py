# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 11:54:42 2024

@author: Campbell
"""

import numpy as np
import matplotlib.pyplot as plt

from PIL import Image, ImageChops

def fft():
    
    fs = "d:/ken/github/campbellmusclelab/models/fibersim/demo_files/x_ray/snapshot/renders/hs_1.png"
    
    im = Image.open(fs).convert('L')
    im = ImageChops.invert(im)
    im_data = np.asarray(im.getdata()).reshape(im.size[1], im.size[0])
    
    print(im.size)
    
    print(im_data)
    
    s = im.size[1]
    s = s/2
    
    # im_data = im_data[0:(2*int(s)), :]
    # print(im_data)
    
    # print(im_data.size)
    
    
    x = np.arange(-s, s)
    y = np.arange(-s, s)
    X, Y = np.meshgrid(x,y)
    
    print(im_data[2][3])
    
    grating = 0 * np.sin(X)
    for (c, xi) in enumerate(x):
        for (r, yi) in enumerate(y):
            grating[r, c] = im_data[r][c]
    
    
    

    plt.set_cmap("gray")
    plt.subplot(131)
    plt.imshow(im_data)
    plt.subplot(132)
    plt.imshow(grating)
    # Calculate Fourier transform of grating
    ft = np.fft.ifftshift(grating)
    ft = np.fft.fft2(ft)
    ft = np.fft.fftshift(ft)
    plt.subplot(133)
    plt.imshow(np.log10(abs(ft)))
    # plt.xlim([480, 520])
    # plt.ylim([520, 480])  # Note, order is reversed for y
    plt.show()
    
if __name__ == "__main__":
    fft()