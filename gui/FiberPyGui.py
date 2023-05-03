# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:38:55 2023

@author: utkug
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import tkinter as tk
import pathlib 
import os
import sys

print(pathlib.Path().absolute())


def set_app_size(app_window):
    
    app_window.title("FiberSim")
    app_window.iconbitmap("favicon.ico")

    screen_width = app_window.winfo_screenwidth()
    screen_height = app_window.winfo_screenheight()

    print("Screen width:", screen_width)
    print("Screen height:", screen_height)

    app_window_width = screen_width * 0.75
    app_window_height = screen_height * 0.75
    
    # app_window_width = 400
    # app_window_height = 600

    x = (screen_width/2) - (app_window_width/2)
    y = (screen_height/2) - (app_window_height/2)

    app_window.geometry('%dx%d+%d+%d' % (app_window_width,app_window_height, x, y))
    
    return app_window      
      

def locate_FiberPy():
    
    folder_name = filedialog.askdirectory()
    folderPath.set(folder_name)
    utku = os.getcwd()
    print(utku)
    os.chdir(folder_name)
    utku = os.getcwd()
    print(utku)


def gen_demo_path(folderPath):
    
    return
    
main_window = tk.Tk()
main_window = set_app_size(main_window)

main_window.columnconfigure(0,weight=1)
main_window.columnconfigure(1,weight=1)

main_window.rowconfigure(0,weight=1)
main_window.rowconfigure(1,weight=1)


frame_1 = tk.LabelFrame(main_window, text="Simulation Input Panel")
frame_2 = tk.Frame(main_window, bg='green')
frame_3 = tk.LabelFrame(main_window, text="Simulation Output Panel")

label_1 = tk.Label(frame_1, text="Simulation Input Panel")



frame_1.grid(row=0,column=0,sticky='WENS')
frame_2.grid(row=1,column=0,sticky='WENS')
frame_3.grid(row=0,column=1,rowspan=2,sticky='WENS')
#main_window = set_app_size(main_window)

folderPath = StringVar()
locate_fiberpy_but = ttk.Button(frame_1, text="Select FiberPy Folder", command=locate_FiberPy)
locate_fiberpy_but.grid(row=0, column=0, padx=10,pady=10, sticky=tk.W)

folder_path_text = ttk.Entry(frame_1,width=50,textvariable=folderPath)
folder_path_text.grid(row=0, column=1)

demo_list = ["Isometric Activation", "Ramp Shortening","Isotonic Shortening","Isometric Twitch"]
demo_selection = ttk.Combobox(frame_1, values=demo_list)
demo_selection.set("Select a demo")
demo_selection.grid(row=1, column=0, padx=10,pady=10)
demo_selection.bind("<<ComboboxSelected>>",gen_demo_path(folderPath))


run_demo_but = ttk.Button(frame_1, text="Run Demo")
run_demo_but.grid(row=1, column=1, padx=10,pady=10, sticky=tk.W)


 








main_window.mainloop()