# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:03:12 2023

@author: utkug
"""

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import tkinter as tk
import pathlib 
import os
import sys
import re


class Main(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.FiberPyPath = StringVar()
        self.SetAppSize()
        self.DivideRowsColumns()
        self.LeftPanel()
        self.RightPanel()
        

    def SetAppSize(self):
        
        self.title("FiberSim")
        self.iconbitmap("favicon.ico")
        
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        print("Screen width:", screen_width)
        print("Screen height:", screen_height)
        
        app_window_width = screen_width * 0.5
        app_window_height = screen_height * 0.5
        
        x = (screen_width/2) - (app_window_width/2)
        y = (screen_height/2) - (app_window_height/2)
        
        self.geometry('%dx%d+%d+%d' % (app_window_width,app_window_height, x, y))
        self.attributes('-topmost',1)
        
    def DivideRowsColumns(self):
            
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1,weight=1)
        self.columnconfigure(2,weight=1)
        
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1,weight=1)
        self.rowconfigure(2,weight=1)
            
        
    def LeftPanel(self):
        
        left_frame = tk.LabelFrame(self, text="Simulation Input Panel")
        left_frame.grid(row=0,column=0,rowspan=3,sticky='WENS')
        
        locate_fiberpy_button = ttk.Button(left_frame,text="Select FiberPy Folder", command=self.LocateFiberPy)
        locate_fiberpy_button.grid(row=0, column=0, padx=10,pady=10, sticky='W')
        
        wid = int(self.winfo_width()/4)
        folder_path_text = ttk.Entry(left_frame, width = wid, textvariable = self.FiberPyPath)
        folder_path_text.grid(row=0,column=1)
        
        demo_list = ["Isometric Activation", "Ramp Shortening","Isotonic Shortening","Isometric Twitch"]
        demo_selection = ttk.Combobox(left_frame, values=demo_list)
        demo_selection.set("Select a demo")
        demo_selection.grid(row=1, column=0, padx=10, pady=10)
        demo_selection.bind("<<ComboboxSelected>>",self.GenerateDemoPath)
        
        run_demo_but = ttk.Button(left_frame, text="Run Demo")
        run_demo_but.grid(row=1, column=1, padx=10,pady=10, sticky=tk.W)
        
    # def MidPanel(self):
    
    def RightPanel(self):
        
        right_frame = tk.LabelFrame(self, text="Simulation Output Panel")
        right_frame.grid(row=0,column=1,rowspan=3,columnspan=2,sticky='WENS')
        
    def LocateFiberPy(self):
        
        folder_name = filedialog.askdirectory()
        self.FiberPyPath.set(folder_name)
        self.GuiPath = os.getcwd()
        print(self.GuiPath)
        os.chdir(folder_name)
        utku = os.getcwd()
        print(utku)

        
    def GenerateDemoPath(self,event):
        
        self.DemoPath = self.GuiPath
        print(self.DemoPath)
        
        
        
    

main = Main()
main.after(50000,lambda:main.destroy())
main.mainloop()    
        