# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:03:12 2023

@author: utkug
"""

from tkinter import *
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog

import tkinter as tk
import pathlib 
import os
import re
import importlib
import sys
import subprocess
from PIL import Image, ImageTk
import json
import collections


class ValueTypes:
    DICT = 1
    LIST = 2
    STR = 3
    FILEPATH = 4
    
class Tags:
    DICT = 'dict'
    LIST = 'list'
    ROOT = 'root'
    LEAF = 'leaf'
    FILE = 'file'


class Main(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        
        self.editor_menu_set = collections.OrderedDict()
        self.editor_menu_set['edit_child'] = {'text': 'Edit',
                                          'action': lambda: self.EditJSON()}
        self.FiberPyPath = StringVar()
        self.SetAppSize()
        self.DivideRowsColumns()
        self.SimSessionPanel()
        self.SimInputPanel()
        self.RightPanel()

    def SetAppSize(self):
        
        self.title("FiberSim")
        self.iconbitmap(default="favicon.ico")
        self.output_label = []

        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        print("Screen width:", screen_width)
        print("Screen height:", screen_height)
        
        app_window_width = screen_width * 0.75
        app_window_height = screen_height * 0.75
        
        x = (screen_width/2) - (app_window_width/2)
        y = (screen_height/2) - (app_window_height/2)
        
        self.geometry('%dx%d+%d+%d' % (app_window_width,app_window_height, 
                                       x, y))
        self.attributes('-topmost',1)
        
    def DivideRowsColumns(self):
            
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1,weight=10)
        #self.columnconfigure(2,weight=1)
        
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1,weight=2)
        self.rowconfigure(2,weight=80)
            
        
    def SimSessionPanel(self):
        
        left_frame_top = tk.LabelFrame(self, text="Simulation Session Panel")
        left_frame_top.grid(row=0,column=0,rowspan=1,
                            sticky='WENS',padx=10,pady=10)

        locate_fiberpy_button = ttk.Button(left_frame_top,
                                           text="Select FiberPy Folder", 
                                           command=self.LocateFiberPy)
        locate_fiberpy_button.grid(row=0, column=0, padx=10,pady=10, sticky='W')
        
        
        wid = int(self.winfo_width())
        folder_path_text = ttk.Entry(left_frame_top, width = 50*wid, 
                                     textvariable = self.FiberPyPath)
        folder_path_text.grid(row=0,column=1)
        
        rad_but = {"Demo Simulations":"demo",
                   "Custom Simulations": "custom"}
        self.radio = tk.StringVar()
        
        session_type = ["demo","custom"]
        self.session = {}
        m = 0
        
        for (text, value) in rad_but.items():
            
            ix = m + 2
            st = session_type[m]
            self.session[st] = ttk.Radiobutton(left_frame_top, text = text, 
                            variable = self.radio, 
                            value=value,
                            command=self.RadioButtonSelected).grid(row=0,
                                                                   column=ix, 
                                                                   padx = 10)
            
            self.session_type = st
            print(self.session_type)
            m += 1
            
        left_frame_top.update()    
        self.sim_sess_wid = left_frame_top.winfo_width()
        self.sim_sess_heg = left_frame_top.winfo_height()
        
    def RadioButtonSelected(self):
        
        self.SimInputPanel()
        
    def SimInputPanel(self):
        
        self.sim_input_panel = tk.LabelFrame(self, 
                                             text="Simulation Input Panel",
                                             width=self.sim_sess_wid,
                                             height=self.sim_sess_heg)
        self.sim_input_panel.grid(row=1,column=0,rowspan=1,
                            sticky='WENS',padx=10,pady=10)
        
        if self.radio.get() == 'demo':
            demo_list = ["Isometric Activation", "Ramp Shortening",
                      "Isotonic Shortening","Isometric Twitch"]
            self.demo_selection = ttk.Combobox(self.sim_input_panel, values=demo_list)
            self.demo_selection.set("Select a demo")
            self.demo_selection.grid(row=1, column=0, padx=10, pady=10)
            self.demo_selection.bind("<<ComboboxSelected>>",self.GenerateDemoPath)
        
            run_demo_but = ttk.Button(self.sim_input_panel, 
                                  text="Run Demo", command = self.RunDemo)
            run_demo_but.grid(row=1, column=1, padx=10,pady=10, sticky='W')
        
        elif self.radio.get() == 'custom':
            
            self.custom_files = {}
            self.but_indicators = {}
            self.custom_buttons = {}
            self.canvas = {}
            
            but_text = ["Batch File","Model File",
                        "Options File","Output File",
                        "Summary File"]
            
            col = ['red','red','red','red','red']
            
            
            
            
            m = 0
            k = 1
            for i in range(len(but_text)):
                
                but_key = but_text[i]
                print(but_key)
                but_key = but_key.lower()
                but_key = re.sub(r"\s","_",but_key)
                
                self.custom_files[but_key[i]]={}
                self.custom_buttons[but_key[i]] = ttk.Button(self.sim_input_panel, 
                                           text=but_text[i],
                                           command=lambda jbut = but_key
                                           :self.UploadJSON(jbut))
                

                
                self.custom_buttons[but_key[i]].grid(row=1, column=i,
                                                     padx=10,pady=10, 
                                                     sticky='W')
                
                m = m + 2
                
                self.canvas[but_key[i]] = Canvas(self.sim_input_panel, 
                                                         width=50, height=50)

                self.canvas[but_key[i]].grid(row=2, column=i, padx=10)
                
                k = k + 2
                self.but_indicators[but_key[i]] = self.canvas[but_key[i]].create_oval(10, 10, 40, 40, 
                                                            width=0, 
                                                            fill=col[i])
                
                print(but_key)
    def LeftPanelBottom(self):
        
        left_frame_bottom = tk.LabelFrame(self,text="JSON File Editor")
        left_frame_bottom.grid(row=2, column=0,rowspan=2,
                               sticky='WENS',padx=10,pady=10)
        
        self.tabs = ttk.Notebook(left_frame_bottom)
        
        tab = {}
        self.tree = {}
        self.dat = {}
        k = 1
        
        self.tab_name=["Batch","Model","Options","Output Handler","Template Summary"]
        for i in range(len(self.tab_name)):
            key = 'tab'+ str(i)
            # print(key)
            tab[key] = ttk.Frame(self.tabs)
            self.tabs.add(tab[key],text=self.tab_name[i])
            # ttk.Label(tab[key]).grid(column=0, row=1)
            self.dat[key] = {}
            self.tree[key] = ttk.Treeview(tab[key], selectmode='browse',height=30)
            self.tree[key].pack(fill=tk.BOTH,padx=10,pady=10)
            self.tree[key].bind('<Button-3>', 
                                 lambda event: self.ShowEditorMenu(event))
            self.editor_menu = tk.Menu(self.tree[key], tearoff=0)
            self.LoadEditorMenu()
            self.key = key
            
            self.SetEditorColumns()

        self.tabs.pack(expand=1, fill="both")
        
        self.tabs.bind('<<NotebookTabChanged>>',self.NavigateJSON)       
    
    def NavigateJSON(self,event):
        
        jbut = self.tabs.tab(self.tabs.select(),"text")
        jix = self.tabs.index(self.tabs.select())
        self.key = 'tab'+str(jix)
        # print(self.key)
        jbut = jbut.lower()
        jbut = re.sub(r"\s","_",jbut)
        
        self.LocateJSON(jbut)        
        
        
    # def JSONEditor(self):
        
    #     json_editor = tk.LabelFrame(self,text="JSON Editor")
    #     json_editor.grid(row=2, column=0,rowspan=1,
    #                            sticky='WENS',padx=10,pady=10)

    #     self.tree = ttk.Treeview(json_editor, selectmode='browse',height=30)
    #     self.tree.pack(fill=tk.BOTH,padx=10,pady=10)
        
    #     self.SetEditorColumns()

    def SetEditorColumns(self, columns=('Key/Parameter', 'Value')):
        
        col_ids = ['#'+str(i) for i in range(len(columns)-1)]
        self.tree[self.key].configure(column=col_ids)
        for i in range(len(columns)):
            self.tree[self.key].heading('#'+str(i), text=columns[i])
            
    def UploadJSON(self,jbut):
        print(jbut)
        self.custom_files[jbut]=filedialog.askopenfilename(filetypes=[("JSON File"
                                                                       ,".json")])
        
        self.canvas[jbut].itemconfig(self.but_indicators[jbut], fill='green')
        
                    
    def LocateJSON(self,jbut):
        self.jbut = jbut
        print(jbut)
        if self.jbut == "batch":
            batch_json_path = self.demo_file
            self.jbut_file = batch_json_path
            # print(batch_json_path)
            self.LoadJSON()
            
        else:
            json_file_name = self.jbut + ".json"
            # print(json_file_name)

            def find(name, path):
                for root, dirs, files in os.walk(path):
                    if name in files:
                        return os.path.join(root, name)
                    
            self.jbut_file = find(json_file_name,self.demo_folder)
            # print(self.jbut_file)
            self.LoadJSON()

            
    def LoadJSON(self):
        
        jf = open(self.jbut_file,'r')
        if not self.dat[self.key]:
            self.dat[self.key] = json.load(jf)
            self.AddItemJSON(jf.name,self.dat[self.key],tags=[Tags.FILE])
        else:
            return
    
    def SaveJSON(self, filepath, data):
        
        # print(data)
        with open(filepath,'w') as jf:
        # print(self.dat[self.key])
            json.dump(data ,jf,indent=4)
    
    def AddItemJSON(self,key,value,parent='',tags=[]):
        
        self.PopulateItem(key,value,parent,tags)
        if parent == '':
            return
        json_root = self.GetJSONRoot(parent)
        if self.tree[self.key].tag_has(Tags.FILE, json_root):
            self.SaveJSON(self.GetJSONFilePath(json_root), 
                                self.GetValue(json_root))
        
    def PopulateItem(self, key, value, node='', tags=[]):
        
        if node == '':
            tags = tags+[Tags.ROOT]
            
        if type(value) is dict:
            # print(type(value))
            node = self.tree[self.key].insert(node, tk.END, text=str(key)+'={}',
                                    tags=tags+[Tags.DICT])
            print(node)
            for k in value:
                # print(value)
                self.PopulateItem(k,value[k],node)
        elif type(value) is list:
            # print(type(value))
            node = self.tree[self.key].insert(node, tk.END, text=str(key)+'=[]',
                                    tags=tags+[Tags.LIST])
            print(node)
            for k in range(len(value)):
                self.PopulateItem(k,value[k],node)
        else:
            print(value)
            # value = float(value)
            self.tree[self.key].insert(node, tk.END, text=str(key), 
                             tags=tags+[Tags.LEAF], values=[value])
    
        # print(self.tree[self.key])
        childrens = self.tree[self.key].get_children()
        for child in childrens:
            self.ExpandItems(child)
        
    def ExpandItems(self,index):

        self.tree[self.key].item(index, open=True)
        childrens = self.tree[self.key].get_children([index])

        if len(childrens) > 0:
            for child in childrens:
                self.ExpandItems(child)
                
    def EditJSON(self):
        
        index = self.tree[self.key].selection()
        is_leaf = self.tree[self.key].tag_has(Tags.LEAF,index)
        is_parent_dict = self.tree[self.key].tag_has(Tags.DICT, self.tree[self.key].parent(index))

        
        if is_leaf and tk.messagebox.askyesno("Edit Value",
                                   "Are you sure to edit value?"):
            value = simpledialog.askstring("Value Input","Enter new value: ")
            print(value)
            if self.VerifyValue(value):
                self.EditItemJSON(index,value=value)
                
    def EditItemJSON(self,index,key=None,value=None):
        
        # print(value)
        if key:
            # print(value)
            self.tree[self.key].item(index, text=key)

        if value:
            # print(value)
            self.tree[self.key].item(index, values=[value])

        if index == '':
            return
        json_root = self.GetJSONRoot(index)
        # print(json_root)
        if self.tree[self.key].tag_has(Tags.FILE, json_root):
            # print(self.GetValue(json_root))
            self.SaveJSON(self.GetJSONFilePath(json_root), self.GetValue(json_root))
            
        
    def VerifyValue(self,value):
        if len(value.encode('utf-8')):
            return True
        return False
    
    def GetJSONRoot(self,index):
        
        if index == '': return None
        if self.tree[self.key].tag_has(Tags.ROOT, index):
            return index
        return self.GetJSONRoot(self.tree[self.key].parent(index))
    
    
    def GetJSONFilePath(self,index):
        return self.GetKey(self.GetJSONRoot(index))
    
    def GetKey(self, index):
        item = self.tree[self.key].item(index)
        key = item['text']
        if self.tree[self.key].tag_has(Tags.DICT, index) or self.tree[self.key].tag_has(Tags.LIST, index):
            key = item['text'].split('=')[0].strip()
        return key
    
    def GetValue(self, index):

        key_text = {}        
        item = self.tree[self.key].item(index)
        value = None
        if self.tree[self.key].tag_has(Tags.DICT, index):
            value = {}
            child_nodes = self.tree[self.key].get_children(index)
            for child in child_nodes:
                # print("utku")
                key_text = self.GetKey(child)
                value[self.GetKey(child)] = self.GetValue(child)
                print(value)
        elif self.tree[self.key].tag_has(Tags.LIST, index):
            value = []
            child_nodes = self.tree[self.key].get_children(index)
            for child in child_nodes:
                key_text = self.GetKey(child)
                value.append(self.GetValue(child))
        else:
            value = item['values'][0]
        
        print(type(value))
        print(key_text)



        if isinstance(value, int):
            value = value
        else:
            try:
                value = float(value)
            except:
                value = value
        return value

        
    def RightPanel(self):
        
        self.right_frame = tk.LabelFrame(self, text="Simulation Output Panel")
        self.right_frame.grid(row=0,column=1,rowspan=3,columnspan=1,sticky='WENS',padx=10,pady=10)
                
    def LocateFiberPy(self):
        
        folder_name = filedialog.askdirectory()
        self.FiberPyPath.set(folder_name)
        self.GuiPath = os.getcwd()
        # print(self.GuiPath)
        os.chdir(folder_name)
        utku = os.getcwd()
        # print(folder_name)
        
    def GenerateDemoPath(self,event):
        
        gui_path = self.GuiPath
        main_folder = os.path.split(gui_path)[0]
        # print(main_folder)
        demo_name = self.demo_selection.get()
        # print(demo_name)
        demo_name = re.sub(r"\s","_",demo_name)
        demo_name = demo_name.lower()
        # print(demo_name)
        self.demo_folder = "..\\..\\..\\" + "demo_files\\getting_started\\" + demo_name
        # print(self.demo_folder)
        files = []
        for file in os.listdir(self.demo_folder):
            if file.endswith('.json'):
                files.append(file)
        self.demo_file = []
        self.demo_file = files[0]
        self.demo_file = self.demo_folder + "\\" + self.demo_file
        # print(self.demo_file)
        self.LeftPanelBottom()
        


        
    def RunDemo(self):
        sys.argv = ["run_batch",self.demo_file]
        # print(os.getcwd())
        # print(sys.argv[1])
        # print('demo starts')
        subprocess.call(["python", "FiberPy.py", *sys.argv])
        # print('demo ends')
        self.OutputDisplay()
        
    def OutputDisplay(self):
        
        output_summary  = self.demo_folder + "\\sim_output\\" + "summary.png"
        print(output_summary)
        im = Image.open(output_summary)
        aspect_ratio = im.height/im.width
        print(aspect_ratio)
        self.right_frame.update()
        height = 0.9*self.right_frame.winfo_height()
        width = height / aspect_ratio
        width = int(width)
        height = int(height)
        im_small = im.resize((width,height))
        im_r = ImageTk.PhotoImage(im_small)
        if not self.output_label:
            self.output_label = Label(self.right_frame, image = im_r,anchor= CENTER)
            self.output_label.image = im_r
        else:
            self.output_label.configure(image = im_r)
            # self.output_label = Label(self.right_frame, image = im_r,anchor= CENTER)
            self.output_label.image = im_r
        #label.grid(row=0,column=1,rowspan=3,columnspan=1,sticky='EW')
        self.output_label.pack()
   
    def ShowEditorMenu(self, event):
        
        # print("hello")
        self.tree[self.key].selection_set(self.tree[self.key].identify_row(event.y))
        if self.tree[self.key].selection():
            item = self.tree[self.key].selection()
            curItem = self.tree[self.key].focus()
            print(self.tree[self.key].item(curItem))
            self.editor_menu.post(event.x_root, event.y_root)
            
    def LoadEditorMenu(self):
        self.editor_menu.delete(0,tk.END)
        for i in self.editor_menu_set:
            self.editor_menu.add_command(label=self.editor_menu_set[i]['text'],
                                         command=self.editor_menu_set[i]['action'])
        
        
        
        
        
        
        


main = Main()


# main.after(25000,lambda:main.destroy())
main.mainloop()    


