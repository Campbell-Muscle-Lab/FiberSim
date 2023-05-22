# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:03:12 2023

@author: utkug
"""

from tkinter import *
from tkinter import ttk
from tkinter import font
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
import pandas as pd
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
        
        self.defaultFont = font.nametofont("TkDefaultFont")
        self.defaultFont.configure(family="Arial")

        self.editor_menu_set = collections.OrderedDict()
        self.editor_menu_set['edit_child'] = {'text': 'Edit',
                                          'action': lambda: self.EditJSON()}
        self.SetAppSize()
        self.DivideRowsColumns()
        self.SimSessionPanel()
        self.SimInputPanel()
        self.BatchJobCounter()
        self.SimOutputPanel()

    def SetAppSize(self):
        
        self.title("FiberSim")
        self.iconbitmap(default="favicon.ico")
        self.output_label = {}

        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        
        app_window_width = screen_width * 0.75
        app_window_height = screen_height * 0.75
        
        x = (screen_width/2) - (app_window_width/2)
        y = (screen_height/2) - (app_window_height/2)
        
        self.geometry('%dx%d+%d+%d' % (app_window_width,app_window_height, 
                                       x, y))
        self.attributes('-topmost',1)
        
        self.tab = {}
        self.tree = {}
        self.dat = {}
        self.FiberPyPath = StringVar()
        self.ProtocolPath = StringVar()
        self.ProtocolName = StringVar()
        self.files = {}
        self.f_struct = []
        self.disp = []
        self.job_tot = 0

        
    def DivideRowsColumns(self):
            
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1,weight=1)
        self.columnconfigure(2,weight=100)

        
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1,weight=3)
        self.rowconfigure(2,weight=60)
            
        
    def SimSessionPanel(self):
        
        left_frame_top = tk.LabelFrame(self, text="Simulation Session Panel")
        left_frame_top.grid(row=0,column=0,rowspan=1,columnspan=2,
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
            m += 1
            
        left_frame_top.update()    
        self.sim_sess_wid = left_frame_top.winfo_width()
        self.sim_sess_heg = left_frame_top.winfo_height()
        
    def RadioButtonSelected(self):

        if self.tree:
            for tabs in self.tree:
                self.dat[tabs] = {}
                for item in self.tree[tabs].get_children():
                    self.tree[tabs].delete(item)
        
        self.disp = []
        self.SimInputPanel()
        
    def SimInputPanel(self):
        
        self.sim_input_panel = tk.LabelFrame(self, 
                                             text="Simulation Input Panel",
                                             width=self.sim_sess_wid,
                                             height=self.sim_sess_heg)
        
        self.sim_input_panel.grid(row=1,column=0,rowspan=1,columnspan = 1,
                            sticky='WENS',padx=10,pady=10)
        
        if self.radio.get() == 'demo':
            self.demo_files = {}
            self.gs_demo_list = ["Isometric Activation", "Ramp Shortening",
                      "Isotonic Shortening","Isometric Twitch"]
            self.demo_list = ["Isometric Activation", "Ramp Shortening",
                      "Isotonic Shortening","Isometric Twitch","pCa Curves","K_tr"]
            self.demo_selection = ttk.Combobox(self.sim_input_panel, values=self.demo_list)
            self.demo_selection.set("Select a demo")
            self.demo_selection.grid(row=1, column=0, padx=10, pady=10)
            self.demo_selection.bind("<<ComboboxSelected>>",self.GenerateDemoPath)
        
            run_demo_but = ttk.Button(self.sim_input_panel, 
                                  text="Run Demo", command = self.RunDemo)
            run_demo_but.grid(row=1, column=1, padx=10,pady=10, sticky='W')
        
        elif self.radio.get() == 'custom':
            
            self.but_indicators = {}
            self.custom_buttons = {}
            self.canvas = {}
            
            but_text = ["Batch File","Model File", "Protocol File",
                        "Options File","Output Handler File",
                        "Template Summary File"]
            ix = [0,0,2,2,4,4]
            lix = [1,1,3,3,5,5]
            rix = [0,1,0,1,0,1]

                        
            for i in range(len(but_text)):
                
                but_key = but_text[i]
                but_key = but_key.lower()
                but_key = re.sub(r"\s","_",but_key)
                
                self.files[but_key]={}
                self.custom_buttons[but_key] = ttk.Button(self.sim_input_panel, 
                                           text=but_text[i],
                                           command=lambda jbut = but_key
                                           :self.UploadJSON(jbut))
                

                
                self.custom_buttons[but_key].grid(row=rix[i], column=ix[i],
                                                     padx=5,pady=10, 
                                                     sticky='W')
                                
                self.canvas[but_key] = Canvas(self.sim_input_panel, 
                                                         width=25, height=25)
                
                self.canvas[but_key].grid(row=rix[i], column=lix[i], padx=5)
                
                self.but_indicators[but_key] = self.canvas[but_key].create_oval(5, 20, 20, 5, 
                                                            width=2, 
                                                            fill='red',
                                                            outline='black')
            run_sim_but = ttk.Button(self.sim_input_panel, 
                                  text="Run Simulation", command = self.RunSimulation)
            
            clear_sim_files_but = ttk.Button(self.sim_input_panel, 
                                  text="Clear Files", command = self.ClearFiles)
            
            run_sim_but.grid(row=0, column=lix[-1]+1, padx=10,pady=10, sticky='W')
            clear_sim_files_but.grid(row=1, column=lix[-1]+1, padx=10,pady=10, sticky='W')
        
        self.sim_input_panel.update()    
        self.sim_sess_wid2 = self.sim_input_panel.winfo_width()
        self.sim_sess_heg2 = self.sim_input_panel.winfo_height()

    def BatchJobCounter(self):
        
        self.batch_counter_panel = tk.LabelFrame(self, 
                                             text="Batch Job Counter",
                                             width=self.sim_sess_wid)
                
        self.batch_counter_panel.grid(row=1,column=1,columnspan = 1,
                            sticky='WENS',padx=10,pady=10)
    
        spin_label = tk.Label(self.batch_counter_panel, text='Batch Job No')
        spin_label.grid(row = 0,column = 0, pady=2,padx=1)
        var = StringVar(self.batch_counter_panel)
        var.set("0")
        self.spin = tk.Spinbox(self.batch_counter_panel, from_=0, to=0,textvariable=var,width=2,wrap=True,command=self.ChangeJob)
        self.spin.grid(row = 1,column = 0, pady=2,padx=5,sticky='WENS')

        sum_figure = ttk.Button(self.batch_counter_panel, text='Summary Figure', command = self.SummaryFigure)

        sum_figure.grid(row=2,column=0,columnspan = 1,
                                        sticky='W',padx=10,pady=10)

    def SummaryFigure(self):
        self.OutputDisplay(self.job_tot)

    def ChangeJob(self):
        
        tab = self.tabs.tab(self.tabs.select(),'text')
        tab = tab.lower()
        tab = re.sub(" ","_",tab)
        job_no = self.spin.get()
        self.disp = []
        job_no = int(job_no)
        self.LocateJSON(job_no)
        self.OutputDisplay(job_no)
        self.tabs.select(self.tab[tab])


        


    def GetProtocolFolder(self):
        folder_name = filedialog.askdirectory()
        self.ProtocolPath.set(folder_name)
        

    def EditorPanel(self):
        
        self.editor_panel_frame = tk.LabelFrame(self,text="Editor Panel")
        self.editor_panel_frame.grid(row=2, column=0,rowspan=2,columnspan=2,
                               sticky='WENS',padx=10,pady=10)
        
        self.tabs = ttk.Notebook(self.editor_panel_frame)
        

        k = 1
        
        self.tab_name=["Batch File","Model File","Protocol File","Options File",
                       "Output Handler File","Template Summary File"]
        
        for i in range(len(self.tab_name)):
           
            key = self.tab_name[i]
            key = key.lower()
            key = re.sub(r"\s","_",key)
            
            self.tab[key] = ttk.Frame(self.tabs)
            self.tabs.add(self.tab[key],text=self.tab_name[i])
            self.dat[key] = {}
            self.tree[key] = ttk.Treeview(self.tab[key], selectmode='browse',height=30)
            self.tree[key].pack(fill=tk.BOTH,padx=10,pady=10)
            self.tree[key].bind('<Button-3>', 
                                 lambda event: self.ShowEditorMenu(event))
            self.editor_menu = tk.Menu(self.tree[key], tearoff=0)
            self.LoadEditorMenu()
            
            self.SetEditorColumns(key)
            
        if self.radio.get() == 'demo':
            self.LocateJSON(key)

        self.tabs.pack(expand=1, fill="both")
        
        self.tabs.bind('<<NotebookTabChanged>>',self.NavigateJSON)
        if self.radio.get() == 'custom':
            self.LocateJSON(self.key)         
    
    def NavigateJSON(self,event):

        jbut = self.tabs.tab(self.tabs.select(),"text")
        jix = self.tabs.index(self.tabs.select())
        jbut = jbut.lower()
        jbut = re.sub(r"\s","_",jbut)
        self.key = jbut

    def SetEditorColumns(self, key, columns=('Key/Parameter', 'Value')):
        
        if key == 'protocol_file':
            columns = ('ix','dt','pCa','dhsl','mode')
            self.tree[key]['show'] = 'headings'

        col_ids = ['#'+str(i) for i in range(len(columns)-1)]
        self.tree[key].configure(column=col_ids)
        for i in range(len(columns)):
            self.tree[key].heading('#'+str(i), text=columns[i])
            
    def UploadJSON(self,jbut):

        self.key = jbut

        if jbut == 'protocol_file':
            text = "Text File"
            ext = ".txt"
        else:
            text = "JSON File"
            ext = ".json"
        
        c_file = re.sub('_',' ',jbut).title()

        dialog_title = 'Open %s' % (c_file)
        self.files[jbut]=filedialog.askopenfilename(filetypes=[(text,ext)],title=dialog_title)
        if self.f_struct:
            self.f_struct = []
        if not self.tree:    
            self.EditorPanel()
        else:
            self.LocateJSON(self.key)
            
    def LocateJSON(self,key):
        def find(name, path):
            for root, dirs, files in os.walk(path):
                if name in files:
                    return os.path.join(root, name)
        self.key = key
        
        with open(self.files['batch_file'], 'r') as f:
            d = json.load(f)

        

        batch_structure = d['FiberSim_batch']
        job_data = batch_structure['job']


        for i,j in enumerate(job_data):
            for f in ['model_file','protocol_file',
                      'options_file','output_handler_file']:
                self.files[f] = j[f]
                if (not j['relative_to']):
                    self.files[f] = os.path.abspath(self.files[f])
                elif (j['relative_to'] == 'this_file'):
                    base_directory = pathlib.Path(self.files['batch_file']).parent.absolute()
                    self.files[f] = os.path.join(base_directory, self.files[f])
                else:
                    base_directory = j['relative_to']
                    self.files[f] = os.path.join(base_directory, self.files[f])

            base_directory = pathlib.Path(self.files['output_handler_file']).parent.absolute()        
            with open(self.files['output_handler_file'],'r') as jf:
                d = json.load(jf)
            temp = d['templated_images']
            self.files['template_summary_file'] = temp[0]['template_file_string']
            self.files['template_summary_file'] = os.path.join(base_directory,self.files['template_summary_file'])
            
            res = temp[0]['output_file_string']
            res = os.path.join(base_directory,res) + ".png"
            self.disp.append(res)
            self.f_struct.append(self.files.copy())

        job_no = int(self.spin.get())
        base_directory = pathlib.Path(self.files['batch_file']).parent.absolute()
        try:
            batch_figures = batch_structure['batch_figures']
            for i in batch_figures:
                try:
                    res = batch_figures[i][0]['output_image_file']
                except:
                    res = batch_figures[i][0]['output_image_file_string']
                res = os.path.join(base_directory,res) + ".png"
                self.disp.append(res)
        except:
            pass
        
        if len(job_data) > 1:
            self.spin.config(to=len(job_data)-1)
        self.LoadJSON(job_no)
        self.tabs.select(self.tab['batch_file'])
        self.job_tot = len(job_data)

    def LoadJSON(self,job_no):

        f_struct = self.f_struct[job_no]
        for files in f_struct:
            self.key = files
            if self.key == 'protocol_file':    
                if self.tree[self.key]:
                    self.dat[self.key] = pd.read_csv(f_struct[self.key],delim_whitespace=True)
                    self.dat[self.key].columns = ['dt','pCa','dhsl','mode']
                for i in range(len(self.dat[self.key]['dt'])):
                    self.tree[self.key].insert('', 'end', values=(self.dat[self.key]['dt'][i],
                                                              self.dat[self.key]['pCa'][i],
                                                              self.dat[self.key]['dhsl'][i],
                                                              self.dat[self.key]['mode'][i]))
                for column in ['#1','#2','#3','#4']:
                    self.tree[self.key].column(column,anchor = CENTER)
                self.tabs.select(self.tab[self.key])

            else:

                if f_struct[self.key]:     

                    for item in self.tree[self.key].get_children():
                        self.tree[self.key].delete(item)
                        self.dat[self.key] = {}
                    with open(f_struct[self.key],'r') as jf:
                        if not self.dat[self.key]:
                            self.dat[self.key] = json.load(jf)
                            self.AddItemJSON(jf.name,self.dat[self.key],tags=[Tags.FILE])
                        else:
                            return
                        
            if self.radio.get() == "custom":
                self.canvas[self.key].itemconfig(self.but_indicators[self.key], fill='green')
                
    def SaveJSON(self, filepath, data):
        
        with open(filepath,'w') as jf:
            json.dump(data ,jf,indent=4)
    
    def AddItemJSON(self,key,value,parent='',tags=[]):
        
        self.tabs.select(self.tab[self.key])
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
            node = self.tree[self.key].insert(node, tk.END, text=str(key)+'={}',
                                    tags=tags+[Tags.DICT])
            for k in value:
                self.PopulateItem(k,value[k],node)
        elif type(value) is list:
            node = self.tree[self.key].insert(node, tk.END, text=str(key)+'=[]',
                                    tags=tags+[Tags.LIST])
            for k in range(len(value)):
                self.PopulateItem(k,value[k],node)
        else:
            self.tree[self.key].insert(node, tk.END, text=str(key), 
                             tags=tags+[Tags.LEAF], values=[value])
    
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
            if self.VerifyValue(value):
                self.EditItemJSON(index,value=value)
                
    def EditItemJSON(self,index,key=None,value=None):
        
        if key:
            self.tree[self.key].item(index, text=key)

        if value:
            self.tree[self.key].item(index, values=[value])

        if index == '':
            return
        json_root = self.GetJSONRoot(index)
        if self.tree[self.key].tag_has(Tags.FILE, json_root):
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
                key_text = self.GetKey(child)
                value[self.GetKey(child)] = self.GetValue(child)
        elif self.tree[self.key].tag_has(Tags.LIST, index):
            value = []
            child_nodes = self.tree[self.key].get_children(index)
            for child in child_nodes:
                key_text = self.GetKey(child)
                value.append(self.GetValue(child))
        else:
            value = item['values'][0]
        
        if isinstance(value, int):
            value = value
        else:
            try:
                value = float(value)
            except:
                value = value
        return value

        
    def SimOutputPanel(self):
        
        self.sim_output_panel = tk.LabelFrame(self, text="Simulation Output Panel")
        self.sim_output_panel.grid(row=0,column=2,rowspan=3,columnspan=1,sticky='WENS',padx=10,pady=10)
        
                
    def LocateFiberPy(self):
        
        folder_name = filedialog.askdirectory()
        self.FiberPyPath.set(folder_name)
        self.GuiPath = os.getcwd()
        os.chdir(folder_name)
        
    def GenerateDemoPath(self,event):

        if self.f_struct:
            self.f_struct = []
        gui_path = self.GuiPath
        main_folder = os.path.split(gui_path)[0]
        demo_name = self.demo_selection.get()

        if demo_name in self.gs_demo_list:
            folder_key = "getting_started\\"
        else:
            folder_key = ""
            
        demo_name = re.sub(r"\s","_",demo_name)
        demo_name = demo_name.lower()

        self.demo_folder = "..\\..\\..\\" + "demo_files\\" + folder_key + demo_name
        print(self.demo_folder)
        files = []
        for file in os.listdir(self.demo_folder):
            if file.endswith('.json'):
                files.append(file)
        print
        self.files['batch_file'] = files[0]
        
        self.files['batch_file'] = self.demo_folder + "\\" + self.files['batch_file']
        if self.disp:
            self.disp = []
        
        self.EditorPanel()

    def ClearFiles(self):
        for files in self.files:
            self.files[files] = {}
            self.dat[files] = {}
            self.canvas[files].itemconfig(self.but_indicators[files], fill='red')
            for item in self.tree[files].get_children():
                self.tree[files].delete(item)
                self.tabs.select(self.tab['batch_file'])

    def RunSimulation(self):

        for i in range(len(self.tab_name)):
            c_file = re.sub(r"\s","_",self.tab_name[i])
            c_file = c_file.lower()

            if not self.files[c_file]:
                text = '%s is missing please upload and run again' % (self.tab_name[i])
                tk.messagebox.showwarning("FiberSim",text)
                return
        sys.argv = ["run_batch",self.files['batch_file']]
        subprocess.call(["python","FiberPy.py",*sys.argv])
        job_no = int(self.spin.get())
        self.OutputDisplay(job_no)
        
    def RunDemo(self):
        sys.argv = ["run_batch",self.files['batch_file']]
        subprocess.call(["python", "FiberPy.py", *sys.argv])
        job_no = int(self.spin.get())
        self.OutputDisplay(job_no)
        
    def OutputDisplay(self,job_no):
        output_summary  = self.disp[job_no]
        im = Image.open(output_summary)
        aspect_ratio = im.height/im.width
        self.sim_output_panel.update()
        
        if aspect_ratio > 1:
            height = 0.9*self.sim_output_panel.winfo_height()
            width = height / aspect_ratio
        else:
            width = 0.9*self.sim_output_panel.winfo_width()
            height = width * aspect_ratio
            
        width = int(width)
        height = int(height)


        im_small = im.resize((width,height))
        im_r = ImageTk.PhotoImage(im_small)
        if job_no < self.job_tot:
            im_text = "Batch Job %s" % (job_no)
        else:
            im_text = "Summary"
        if not self.output_label:
            self.output_label = Label(self.sim_output_panel, text=im_text,
                                        image = im_r,anchor= CENTER,compound='bottom')
            self.output_label.image = im_r
        else:
            self.output_label.configure(image = im_r,text=im_text)
            self.output_label.image = im_r
        self.output_label.pack()
   
    def ShowEditorMenu(self, event):
        
        self.tree[self.key].selection_set(self.tree[self.key].identify_row(event.y))
        if self.tree[self.key].selection():
            item = self.tree[self.key].selection()
            curItem = self.tree[self.key].focus()
            self.editor_menu.post(event.x_root, event.y_root)
            
    def LoadEditorMenu(self):
        self.editor_menu.delete(0,tk.END)
        for i in self.editor_menu_set:
            self.editor_menu.add_command(label=self.editor_menu_set[i]['text'],
                                         command=self.editor_menu_set[i]['action'])
        
main = Main()
main.mainloop()    