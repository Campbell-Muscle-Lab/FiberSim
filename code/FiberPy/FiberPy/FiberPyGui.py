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
import time
import webbrowser
from tkinter.scrolledtext import ScrolledText



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
        s = ttk.Style()
        s.configure('Treeview', rowheight=30)   

        self.editor_menu_set = collections.OrderedDict()
        self.editor_menu_set['edit_child'] = {'text': 'Edit',
                                          'action': lambda: self.EditJSON()}
        

        self.SetAppSize()
        self.DivideRowsColumns()
        self.SimSessionPanel()
        self.SimOutputPanel()

    def SetAppSize(self):
        
        self.title("FiberSim")
        self.iconbitmap(default="favicon.ico")
        self.output_label = {}

        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        
        app_window_width = screen_width * 0.7
        app_window_height = screen_height * 0.75
        
        x = (screen_width/2) - (app_window_width/2)
        y = (screen_height/2) - (app_window_height/2)
        
        self.geometry('%dx%d+%d+%d' % (app_window_width,app_window_height, 
                                       x, y))
        #self.attributes('-topmost',1)
        
        self.tab = {}
        self.tree = {}
        self.dat = {}
        self.FiberPyPath = StringVar()
        self.BatchPath = StringVar()
        self.ProtocolPath = StringVar()
        self.ProtocolName = StringVar()
        self.progress_completion = StringVar()
        self.files = {}
        self.f_struct = []
        self.disp = []
        self.job_tot = 0
        self.sim_res_but_visible = 0
        self.running_message = {}
        self.im_labels = []

        
    def DivideRowsColumns(self):
            
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1,weight=8)
        self.columnconfigure(2,weight=8)
        self.columnconfigure(3,weight=8)


        
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1,weight=1)
        self.rowconfigure(2,weight=60)
            
        
    def SimSessionPanel(self):
        
        self.session_panel = tk.LabelFrame(self, text="Simulation Session")
        self.session_panel.grid(row=0,column=0,
                            sticky='WENS',padx=10,pady=10,columnspan=2)

        locate_fiberpy_button = ttk.Button(self.session_panel,
                                           text="Select FiberPy Folder", 
                                           command=self.LocateFiberPy)
        locate_fiberpy_button.grid(row=0, column=0, padx=10,pady=10,sticky='WENS')
        
        load_batch_file_button = ttk.Button(self.session_panel,
                                           text="Load Batch File", 
                                           command=lambda jbut = 'batch_file'
                                           :self.UploadJSON(jbut))
        
        load_batch_file_button.grid(row=1, column=0, padx=10,pady=10, sticky='WENS')

        wid = int(self.winfo_width())
        folder_path_text = ttk.Entry(self.session_panel, 
                                     textvariable = self.FiberPyPath,width=100)
        folder_path_text.grid(row=0,column=1,sticky='WE')

        load_batch_file_text = ttk.Entry(self.session_panel, 
                                     textvariable = self.BatchPath,width=100)
        load_batch_file_text.grid(row=1,column=1,sticky='WE')

        self.run_sim_but = ttk.Button(self.session_panel, 
                                  text="Run FiberSim", command = self.RunSimulation)
        
        self.sim_res_but = ttk.Button(self.session_panel, 
                                  text="FiberSim Results", command = self.ResultsFolder)
            
        self.run_sim_but.grid(row=2, column=0, padx=10,pady=10, sticky='WENS')
        self.sim_sess_wid = self.session_panel.winfo_width()
        self.sim_sess_heg = self.session_panel.winfo_height()

        spin_label = tk.Label(self.session_panel, text='Batch Job No')
        # spin_label.grid(row = 1,column = 2, pady=5,padx=5)
        var = StringVar(self.session_panel)
        var.set("1")
        self.spin = tk.Spinbox(self.session_panel, from_=1, to=1,textvariable=var,width=2,wrap=True,command=self.ChangeJob)
        # self.spin.grid(row = 1,column = 3, pady=5,padx=1,sticky='W')
        self.session_panel.update()    
        self.sim_sess_wid2 = self.session_panel.winfo_width()
        self.sim_sess_heg2 = self.session_panel.winfo_height()

    def JobNavigator(self):

        self.job_list_panel = tk.LabelFrame(self, text="Batch Job List")
        self.job_list_panel.grid(row=1, column=0,
                            sticky='WENS',padx=5,pady=10,columnspan=1,rowspan=2)
        self.job_list = tk.Listbox(self.job_list_panel)
        self.job_list.pack(fill=BOTH,expand=True)
        self.job_list.bind("<Double-Button-1>", lambda event: self.ChangeJob(event))
        for i in range(self.job_tot):
            t = "Job %i" % (i+1)
            self.job_list.insert(i+1,t)
        print(self.job_tot)

    def ChangeJob(self,event):
        
        tab = self.tabs.tab(self.tabs.select(),'text')
        tab = tab.lower()
        tab = re.sub(" ","_",tab)
        # self.disp = []
        # job_no = int(job_no)-1
        self.ContainerReset()
        self.LocateJSON(tab)
        # self.OutputDisplay(job_no)
        print(tab)
        self.tabs.hide(self.tab['batch_file'])
        self.tabs.select(self.tab[tab])

    def EditorPanel(self):
        
        self.editor_panel_frame = tk.LabelFrame(self,text="Editor")
        self.editor_panel_frame.grid(row=1, column=1,rowspan=2,columnspan=1,
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

        self.tabs.pack(expand=1, fill="both")
        
        self.tabs.bind('<<NotebookTabChanged>>',self.NavigateJSON)
        self.LocateJSON(self.key)
        self.tabs.hide(self.tab['batch_file'])         

    
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

        try:
            self.results_listbox.delete('1.0', END) 
            self.results_listbox.images.clear()        

        except:
            pass

        try:
            self.job_list.delete(0,END)
        except:
            pass

        if self.job_tot != 0:
            self.job_tot = 0



        self.disp = []
        self.im_labels = []
        if self.sim_res_but_visible == 1:
            self.sim_res_but.grid_forget()

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
        self.BatchPath.set(self.files[jbut])

        if self.f_struct:
            self.f_struct = []
        if not self.tree:    
            self.EditorPanel()
            self.JobNavigator()
        else:
            self.ContainerReset()
            self.LocateJSON(self.key)
            print(self.job_tot)
            self.JobNavigator()
            self.tabs.hide(self.tab['batch_file'])

            
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

        try:
            print('trying list')
            selection = self.job_list.curselection()
            for item in selection:
                job_no = item
            if not job_no:
                job_no = 0
        except:
            print('trying list exception')
            job_no = 0
            
        base_directory = pathlib.Path(self.files['batch_file']).parent.absolute()
        try:
            batch_figures = batch_structure['batch_figures']
            self.single_res_folder = {}
            for i in batch_figures:
                try:
                    res = batch_figures[i][0]['output_image_file']
                except:
                    res = batch_figures[i][0]['output_image_file_string']
                res = os.path.join(base_directory,res) + ".png"
                res1 = res + ".png"
                self.im_labels.append(res1)
                self.disp.append(res)
                self.batch_res_folder = os.path.dirname(res)
        except:
            pass

        for i,j in enumerate(job_data):
            for f in ['model_file','protocol_file',
                      'options_file','output_handler_file','progress_file']:
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
            res1 = res + ".png"
            self.im_labels.append(res1)
            res = os.path.join(base_directory,res) + ".png"
            self.single_res_folder = os.path.dirname(res)
            self.disp.append(res)
            self.f_struct.append(self.files.copy())
        
        if len(job_data) > 1:
            self.spin.config(to=len(job_data))
        self.LoadJSON(job_no)
        self.tabs.select(self.tab['batch_file'])
        self.job_tot = len(job_data)
        print(self.job_tot)

    def ContainerReset(self):
        for i,j in enumerate(self.tab):
            for item in self.tree[j].get_children():
                self.tree[j].delete(item)
                self.dat[j] = {} 


    def LoadJSON(self,job_no):

        f_struct = self.f_struct[job_no]
        print(f_struct['template_summary_file'])
        for files in f_struct:
            self.key = files
            print(files)
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
            elif self.key == 'progress_file':
                pass
            else:
                if f_struct[self.key]:
                    with open(f_struct[self.key],'r') as jf:
                        if not self.dat[self.key]:
                            self.dat[self.key] = json.load(jf)
                            self.AddItemJSON(jf.name,self.dat[self.key],tags=[Tags.FILE])
                        else:
                            return
                
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
        
        self.sim_output_panel = tk.LabelFrame(self, text="Simulation Output")
        self.sim_output_panel.grid(row=0,column=2,rowspan=3,columnspan = 2, sticky='WENS',padx=5,pady=10)

        self.results_listbox = ScrolledText(self.sim_output_panel,width=10,wrap=WORD)
        self.results_listbox.pack(fill=tk.BOTH, side=tk.LEFT, expand=True)

        # self.results_listbox = tk.Listbox(self.sim_output_panel)
        # self.results_listbox.pack(side=LEFT, fill=BOTH,expand=True)
                
    def LocateFiberPy(self):
        
        folder_name = filedialog.askdirectory()
        self.FiberPyPath.set(folder_name)
        self.GuiPath = os.getcwd()
        os.chdir(folder_name)

    def RunSimulation(self):
        self.FiberSimInfo()
        if self.sim_res_but_visible == 1:
            self.sim_res_but.grid_forget()
        self.progress = ttk.Progressbar(self.session_panel,orient = HORIZONTAL,
            length = 100, mode = 'determinate')
        self.progress.grid(row = 2, column = 1,sticky='WE')

        self.progress_label = ttk.Label(self.session_panel, text="Simulation Status: 0%")
        self.progress_label.grid(column=2, row=2)
        for i in range(len(self.tab_name)):
            c_file = re.sub(r"\s","_",self.tab_name[i])
            c_file = c_file.lower()

            if not self.files[c_file]:
                text = '%s is missing please upload and run again' % (self.tab_name[i])
                tk.messagebox.showwarning("FiberSim",text)
                return

        self.PassToTerminal()

    def PassToTerminal(self):
        try:
            self.results_listbox.delete('1.0', END) 
            self.results_listbox.images.clear()        
        except:
            pass
        sys.argv = ["run_batch",self.files['batch_file']]
        self.p = subprocess.Popen(["python","FiberPy.py",*sys.argv])
        self.sim_res_but_visible = 1
        self.IsFiberSimRunning()
        
    def ResultsFolder(self):
        

        if self.single_res_folder:
            folder = self.single_res_folder
            os.startfile(folder)
        else:
            folder = self.batch_res_folder
            os.startfile(folder)
        
   
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
            
    def IsFiberSimRunning(self):
        try:
            token = self.p.poll()
            if token is None: 
                if not self.running_message:
                    print(self.f_struct['progress_file'])

                prog = 0
                for job in range(self.job_tot):
                    print(self.f_struct[job]['progress_file'])
                    f = open(self.f_struct[job]['progress_file'],'r')
                    p_j = f.read()
                    prog = prog + float(p_j)
                    print(prog)
                    

                prog = prog/self.job_tot
                self.progress['value'] = prog
                self.progress_label['text'] = 'Simulation Status: %.2f %% completed' % prog
                self.session_panel.update_idletasks()
                main.after(3000, self.IsFiberSimRunning)
            elif token == 0:
                self.p = {}
                self.sim_res_but.grid(row=2, column=1, padx=5, pady=5,sticky='W')
                self.running_message.destroy()
                self.running_message = {}
                self.progress_label.grid_forget()
                self.progress.grid_forget()
                self.OutputImages()
        except:
            pass

    def FiberSimInfo(self):
        def callback(url):
            webbrowser.open_new_tab(url)
            
        self.running_message = Toplevel()
        self.running_message.geometry('472x60')
        self.running_message.title('FiberSim')
        t1 = Label(self.running_message, text="FiberSim is running",font=('Arial', 11))
        t2 = Label(self.running_message, text="Visit our website in the meantime: https://www.campbellmusclelab.org",font=('Arial', 11),
                   fg="blue", cursor="hand2")
        t1.grid(row=0,column=0,padx=10)
        t2.grid(row=1,column=0,padx=10)
        t2.bind("<Button-1>", lambda e:callback("https://sites.google.com/g.uky.edu/cml/home"))

    def OutputImages(self):
        self.results_listbox.images = []
        self.results_listbox.delete('1.0', END)
        for i in range(len(self.disp)+1):
            img = Image.open(self.disp[i])

            aspect_ratio = img.height/img.width
            self.sim_output_panel.update()
        
            width = 0.6*self.results_listbox.winfo_width()
            height = width * aspect_ratio
            
            width = int(width)
            height = int(height)
            im_small = img.resize((width,height))
            im_r = ImageTk.PhotoImage(im_small)

            self.results_listbox.insert(INSERT, self.im_labels[i])
            self.results_listbox.insert(INSERT, '\n')
            self.results_listbox.image_create(INSERT, padx=5, pady=5, image=im_r)
            self.results_listbox.images.append(im_r)
            self.results_listbox.insert(INSERT, '\n')
            
            
main = Main()
main.mainloop()    