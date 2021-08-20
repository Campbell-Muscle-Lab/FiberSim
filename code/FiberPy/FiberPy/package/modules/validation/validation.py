import os, sys
import json

# Add the validation package

ROOT = os.path.dirname(__file__)
MODULES_ROOT = os.path.realpath(os.path.join(ROOT, "..", ".."))
sys.path.append(MODULES_ROOT)

from modules.validation import myosin_kinetics
from modules.validation import mybpc_kinetics
from modules.validation import actin_kinetics
from modules.validation import force_balance


def run_validation(kinetic_data, batch_file_string):
    """Entrance function for the FiberSim testing suite"""
    
    print("Run_validation starting")
    
    base_folder = os.path.dirname(batch_file_string)
        
    # Get  model/protocol/option files depending on the validation type
    
    val_type = kinetic_data["validation_type"]
    
    if (val_type == 'force_balance'):  # option files is a list for force-balance check
    
        if isinstance(kinetic_data['options_file'], list): # check if a list of options is provided        
            
            dump_folder_list = []
            dump_precision_list = []
        
            for elmt in kinetic_data['options_file']:
                
                if (kinetic_data['relative_to'] == 'this_file'):
                
                    options_file = os.path.join(base_folder, elmt)
                    output_folder = os.path.join(base_folder, kinetic_data['output_data_folder'])
                    
                else:
                    
                    options_file = elmt
                    output_folder = kinetic_data['output_data_folder']                   
                    
                with open(options_file, 'r') as f:
                    opt = json.load(f)
                    
                option_dir = os.path.dirname(options_file)
                    
                if 'status_files' in opt['options']:
                    
                    if opt['options']['status_files']['relative_to'] == 'this_file':
                        dump_folder = os.path.join(option_dir , opt['options']['status_files']['status_folder'])
                    else:
                        dump_folder = opt['options']['status_files']['status_folder']
                                    
                else:
                    raise RuntimeError("No dump folder found to run the test")
                    
                dump_precision_list.append(opt["options"]["x_pos_rel_tol"])
                dump_folder_list.append(dump_folder)
                
        else:
            raise RuntimeError("Option file(s) not provided in a list form")          
            
        # Create output folder if it does not exist
        
        if not os.path.exists(output_folder):
            os.makedirs(output_folder) 
                                    
        print("force_balance check")
        
        force_balance.compute_force_balance(dump_precision_list, dump_folder_list, output_folder) 
        
    else:  # it is a kinetic test    
    
        if (kinetic_data['relative_to'] == 'this_file'):  
            
            model_file = os.path.join(base_folder, kinetic_data['model_file'])
        
            protocol_file = os.path.join(base_folder, kinetic_data['protocol_file'])
            
            options_file = os.path.join(base_folder, kinetic_data['options_file'])
            
            output_folder = os.path.join(base_folder, kinetic_data['output_data_folder'])
            
        else:            
            model_file = kinetic_data['model_file']
            protocol_file = kinetic_data['protocol_file']
            options_file = kinetic_data['options_file']
            output_folder = kinetic_data['output_data_folder']
            
        # Create output folder if it does not exist

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)   
                                       
        # Get dump_folder from option file
        
        with open(options_file, 'r') as f:
            opt = json.load(f)
        
        option_dir = os.path.dirname(options_file)
                
        if 'status_files' in opt['options']:
            if opt['options']['status_files']['relative_to'] == 'this_file':
                dump_folder = os.path.join(option_dir, opt['options']['status_files']['status_folder'])
            else:
                dump_folder = opt['options']['status_files']['status_folder']
                
        else:
            raise RuntimeError("No dump folder found to run the test")
            
        # Get number of adjacent bs from option file
        
        if 'adjacent_bs' in opt['options']:
            adj_bs = opt["options"]["adjacent_bs"]        
        else:
            adj_bs = 0
        
        # Now run the kinetics test
        
        val_type = kinetic_data["validation_type"]
        
        if (val_type == 'm_kinetics'): 
            
            print("Myosin kinetics check")
            
            myosin_kinetics.compute_rate(model_file, protocol_file, dump_folder, output_folder, adj_bs)
            
        if (val_type == 'c_kinetics'): 
                  
            print("Mybpc kinetics check")
            
            mybpc_kinetics.compute_rate(model_file, protocol_file, dump_folder, output_folder, adj_bs)
            
        if (val_type == 'a_kinetics'): 
                    
            print("Actin kinetics check")
            
            actin_kinetics.compute_rate(model_file, protocol_file, dump_folder, output_folder)   
            

        
        