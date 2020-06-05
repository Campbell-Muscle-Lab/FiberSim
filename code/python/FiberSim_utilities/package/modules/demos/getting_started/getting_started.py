"""
Entry point for FiberSim_utilities
@author: Ken Campbell
"""

from package.modules.demos.getting_started.single_run.single_run \
    import demo_single_run
from package.modules.demos.getting_started.single_run_with_log.single_run_with_log \
    import demo_single_run_with_log
from package.modules.demos.getting_started.multiple_run.multiple_run \
    import demo_multiple_run

def parse_inputs(inputs):
    # Get the number of arguments
    no_of_arguments = len(inputs)
    
    if (no_of_arguments > 2):
        if (inputs[3] == "single_run"):
            demo_single_run()
        if (inputs[3] == "single_run_with_log"):
            demo_single_run_with_log()            
        if (inputs[3] == "multiple_run"):
            demo_multiple_run()            

if __name__ == "__main__":
    parse_inputs()
           
    



# import sys
# import subprocess

# import numpy as np
# import json

# from modules.FiberSim_analysis.FiberSim_analysis import FiberSim_analysis as FiberSim_anal
# from modules.half_sarcomere import half_sarcomere as hs

# from modules.protocol.protocol import write_protocol_to_file
# from modules.batch.batch import run_batch
# from modules.optimize.optimize import run_optimization

# if __name__ == "__main__":
    
#     # Get the number of arguments
#     no_of_arguments = len(sys.argv)
    
#     if (no_of_arguments == 1):
        
#         # run_batch("..\\demo_files\\test_stretch\\batch\\batch_known_single.json")
        
#         # dfs = "..\\demo_files\\test_stretch\\results\\10\\results.txt"
#         # fs_anal = FiberSim_anal()
#         # fs_anal.display_data(dfs)
        
#         ofs = "..\\demo_files\\optimization_files\\known_optimization.json"
#         run_optimization(ofs)
        
        
#         # fs_anal = FiberSim_anal()
#         # dfs = "c:\\temp\\log_58b\\hs_status"
#         # fs_anal.animate_filament_states(dfs,
#         #                                 "c:\\temp\\thick_states.gif")
        
#         # subprocess.call("C:\\ken\\GitHub\\CampbellMuscleLab\\Models\\FiberSim_VS2019\\C++\\FiberSim\\x64\Release\\FiberSim.exe")
        
        
#         # fs_anal = FiberSim_anal();
#         # dfs = "C:\\temp\\ken.txt"
#         # fs_anal.display_data(dfs)
        
#         # sfs = "c:\\temp\\log\\hs_status\\hs_0_time_step_249.json"
        
#         # hs_0 = hs.half_sarcomere(sfs)
#         # hs_0.draw_filaments()
        
#         # fs_anal.animate_cb_distributions("c:\\temp\\log\\hs_status",
#         #                                  "c:\\temp\cb_dist.gif")
        
#         # hs_0 = hs.half_sarcomere(sfs)
        
#         # hs_0.draw_cb_distributions("c:\\temp\\cb_dist.png")
        
        
#         # hs_0.draw_filament_lengths()
        
#         # fs_anal.animate_filament_lengths("c:\\temp\log\\hs_status",
#         #                                  "c:\\temp\\test.gif",
#         #                                   frames=[])
        
#         """print("FiberSim_utilities called with no inputs")
#         print("Displaying data")
#         sfs = "c:\\temp\\log\\hs_status\\hs_0_time_step_4.json"
#         fs_anal.load_hs_status_file(sfs)"""
        
#     if (no_of_arguments == 2):
#         if (sys.argv[1] == 'display_FiberSim_data'):
#             print("Displaying data")
#             dfs = "C:\\temp\\pCa_test\\result_pCa45.txt"
#             fs_a = FiberSim_anal()
#             fs_a.display_data(dfs,output_file_string='c:\\temp\\test_display.png')

#         if (sys.argv[1] == "generate_protocol"):
#             pfs = "C:\\ken\\GitHub\\CampbellMuscleLab\\Models\\FiberSim_VS2019\\demo_files\\protocol_files\\pCa"
#             for pCa_value in [9.0, 6.5, 6.3, 6.2, 6.1, 6.0, 5.9, 5.8, 5.7, 5.5, 4.5]:
#                 print([pCa_value])
#                 fs = ("%s_%.0f_protocol.txt") % (pfs, 10*pCa_value)
#                 print(fs)
#                 write_protocol_to_file(fs, pCa = [pCa_value])
        
#         if (sys.argv[1] == "generate_protocol_2"):
#             pfs = "C:\\ken\\GitHub\\CampbellMuscleLab\\Models\\FiberSim_VS2019\\demo_files\\protocol_files\\pCa58_step.txt"
#             dhsl = np.zeros(2000)
#             dhsl[1300] = 15.0
#             write_protocol_to_file(pfs, pCa = [5.8], dhsl = dhsl)

#     if (no_of_arguments == 3):
#         if (sys.argv[1] == 'run_batch'):
#             run_batch(sys.argv[2])
#             # print("Running batch: %s" % sys.argv[2])
#             # with open(sys.argv[2]) as json_file:
#             #     json_data = json.load(json_file)
#             #     FiberSim_batch = json_data['FiberSim_batch']
#             #     job_data = FiberSim_batch['job']
#             #     for i, j in enumerate(job_data):
#             #         print('Running job %i' % (i))
#             #         model_file_string = j['model_file_string']
#             #         options_file_string = j['options_file_string']
#             #         protocol_file_string = j['protocol_file_string']
#             #         results_file_string = j['results_file_string']

#             #         # Generate a command
#             #         exe_string = "C:\\ken\\GitHub\\CampbellMuscleLab\\Models\\FiberSim_VS2019\\C++\\FiberSim\\x64\Release\\FiberSim.exe"
#             #         command_string = "%s %s %s %s %s" % (exe_string,
#             #                                           model_file_string,
#             #                                           options_file_string,
#             #                                           protocol_file_string,
#             #                                           results_file_string)
                    
#             #         subprocess.call(command_string)
                    
#         if (sys.argv[1] == 'display_pCa_data'):
#             print("Displaying pCa data for: %s" % sys.argv[2])
#             fs_a = FiberSim_anal()
#             fs_a.display_pCa_data(sys.argv[2],
#                                   output_file_string="c:\\temp\\f_test")
