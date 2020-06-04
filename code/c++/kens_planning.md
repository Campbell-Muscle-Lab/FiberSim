## Planning

## FiberSim.exe
+ read an input file
+ launches >=1 muscles as separate thread
+ collates results by lumping FiberSim_results objects into one
+ writes data to file

## Muscle class
+ takes as input
  + FiberSim_model_json_file_string
  + FiberSim_protocol_json_file_string
  + FiberSim_options_json_file_string
+ saves as a FiberSim_data object
  + mechanics data
  + dump files

## FiberSim protocol
+ handles
  + mode
  + delta hsl
  + pCa

## Options
+ Sets a debug mode
  + writes debug information to file in a folder
+ Sets a dump mode
  + dumps map files to a folder

### Optimization

Top-level wrapper calls FiberSim with different input files

### Next steps

+ add in support for mybpc to FiberSim_data
+ add in count mybpc to half-sarcomere class
+ check for mybpc binding
+ then check for sl bumps
