#!/usr/bin/python3 -i
'''
Initialization python file ot run PEST
  1) takes an input model name
  2) creates SEAWAT input files
  3) runs SEAWAT 1x
  4) creates:
     template file
     instruction filef
     control file
'''
################
#1 Take input model name from cmd
################

import sys
print('Running ' + sys.argv[0])
modelname = sys.argv[1]

################
#2 Create SEAWAT input files
################

import SGD_writeSWinput
m, ocean_col = SGD_writeSWinput.write_swt_input(modelname)
m.write_ref_file()

################
#3 Run SEAWAT model
################

import SGD_readrun_model

################
#4 Write .tpl , .ins , .pst
################
ins_data = m.write_ins()
tpl_data = m.write_tpl()
m.copy_output2obs()
m.write_pst(tpl_data,ins_data)
#import SGD_write_pest_files
#ins_data = SGD_write_pest_files.write_ins(m)
#tpl_data = SGD_write_pest_files.write_tpl(m)
#SGD_write_pest_files.write_pst(m,tpl_data,ins_data)

##Plot the initial K-field
#SGD_writeSWinput.plot(m)
