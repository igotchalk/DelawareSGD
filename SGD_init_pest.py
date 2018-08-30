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

import SGD
import SGD_writeSWinput
print('Defining model and preparing Seawat files...')
m, ocean_col = SGD_writeSWinput.write_swt_input(modelname)
m.write_ref_file(m.storage_dict)
<<<<<<< HEAD
m.plot_hk_ibound()
=======
>>>>>>> 193c2bd82e3f8606ec55c4e108b3fda97dfdd6b5

################
#3 Run SEAWAT model
################

print('Model size: ' + '{:d}x{:d}x{:d} = {:d}'.format(m.nlay,m.nrow,m.ncol,m.nlay*m.nrow*m.ncol) + ' cells')
print('Stress periods: ' + str(m.nper))
print('Running model...')
import SGD_readrun_model

################
#4 Write .tpl , .ins , .pst
################
print('Writing instruction file...')
ins_data = m.write_ins()

print('Writing TPL file...')

tpl_data = m.write_tpl()

print('Copying output to observation file...')
m.copy_output2obs()

print('Writing control (PST) file..')
m.write_pst(tpl_data,ins_data)

print('...Done!')
