
# coding: utf-8
import flopy
import os
import sys
import config
import SGD
import re
import ast
from pathlib import Path,PureWindowsPath

#Load ref file to get info about model
reffile = os.path.join('.','ref_file.txt')
reftext = open(reffile, 'r').read()    
beg = [m.start() for m in re.finditer('<<<', reftext)]
betw = [m.start() for m in re.finditer('>>>', reftext)]
end = [m.start() for m in re.finditer('\n', reftext)]

#define model_ws and modelname in the ref file
for i in range(len(beg)):
    exec(reftext[beg[i]+3:betw[i]] + ''' = "''' + reftext[betw[i]+3:end[i]] + '''"''')
print(model_ws)
ws = Path(model_ws)
fname = [f for f in ws.iterdir() if f.suffix == '.nam'][0]
m = flopy.seawat.Seawat.load(str(fname),exe_name = config.swexe, model_ws = str(ws))
SGD.ModelSGD.Seawat2SGD(m)
#Run model
v = m.run_model(silent=True, report=True)

#Write output
ocean_col = ast.literal_eval(ocean_col) #covert ocean_col from string to list
d = {'ocean_col':ocean_col}
m.set_storage_dict(d)
m.write_output()
