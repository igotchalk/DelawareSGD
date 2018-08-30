
# coding: utf-8
import flopy
import os
import sys
import config
import SGD
import re
import ast
from pathlib import Path,PureWindowsPath
import utils


d = utils.read_ref()

ws = Path(d['model_ws'])
fname = [f for f in ws.iterdir() if f.suffix == '.nam'][0]
m = flopy.seawat.Seawat.load(str(fname),exe_name = config.swexe, model_ws = str(ws.as_posix()))
SGD.ModelSGD.Seawat2SGD(m)
#Run model
v = m.run_model(silent=True, report=True)

#Write output
#ocean_col = ast.literal_eval(d['ocean_col']) #covert ocean_col from string to list
#m.set_storage_dict(d)
m.write_output()
