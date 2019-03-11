# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 14:42:39 2017

@author: ianpg
"""
import os
dir_path = os.path.dirname(os.path.realpath(get_params.py))
def listparams(filename):
    tempfile = cd
    with open(filename,'r') as f, open(tempfile,'w') as pfile:
        for line in f:
            if line=='#Parameters\n':
                print('start recording parameters')
                params = ''
            elif line=='#/Parameters\n':
                print('end of parameter list')
                break
            else:
                params = params + line.strip() + '\n'
    print(params)

