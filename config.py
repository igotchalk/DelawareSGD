# setup absolute paths to the executables based on the OS

import sys
import os

exeext = ''

if sys.platform.lower() == 'darwin':
    bindir = os.path.abspath('/Users/ianpg/Documents/ProjectsLocal/SyntheticSWI/bin')
    exepth = os.path.join(bindir, 'modflow_macOS')
elif 'win' in sys.platform.lower():
    bindir = os.path.join('E:\Projects\DelawareSGD','bin')
    exeext = '.exe'
    winarch = 'win64'  #change to 'win32' for a 32-bit system
    exepth = os.path.join(bindir, winarch)
else:
    raise Exception('Could not find binaries for {}'.format(sys.platform))    

mfexe = os.path.abspath(os.path.join(exepth, 'mf2005{}'.format(exeext)))
mpexe = os.path.abspath(os.path.join(exepth, 'mp6{}'.format(exeext)))
mtexe = os.path.abspath(os.path.join(exepth, 'mt3dms{}'.format(exeext)))
swexe = os.path.abspath(os.path.join(exepth, 'swt_v4{}'.format(exeext)))

mf6exe =       os.path.abspath(os.path.join(exepth, 'mf6{}'.format(exeext)))
mf2000exe =    os.path.abspath(os.path.join(exepth, 'mf2000{}'.format(exeext)))
mflgrexe =     os.path.abspath(os.path.join(exepth, 'mflgr{}'.format(exeext)))
mfnwtexe =     os.path.abspath(os.path.join(exepth, 'mfnwt{}'.format(exeext)))
mfusgexe =     os.path.abspath(os.path.join(exepth, 'mfusg{}'.format(exeext)))
mt3dusgsexe =  os.path.abspath(os.path.join(exepth, 'mt3dusgs{}'.format(exeext)))

exelist = [mfexe, mpexe, mtexe, swexe, mf6exe, mf2000exe, mflgrexe,
           mfnwtexe, mfusgexe, mt3dusgsexe]
for e in exelist:
    if not os.path.isfile(e):
        pass
        #print('Executable file could not be found: {}'.format(e))
