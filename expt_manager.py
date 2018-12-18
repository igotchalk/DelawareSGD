import os
import numpy as np
import flopy
import scipy.stats as stats
from pathlib import Path
from SGD import ModelSGD
'''
This class is for managing Monte Carlo experiments. Built to be used with
Seawat FloPy simulations.
Initiated Dec 1, 2018
'''

class Expt(ModelSGD):

    #pars_available = namedtuple("Pars", ["left", "right"])

    '''
    This class can load, create, and save experiments for running
    '''
    def __init__(self,filename= None, expt_ws = None, expt_dict=dict()):
        super(ModelSGD, self).__init__()  #inherit ModelSGD properties
        self.filename = filename
        self.expt_dict = expt_dict
        print(self.expt_dict)
        self.filename = filename
        self.expt_ws = expt_ws
        if self.expt_ws is None and self.filename is not None:
            self.expt_ws = Path(filename).parent

    @classmethod
    def SGD2expt(self, objSGD):
        objSGD.load
        objSGD.__class__ = Expt
        objSGD.__init__() #initialize to get all default properties
        return

    def add_par(self,parname,packname,bounds):

        #right now uniform sampling between bounds
        pardict = {'parname':parname,
                   'packname':packname,
                   'bounds':bounds}
        self.expt_dict[parname] = pardict
        return

    def load(self,filename=None):
        #Navigate to file, load experiment data stored in directory
        if filename is None:
            filename=self.filename
        return