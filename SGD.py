import os
from pathlib import Path,PureWindowsPath
import numpy as np
import flopy
from flopy.seawat import Seawat
import sys
import utils

#Subclass of flopy.seawat.Seawat
#    Carries relevant info about the ocean boundary and other info
class ModelSGD(Seawat):
    def __init__(self,ocean_arr = None,storage_dict=None, ocean_col=[30,69]):
        super(Seawat, self).__init__()  #inherit Seawat properties
        
        #Set other properties
        self.ocean_arr = ocean_arr
        self.storage_dict = storage_dict
        self.ocean_col = ocean_col
        print("made an SGD model!")
    #Take a Seawat object and turn it into an SGD object
    @classmethod
    def Seawat2SGD(self, objSeawat):
        objSeawat.load
        objSeawat.__class__ = ModelSGD
        return
    def set_ocean_arr(self,arr):
        self.ocean_arr = arr
        return
    def set_storage_dict(self,storage_dict):
        self.storage_dict = storage_dict
        return
    def write_output(self,fname='flux.smp'):
        #Get flux at ocean boundary
        try:
            ocean_col = self.storage_dict['ocean_col']
        except KeyError:
            ocean_col = self.ocean_col
        ocean_outflow = utils.get_ocean_outflow_chd(self,ocean_col)
        ocean_col_vec = np.arange(ocean_col[0],ocean_col[1]+2)
        print(type(ocean_outflow))
        print(ocean_outflow)
        #Print out coordinates and flux to text file
        fout= open(os.path.join(self.model_ws,fname),"w")
        fout.write('Values are zero-based \n')
        fout.write('{:14s} {:4s} {:4s} {:4s} \n'
                   .format("flux", "lay","row", "col"))
        for i in range(len(ocean_outflow)):
             fout.write('{:14.4e} {:4d} {:4d} {:4d}\n'
                        .format(ocean_outflow[i],0,0,ocean_col_vec[i]))
        fout.close()
        print('output FILE WRITTEN: ' + os.path.join(self.model_ws, fname))
        return
    
    #Copy the output file to create a synthetic observation file
    def copy_output2obs(self,suffix='.smp'):
        import shutil
        
        def get_fname(model_ws,ext):
            fname = [os.path.join(model_ws,f) for f in os.listdir(model_ws) if f.endswith(ext)][0];
            return fname
        f_output = get_fname(self.model_ws,suffix)
        f_obs = f_output[:-4] + '.obs'
        shutil.copy2(f_output,f_obs)
        print('observation file copied from: ' + f_output + '\nTo: ' + f_obs)
        
    #Write instruction file
    def write_ins(self,fname = 'flux.ins',obs_name = 'flux'):
        #Get flux at ocean boundary so you know how many lines to write
        ocean_col = self.storage_dict['ocean_col']
        ocean_col_vec = np.arange(ocean_col[0],ocean_col[1]+1)
        ocean_outflow = utils.get_ocean_outflow_chd(self,ocean_col)

        #Write an instruction file
        finst = open(os.path.join(self.model_ws,fname),"w")
        finst.write('pif #\n')
        finst.write('#flux#\n')
        for i in range(len(ocean_outflow)):
             finst.write('l1 w !{:s}!\n'.format(obs_name + '_' + str(i)))
        finst.close()
        print('.ins FILE WRITTEN: ' + os.path.join(self.model_ws, fname))
        nobs = len(ocean_outflow)
        ins_data = [obs_name,nobs,ocean_outflow]
        return ins_data
    
    #Make template file
    def write_tpl(self,zonearray=None, parzones = [2],lbound = 0.001,
                  ubound = 1000.,transform='log'):
        mfpackage = 'lpf'
        partype = 'hk'
        if zonearray == None:
            zonearray = np.ones((self.nlay, self.nrow, self.ncol), dtype=int)
            zonearray[2] = 2 #make layer 3 the zone (zone # 2) that will be parameterized
        lpf = self.get_package(mfpackage)
        parvals = [np.mean(lpf.hk.array)]
        
        plist = flopy.pest.zonearray2params(mfpackage, partype, parzones, lbound,
                                              ubound, parvals, transform, zonearray)
        tpl_data = [mfpackage, partype, parzones, lbound,
                                              ubound, parvals, transform, zonearray]
        # Write the template file
        tw = flopy.pest.templatewriter.TemplateWriter(self, plist)
        tw.write_template()
        print('.tpl FILE WRITTEN: ' + os.path.join(self.model_ws, self.name +'.' + mfpackage + '.tpl'))
        npar = len(parzones)
        return tpl_data
    
    def write_pst(self,tpl_data,ins_data):
        import pandas

        fname = self.name + '.pst'
        f = open(os.path.join(self.model_ws,fname),"w")
        f.write('pcf\n')

        #Control data:
        '''
        RSTFLE PESTMODE
        NPAR NOBS NPARGP NPRIOR NOBSGP [MAXCOMPDIM] [DERZEROLIM]
        NTPLFLE NINSFLE PRECIS DPOINT [NUMCOM JACFILE MESSFILE] [OBSREREF]
        RLAMBDA1 RLAMFAC PHIRATSUF PHIREDLAM NUMLAM [JACUPDATE] [LAMFORGIVE] [DERFORGIVE]
        RELPARMAX FACPARMAX FACORIG [IBOUNDSTICK UPVECBEND] [ABSPARMAX]
        PHIREDSWH [NOPTSWITCH] [SPLITSWH] [DOAUI] [DOSENREUSE] [BOUNDSCALE]
        NOPTMAX PHIREDSTP NPHISTP NPHINORED RELPARSTP NRELPAR [PHISTOPTHRESH] [LASTRUN] [PHIABANDON]
        ICOV ICOR IEIG [IRES] [JCOSAVE] [VERBOSEREC] [JCOSAVEITN] [REISAVEITN] [PARSAVEITN] [PARSAVERUN]
        '''
        npar = len(tpl_data[2]) #length of parzones
        nobs = ins_data[1]
        npargp = 1 #number of param groups
        nprior = 0 #number of articles of prior info
        nobsgp = 1 #number of obs groups

        ntplfle = 1 #num tpl files
        ninsfle = 1 #num ins files
        precis = 'double'
        dpoint = 'point'

        rlambda1 = 10.0
        rlamfac = -3.0
        phiratsuf = 0.3
        phiredlam = 0.03
        numlam = 10

            #im not writing any more of these variables for now...
        f.write('* control data\n')
        f.write('restart estimation\n')
        f.write('{:d} {:d} {:d} {:d} {:d}\n'
                .format(npar,nobs,npargp,nprior,nobsgp))
        f.write('{:d} {:d} {:s} {:s}\n'
                .format(ntplfle,ninsfle,precis,dpoint))
        f.write('{:f} {:f} {:f} {:f} {:d}\n'
                .format(rlambda1,rlamfac,phiratsuf,phiredlam,numlam))
        f.write('10.0  10.0  0.001\n')
        f.write('0.1\n')
        f.write('30  0.005  4  4  0.005  4  1E-5\n')
        f.write('1  1  1\n')


        #LSQR
        '''
        LSQRMODE
        LSQR_ATOL LSQR_BTOL LSQR_CONLIM LSQR_ITNLIM
        LSQRWRITE
        '''   
        #f.write('* lsqr\n')

        #parameter groups
        '''
        PARGPNME INCTYP DERINC DERINCLB FORCEN DERINCMUL DERMTHD [SPLITTHRESH SPLITRELDIFF SPLITACTION]
        '''
        f.write('* parameter groups\n')
        f.write('hk relative 0.01 0.0 switch 2.0 parabolic\n')

        #parameter data
        '''
        PARNME PARTRANS PARCHGLIM PARVAL1 PARLBND PARUBND PARGP SCALE OFFSET DERCOM 
        '''
        parname = tpl_data[1] + '_' + str(tpl_data[2][0])
        partrans = 'log'
        parchglim = 'factor'
        parval1 = tpl_data[5][0]
        parlbnd = tpl_data[3]
        parubnd = tpl_data[4]
        parg = 'hk'
        scal = 1.0
        offset = 0.0
        dercom = 1

        f.write('* parameter data\n')
        f.write('{:s} {:s} {:s} {:f} {:e} {:e} {:s} {:f} {:f} {:d}\n'.
               format(parname,partrans,parchglim,parval1,parlbnd,parubnd,parg,scal,offset,dercom))

        #observation groups
        '''
        OBGNME [GTARG] [COVFLE]
        '''
        obgnme = 'heads'
        f.write('* observation groups\n')
        f.write('{:s}\n'.format(obgnme))
        
        #observation data
        import pandas
        
        def get_fname(model_ws,ext):
            fname = [os.path.join(model_ws,f) for f in os.listdir(model_ws) if f.endswith(ext)][0];
            return fname
        
        #"True" observations: should be same format as the .ins file and the .pts file
        fname_obs = get_fname(self.model_ws,'.obs')
        obs_array = pandas.read_csv(fname_obs,delim_whitespace=True,skiprows=1)
        obsval = obs_array.loc[:,'flux']
        weight = 10

        f.write('* observation data\n')    
        for i in range(nobs):
            obsname = ins_data[0] + '_' + str(i)
            f.write('{:s} {:f} {:f} {:s}\n'.format(obsname,obsval[i],weight,obgnme))

        #model command line
        f.write('* model command line\n')
        f.write('SGD_runmodel.bat\n')

        #model input/output
        '''
        TEMPFLE INFLE
        INSFLE OUTFLE
        '''


        tempfle = get_fname(self.model_ws,'.tpl')
        infle = get_fname(self.model_ws,'.' + tpl_data[0])
        insfle = get_fname(self.model_ws,'.ins')
        outfle = get_fname(self.model_ws,'.smp')

        f.write('* model input/output\n')
        f.write('{:s} {:s}\n'.format(tempfle,infle))
        f.write('{:s} {:s}\n'.format(insfle,outfle))

        #DONE
        f.close()
        print('.pst FILE WRITTEN: ' + os.path.join(self.model_ws, fname))
        return

    def write_ref_file(self,d=None):
        from pathlib import Path
        if d==None:
            d = {'model_ws':Path(self.model_ws),
            'modelname': self.name,
            'ocean_col': self.storage_dict['ocean_col']
            }
        fname = Path(os.path.join(self.model_ws,'..','..','ref_file.txt'))
        #write the file
        fo = open(fname, "w")
        for k, v in d.items():
            fo.write('<<<' + str(k) + '>>>'+ str(v) + '\n')
        fo.close()
        print('reference FILE WRITTEN: ' + str(fname))
        return