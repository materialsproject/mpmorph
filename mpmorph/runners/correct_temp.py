#! /usr/bin/env python

import glob
from pymatgen.analysis.mddiff.correct_volume import CorrectVolume
from pymatgen.io.vasp.sets import MITMDVaspInputSet
from pymatgen.io import vasp

import numpy as np
import os

scratch = os.environ['SCRATCH']
old_temp_dirs = glob.glob(scratch + "/diff/runs/*/2200/data_2")
new_temp = 1800
old_temp = 2200

vaspmd = MITMDVaspInputSet(old_temp, new_temp, 2000, time_step=2, 
           user_incar_settings={'LWAVE':'.False.', 'NCORE':24, 'EDIFF':1e-5,
                                'ISIF':1, 'SMASS': 3, 'ENCUT': 400})
for dir in old_temp_dirs:
    new_dir = dir + "/../../" + str(old_temp) +"_to_" + str(new_temp)
    xdatcar_path = dir + '/XDATCAR'
    print new_dir, xdatcar_path
    # create a new poscar to approximate_the_new_temp
    xdatcar = vasp.Xdatcar(xdatcar_path)
    corr_vol = CorrectVolume(xdatcar.structures[-1], temperature=old_temp)
    print corr_vol
    corr_vol.rescale_vol(scale='temperature', target_temperature=new_temp)

    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
        os.makedirs(new_dir+"/data_1/")
    
    vaspmd.write_input(corr_vol.structure, new_dir+"/data_1/") 
