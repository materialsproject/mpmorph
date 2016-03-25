#! /usr/bin/env python

import glob
from pymatgen.analysis.mddiff import md_data
from pymatgen.analysis.mddiff.correct_volume import CorrectVolume
from pymatgen.io.vasp.sets import MITMDVaspInputSet
from pymatgen.io import vasp

import numpy as np
import os

scratch = os.environ['SCRATCH']
old_temp_dirs = glob.glob(scratch + "/diff/runs/*/2200_to_1800/data_1")
new_temp = 1800
old_temp = 1800

vaspmd = MITMDVaspInputSet(old_temp, new_temp, 2000, time_step=2, 
           user_incar_settings={'LWAVE':'.False.', 'LCHARG':'.True.', 'NCORE':24, 'EDIFF':1e-5,
                                'ISIF':1, 'SMASS': 3, 'ENCUT': 400,
                                'SIGMA':0.1}, settings_file = "MDAmorphousVaspInputSet.yaml")
for dir in old_temp_dirs:
    new_dir = dir + "/../../" + str(new_temp)
    contcar_path = dir + '/CONTCAR'
    print new_dir, contcar_path
    contcar = vasp.Poscar.from_file(contcar_path)
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
        os.makedirs(new_dir+"/data_1/")
    
    vaspmd.write_input(contcar.structure, new_dir+"/data_1/") 
