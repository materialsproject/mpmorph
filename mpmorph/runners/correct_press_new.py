#! /usr/bin/env python

import glob
from pymatgen.analysis.mddiff import md_data
from pymatgen.analysis.mddiff.correct_volume import CorrectVolume
from pymatgen.io.vasp.sets import MITMDVaspInputSet
from pymatgen.io import vasp

import numpy as np
import os

scratch = os.environ['SCRATCH']
old_pres_dirs = glob.glob(scratch + "/diff/runs/*/1800/")
new_temp = 1800
old_temp = 1800
initial_offset = 0.5 # skips this fraction of data points


vaspmd = MITMDVaspInputSet(old_temp, new_temp, 2000, time_step=2, 
           user_incar_settings={'LWAVE':'.False.', 'LCHARG':'.True.', 'NCORE':24, 'EDIFF':1e-5,
                                'ISIF':1, 'SMASS': 3, 'ENCUT': 400,
                                'SIGMA':0.1}, settings_file = "MDAmorphousVaspInputSet.yaml")

for dir in old_pres_dirs:
    # check if average pressure is below 5 kb
    # a) if one data point -> extrapolate from compressibility
    # b) if two data points -> extrapolate  from linear fit
    # c) if three data points -> fit polynomial
    # Check if the new predicted volume is more than 30% different than the initial volumes

    prev_press = glob.glob(dir+"/data_1")
    n_dirs = len(prev_press)   
    # folder naming convention : _1, _1_1, _1_1_1

    data_dir_dict = {}

    for dir in prev_press:
        outcar_path = dir+"/OUTCAR"
        print "reading outcar"
        outcar = vasp.outputs.Outcar(outcar_path)
        print "done"
        outcar.read_pattern({"pressure": "external\s+pressure\s+=\s+([\d\-\.]+)"})
        p = np.array([ float(i[0]) for i in outcar.data['pressure']])
        p_avg = p[ int(len(p)*initial_offset): ].mean() 
        print p_avg 
        data_dir_dict[dir] = {}
        data_dir_dict[dir]['p_avg'] = p_avg
        data_dir_dict[dir]['outcar'] = outcar
     

    if n_dirs == 1:      
        dir = prev_press[0]
        new_dir = dir + "_1"
        xdatcar_path = dir + '/XDATCAR'
        print new_dir, xdatcar_path
        # create a new poscar to approximate_the_new_temp
        xdatcar = vasp.Xdatcar(xdatcar_path)
        corr_vol = CorrectVolume(xdatcar.structures[-1], pressure=data_dir_dict[dir]['p_avg']*1000, temperature=old_temp)
        print corr_vol
        corr_vol.rescale_vol(scale='pressure', target_pressure=0.0)

        vaspmd.write_input(corr_vol.structure, new_dir) 
