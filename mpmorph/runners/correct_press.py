#! /usr/bin/env python

import glob
from pymatgen.analysis.mddiff import md_data
from pymatgen.analysis.mddiff.correct_volume import CorrectVolume
import numpy as np
import os

dirs = glob.glob("../runs/*/2200/data_1")

print dirs

for dir in dirs:
    print dir
    try:
        os.system("cp -r " + dir + " " + dir + "_1")
    except:
        print ("exception occured")
        continue

    # Find inf the outcar and getting the external pressrure     
    outcar_path = dir + '_1/OUTCAR'
    contcar_path = dir + '_1/CONTCAR'
    print outcar_path, contcar_path
    outcar_data = md_data.get_MD_data(outcar_path)
    outcar_stats = md_data.get_MD_stats(outcar_data)
    external_pressure = outcar_stats[0][0]
    print external_pressure
   
    # create a new poscar to make pressure zero
    corr_vol = CorrectVolume.from_poscar(contcar_path, pressure = external_pressure * 1000, temperature = 2200.0)
    print corr_vol
    corr_vol.rescale_vol()
    corr_vol.poscar.write_file(dir+'_1/POSCAR')
  
