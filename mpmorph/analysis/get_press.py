#! /usr/bin/env python

import glob
from pymatgen.analysis.mddiff import md_data
from pymatgen.analysis.mddiff.correct_volume import CorrectVolume
from pymatgen.io.vasp import Poscar
import numpy as np
import os
import matplotlib.pyplot as plt

import itertools
marker = itertools.cycle((',', '+', '.', 'o', '*', 'o', '^')) 

dirs = glob.glob(os.environ['SCRATCH']+"/diff/runs/*/2200/data")
colors = np.random.rand(len(dirs))
print colors
plt.figure()

for di in dirs:
    print di
    p = []
    v = []
    for label in ['_1_1','_1_1_1']:
        d = di +label
        try:
            outcar_data = md_data.get_MD_data(d+"/OUTCAR")[400:]
            poscar_path = d+ "/POSCAR"
        except:
            print "problem with " + d
            continue
        outcar_stats = md_data.get_MD_stats(outcar_data);
        poscar = Poscar.from_file(poscar_path)
        print d[7:14]+d[-5:], outcar_stats[0][0], poscar.structure.volume
        p.append(outcar_stats[0][0])
        v.append(poscar.structure.volume)
    l = sorted( zip(p,v), key = lambda i: i[1])
    vec = [list(i) for i in zip(*l)]
    plt.plot(vec[0],vec[1], linestyle='-', marker=marker.next(),label=di)
plt.show()
 #   external_pressure = outcar_stats[0][0]
 #   print external_pressure
