#! /usr/bin/env python

import glob
from pymatgen.analysis.mddiff import md_data
import matplotlib.pyplot as plt
import numpy as np
import os


scratch = os.environ['SCRATCH']
outcars = glob.glob(scratch + "/diff/runs/*/2200/data_2/OUTCAR")
#outcars = glob.glob("../Cr28Fe4O48/2600/data/OUTCAR")
print outcars

plt.figure()
for outcar in outcars:
#    plt.figure()
    print outcar
    data = md_data.get_MD_data(outcar,
#            search_keys = ['kinetic energy EKIN', 'ETOTAL'], search_data_column = [4, 4]
           )
    print md_data.get_MD_stats(data[1:])
    data = np.array(data)
    plt.plot(data[:,0])
    plt.ylabel("Pressure (kb)")
    plt.xlabel("MD step")
plt.show()

