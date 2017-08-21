# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This Drone tries to produce a more sensible task dictionary than the default VaspToDbTaskDrone.
Some of the changes are documented in this thread:
https://groups.google.com/forum/#!topic/pymatgen/pQ-emBpeV5U
"""

import os
import re
import datetime
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob

from monty.io import zopen

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.io.vasp.inputs import Poscar, Potcar, Incar, Kpoints
from atomate.vasp.drones import VaspDrone

from matgendb.creator import get_uri

from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew, Shyue Ping Ong, Shyam Dwaraknath, Anubhav Jain'
__email__ = 'kmathew@lbl.gov'
__date__ = 'Mar 27, 2016'
__version__ = "0.1.0"

logger = get_logger(__name__)

class VaspMDDrone(VaspDrone):
    schema = {
        "root": {
            "schema", "dir_name", "chemsys", "composition_reduced",
            "formula_pretty", "formula_reduced_abc", "elements",
            "nelements", "formula_anonymous", "calcs_reversed", "completed_at",
            "nsites", "composition_unit_cell", "input", "output", "state",
            "analysis", "run_stats"
        },
        "input": {'is_lasph', 'is_hubbard', 'xc_override', 'potcar_spec',
                  'hubbards', 'structure', 'pseudo_potential'},
        "output": {'structure', 'spacegroup', 'density', 'energy',
                   'energy_per_atom', 'is_gap_direct', 'bandgap', 'vbm',
                   'cbm', 'is_metal'},
        "calcs_reversed": {
            'dir_name', 'run_type', 'elements', 'nelements',
            'formula_pretty', 'formula_reduced_abc', 'composition_reduced',
            'vasp_version', 'formula_anonymous', 'nsites',
            'composition_unit_cell', 'completed_at', 'task', 'input', 'output',
            'has_vasp_completed',
        },
        "analysis": {'delta_volume_percent', 'delta_volume', 'max_force',
                     'errors',
                     'warnings'}
    }
    def process_vasprun(self, dir_name, taskname, filename):
        """
        Adapted from matgendb.creator
        Process a vasprun.xml file.
        """
        vasprun_file = os.path.join(dir_name, filename)
        if self.bandstructure_mode:
            vrun = Vasprun(vasprun_file, parse_eigen=True, parse_projected_eigen=True)
        else:
            vrun = Vasprun(vasprun_file)

        d = vrun.as_dict()

        # rename formula keys
        for k, v in {"formula_pretty": "pretty_formula",
                     "composition_reduced": "reduced_cell_formula",
                     "composition_unit_cell": "unit_cell_formula"}.items():
            d[k] = d.pop(v)

        for k in ["eigenvalues", "projected_eigenvalues"]:  # large storage space breaks some docs
            if k in d["output"]:
                del d["output"][k]

        comp = Composition(d["composition_unit_cell"])
        d["formula_anonymous"] = comp.anonymized_formula
        d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = str(datetime.datetime.fromtimestamp(os.path.getmtime(vasprun_file)))
        d["density"] = vrun.final_structure.density

        # replace 'crystal' with 'structure'
        d["input"]["structure"] = d["input"].pop("crystal")
        d["output"]["structure"] = d["output"].pop("crystal")
        for k, v in {"energy": "final_energy", "energy_per_atom": "final_energy_per_atom"}.items():
            d["output"][k] = d["output"].pop(v)

        if self.parse_dos and self.parse_dos != 'final':
            try:
                d["dos"] = vrun.complete_dos.as_dict()
            except:
                raise ValueError("No valid dos data exist in {}.".format(dir_name))

        if self.bandstructure_mode:
            bs = vrun.get_band_structure(line_mode=(self.bandstructure_mode.lower() == "line"))
        else:
            bs = vrun.get_band_structure()

        d["bandstructure"] = bs.as_dict()

        d["output"]["vbm"] = bs.get_vbm()["energy"]
        d["output"]["cbm"] = bs.get_cbm()["energy"]
        bs_gap = bs.get_band_gap()
        d["output"]["bandgap"] = bs_gap["energy"]
        d["output"]["is_gap_direct"] = bs_gap["direct"]
        d["output"]["is_metal"] = bs.is_metal()
        d["task"] = {"type": taskname, "name": taskname}

        if hasattr(vrun, "force_constants"):
            # phonon-dfpt
            d["output"]["force_constants"] = vrun.force_constants.tolist()
            d["output"]["normalmode_eigenvals"] = vrun.normalmode_eigenvals.tolist()
            d["output"]["normalmode_eigenvecs"] = vrun.normalmode_eigenvecs.tolist()
        return d
