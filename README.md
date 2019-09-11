# MPmorph

MPmorph is a collection of tools to run and analyze ab-initio molecular dynamics (AIMD) calculations run with VASP, 
and is currently under development.
It relies heavily on tools developed by the Materials Project ([pymatgen](http://www.pymatgen.org), 
[custodian](https://github.com/materialsproject/custodian), 
[fireworks](https://github.com/materialsproject/fireworks)) and [atomate](https://github.com/hackingmaterials/atomate).

MPmorph provides:
* Infrastructure for dynamic VASP MD workflows:
  * Tools to create dynamic MD workflows using atomate
  * E.g. generation of new MD runs based on the given criterion (currently pressure, for liquid/amorphous phase density estimation)
* Tools for statistical analysis of:
  * Static observables:
    * Radial distribution functions
    * Coordination numbers
    * Voronoi analysis
    * Polyhedra connectivities
    * Thermodynamic quantities
  * Diffusion coefficients:  
    * Robust calculation of diffusion coefficients (D) and activation energies (Q).
    * Rigorous error analysis for D and Q.

# Installation

Before installing mpmorph, install the latest version of [pymatgen](http://www.pymatgen.org), 
[custodian](https://github.com/materialsproject/custodian), 
[fireworks](https://github.com/materialsproject/fireworks)) and [atomate](https://github.com/hackingmaterials/atomate)

clone the repository to your computer and install using 
```bash
python setup.py install
```

If you wish to make amorphous structures, please install [packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml) on your machine and add the following line to your bash profile
```bash
export PACKMOL_PATH="path_to_packmol_here"
```

# Using MPmorph

A sample of using mpmorph to run an AIMD simulation at 1500K for 200ps (100k steps at 2fs/step) is shown below:

```python
from mpmorph.workflow.converge import get_converge_wf
from pymatgen import MPRester
from fireworks import LaunchPad

mpr = MPRester()
structure = mpr.get_structure_by_material_id('mp-1143')
structure.make_supercell([3, 3, 3])

wf = get_converge_wf(structure, temperature = 1500, max_steps = 5000, target_steps = 100000)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
```

To generate an amorphous structure, run the following code:

```python
from mpmorph.runners.amorphous_maker import get_random_packed
from mpmorph.workflow.converge import get_converge_wf
from pymatgen import MPRester
from fireworks import LaunchPad

structure = get_random_packed('Li', target_atoms=100)

wf = get_converge_wf(structure, temperature = 5000, max_steps = 5000, target_steps = 10000)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
```

# Citation

If you use mpmorph, please cite the following papers:
 * Aykol, M., Dwaraknath, S.S., Sun, W. and Persson, K.A., 2018. Thermodynamic limit for synthesis of metastable inorganic materials. Science advances, 4(4), p.eaaq0148.
 * Aykol, M. and Persson, K.A., 2018. Oxidation Protection with Amorphous Surface Oxides: Thermodynamic Insights from Ab Initio Simulations on Aluminum. ACS applied materials & interfaces, 10(3), pp.3039-3045.
