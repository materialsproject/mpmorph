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
python setup.py develop
```

If you wish to make amorphous structures, please install [packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml) on your machine and add the following line to your bash profile. 
```bash
export PACKMOL_PATH="path_to_packmol_executable_here"
```

## Installation with M3Gnet

This gets tricky. m3gnet is difficult to install on Apple silicon (M1) because of the tensorflow dependency. The most reliable method is as follows:

1) Create and activate a new conda environment/virtualenv. Python 3.10.6 works well.
2) ```conda install -c apple tensorflow-deps```
3) ```pip install tensorflow-macos```
4) Perform a test here by opening a python REPL and running ```import tensorflow```. If you see an issue regarding a numpy datatype not being convertible to a python datatype, reinstall numpy:
```
conda uninstall numpy
conda install numpy
```
5) Install mpmorph
```
pip install mpmorph
```
6) Next, install m3gnet without dependencies (because it will break when it reaches tensorflow)
```
pip install --no-deps m3gnet
```
7) Finally, install the other m3gnet dependencies manually:
```
pip install protobuf==3.19.6 pymatgen ase cython
```
8) Test your installation by opening a python REPL and running:
```import m3gnet```


# Using MPmorph
Before diving headfirst into mpmorph, one should get familiar with how to run fireworks before jumping into mpmorph. Relevant features of fireworks include setting up a fireworks database, adding workflows to the database, configuring a conda environment for running fireworks on a supercomputer, launching jobs to a workload manager (using qlaunch), and monitoring fireworks jobs. For a quick tutorial, check out [beginner](https://www.youtube.com/watch?v=-MalOMJt34U) and [intermediate](https://www.youtube.com/watch?v=zYA_BbKwVO4) atomate/fireworks tutorials from our 2019 Materials project workshop. 

A sample of using mpmorph to run an AIMD simulation at 1500K for 200ps (100k steps at 2fs/step) is shown below:

```python

from mpmorph.workflows.converge import get_converge_wf
from pymatgen.ext.matproj import MPRester

from fireworks import LaunchPad

mpr = MPRester()
structure = mpr.get_structure_by_material_id('mp-1143')
structure.make_supercell([3, 3, 3])

wf = get_converge_wf(structure, temperature = 750, target_steps = 100000)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
```

To generate an amorphous structure, run the following code:

```python
from mpmorph.runners.amorphous_maker import get_random_packed
from mpmorph.workflows.converge import get_converge_wf
from fireworks import LaunchPad

structure = get_random_packed('Li', target_atoms=100)

wf = get_converge_wf(structure, temperature = 5000, target_steps = 10000)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
```

## Customizing runs
At the simplest level of customizing the workflow, one can change the temperature, total number of steps by changing the args passed to get_converge_wf().

For more advanced changes, pass a dictionary for "converge_args", "spawner_args", or "prod_args". This allows customization of vasp inputs, starting and ending temperatures, number of steps, and convergence crieteria.

### Examples:

Changing the rescale parameter, to make for larger volume changes at each stage:
```python
from mpmorph.runners.amorphous_maker import get_random_packed
from mpmorph.workflows.converge import get_converge_wf
from fireworks import LaunchPad

structure = get_random_packed('Li', target_atoms=100)

spawner_args = {'rescale_params': {'beta': 5e-6}}
wf = get_converge_wf(structure, temperature = 5000, target_steps = 10000, spawner_args=spawner_args)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
```

Changing precision of production runs to PREC=Low:
```python
from mpmorph.runners.amorphous_maker import get_random_packed
from mpmorph.workflows.converge import get_converge_wf
from fireworks import LaunchPad

structure = get_random_packed('Li', target_atoms=100)

prod_args = {'optional_fw_params': {'override_default_vasp_params': {'user_incar_settings': {'PREC': 'Normal'}}}}
wf = get_converge_wf(structure, temperature = 5000, target_steps = 10000, prod_args = prod_args)

lp = LaunchPad.auto_load()
lp.add_wf(wf)
```


# Citation

If you use mpmorph, please cite the following papers:
 * Aykol, M., Dwaraknath, S.S., Sun, W. and Persson, K.A., 2018. Thermodynamic limit for synthesis of metastable inorganic materials. Science advances, 4(4), p.eaaq0148.
 * Aykol, M. and Persson, K.A., 2018. Oxidation Protection with Amorphous Surface Oxides: Thermodynamic Insights from Ab Initio Simulations on Aluminum. ACS applied materials & interfaces, 10(3), pp.3039-3045.
