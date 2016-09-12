# MPmorph

MPmorph is a collection of tools to run and analyze ab-initio molecular dynamics (AIMD) calculations run with VASP, 
and is currently under development.
It relies heavily on tools developed by the Materials Project ([pymatgen](http://www.pymatgen.org), 
[custodian](https://github.com/materialsproject/custodian), 
[fireworks](https://github.com/materialsproject/fireworks)) and [MatMethods](https://github.com/hackingmaterials/MatMethods).

MPmorph provides:
* Infrastructure for dynamic VASP MD workflows:
  * Tools to create dynamic MD workflows using MatMethods
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