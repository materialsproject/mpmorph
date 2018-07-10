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
    
# Citation

If you use mpmorph, please cite the following papers:
 * Aykol, M., Dwaraknath, S.S., Sun, W. and Persson, K.A., 2018. Thermodynamic limit for synthesis of metastable inorganic materials. Science advances, 4(4), p.eaaq0148.
 * Aykol, M. and Persson, K.A., 2018. Oxidation Protection with Amorphous Surface Oxides: Thermodynamic Insights from Ab Initio Simulations on Aluminum. ACS applied materials & interfaces, 10(3), pp.3039-3045.
