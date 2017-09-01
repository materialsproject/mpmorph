# Create Adsorbtion Structures
from fireworks import FireTaskBase, Firework, Workflow, FWAction, explicit_serialize
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.io.vasp import Poscar
from pymatgen import Structure

from pymatgen.analysis.adsorption import AdsorbateSiteFinder


@explicit_serialize
class AdsorbateGeneratorTask(FireTaskBase):
    required_params = ["run_specs", "molecule"]
    optional_params = ["adsorbate_gen_args"]

    def run_task(self, fw_spec):
        run_specs = self["run_specs"]
        molecule = self["molecule"]
        adsorbate_gen_args = self.get("adsorbate_gen_args", {"find_args": {"distance": 2.0}, "repeat": [1, 1, 1]})

        # Load structure from file
        _poscar = Poscar.from_file("CONTCAR.gz")
        structure = _poscar.structure

        adsorb_finder = AdsorbateSiteFinder(slab=structure, selective_dynamics=True)
        ad_structs = adsorb_finder.generate_adsorption_structures(molecule=molecule, **adsorbate_gen_args)

        opt_fws = []
        stat_fws = []
        for ad_struct in ad_structs:
            _struct_dict = ad_struct.as_dict()
            for site in __struct_dict["sites"]:
                if "velocities" in site["properties"].keys():
                    del site["properties"]["velocities"]

            _ad_struct = Structure.from_dict(_struct_dict)
            opt_fws.append(
                OptimizeFW(structure=_ad_struct, name=structure.composition.reduced_formula + "_slab_optimize",
                           **run_specs))
            stat_fw = StaticFW(structure=_ad_struct, name=structure.composition.reduced_formula + "_slab_static",
                               prev_calc_loc=True, **run_specs, parents=[opt_fws[-1]])
            stat_fw.tasks.append(SpawnVibrationalFWTask(run_specs=run_specs))
            stat_fws.append(stat_fw)
        fws = opt_fws
        fws.extend(stat_fws)
        wf = Workflow(fws)

        return FWAction(additions=wf)


@explicit_serialize
class SpawnVibrationalFWTask(FireTaskBase):
    required_params = ["run_specs"]
    optional_params = ["incar_updates"]

    def run_task(self, fw_spec):
        run_specs = self["run_specs"]
        incar_updates = self.get("incar_updates", {})

        # Load structure from file
        _poscar = Poscar("CONTCAR")
        structure = _poscar.structure
        fw = VibrationalFW(structure, incar_updates=incar_updates)
        wf = Workflow(fw)
        return FWAction(additions=wf)

from mpmorph.fireworks.core import VibrationalFW