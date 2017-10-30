from fireworks import Workflow
from mpmorph.fireworks.core import MDFW, OptimizeFW, StaticFW
from mpmorph.fireworks import powerups
from pymatgen.core.surface import SlabGenerator
from mpmorph.fireworks.powerups import add_adsorbate_task
from mpmorph.util import recursive_update


def get_phonon_frequency_wf(structure, molecule, slab_gen_args, adsorbate_gen_args, structure_is_slab, priority=None,
                            converge=False, converge_args={}, converge_type=("density", 5), spawner_args={}, name = "phonon_frequency",
                            vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", spec={}):

    fw_list = []

    if not structure_is_slab:
        slab_args = {}
        recursive_update(slab_args, slab_gen_args)
        slab_gen = SlabGenerator(initial_structure=structure, **slab_args)
        structure = slab_gen.get_slab()

    if converge:
        #Converge slab
        run_args = {"md_params": {"start_temp": 3000, "end_temp": 3000, "nsteps": 2000},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                  "wall_time": 86400},
                    "optional_fw_params": {"override_default_vasp_params": {}, "copy_vasp_outputs": False, "spec": {}}}

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False}})
        run_args = recursive_update(run_args, converge_args)
        run_args["optional_fw_params"]["spec"]["_priority"] = priority

        fw = MDFW(structure=structure, name="run0", previous_structure=False, insert_db=True, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"])

        _spawner_args = {"converge_params": {"converge_type": [converge_type], "max_rescales": 10, "spawn_count": 0},
                         "rescale_args": {"beta": 0.0000004},
                         "run_specs": run_args["run_specs"], "md_params": run_args["md_params"],
                         "optional_fw_params": run_args["optional_fw_params"]}
        _spawner_args["md_params"].update({"start_temp": run_args["md_params"]["end_temp"]})
        _spawner_args = recursive_update(_spawner_args, spawner_args)

        fw = powerups.add_converge_task(fw, **_spawner_args)

        fw_list.append(fw)

    #Create fireworks to relax slab
    run_specs = {"vasp_cmd":vasp_cmd, "db_file":db_file, "spec":spec,
                 "optional_fw_params": {"override_default_vasp_params": {}, "copy_vasp_outputs": False, "spec": {}}}
    run_specs["optional_fw_params"]["spec"]["_priority"] = priority
    fw1 = OptimizeFW(structure=structure, name=structure.composition.reduced_formula + "_slab_optimize",
                     previous_structure=converge, **run_specs)
    fw2 = StaticFW(structure=structure, name=structure.composition.reduced_formula + "_slab_static",
                   previous_structure=True, parents=[fw1], **run_specs)


    fw2 = add_adsorbate_task(fw2, molecule=molecule, adsorbate_gen_args=adsorbate_gen_args, run_specs=run_specs)
    fw_list.extend([fw1, fw2])

    wf = Workflow(fw_list, name=name)
    return wf