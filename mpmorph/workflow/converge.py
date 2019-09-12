import uuid
from copy import deepcopy

from fireworks import Workflow
from mpmorph.fireworks import powerups
from mpmorph.fireworks.core import MDFW
from mpmorph.util import recursive_update

__author__ = 'Eric Sivonxay, Jianli Cheng, and Muratahan Aykol'
__maintainer__ = 'Eric Sivonxay'
__email__ = 'esivonxay@lbl.gov'


def get_converge_wf(structure, temperature, converge_scheme='EOS', priority=None,
                    max_steps=5000, target_steps=10000, preconverged=False,
                    notes=None, save_data="all", aggregate_trajectory=True, **kwargs):
    """
    :param structure: Starting structure for the run
    :param temperature: Temperature for which to obtain a liquid
    :param images: Perturbations to the volume for the fit (fraction of input structure volume)
    :param preconverged: Is the structure already converged (i.e. Pressure 0bar) or volume rescaling not desired?
    :param max_steps: Maximum number of steps per chunk of production run MD simulation
    :param target_steps: Target number of steps for production MD run
    :param priority: Priority of all fireworks in the workflow
    :param kwargs: arguments such as spawner_args, converge_args, convergence_criteria, tag_id, prod_count, etc.
    :param notes: Any additional comments to propagate with this run
    :param save_data: Level to save job outputs. Options are "all", 'production', and None

    :return:
    """
    # Generate a unique identifier for the fireworks belonging to this workflow
    tag_id = kwargs.get('tag_id', uuid.uuid4())
    prod_count = kwargs.get('prod_count', 0)

    # To aggregate trajectory, the job output from the production runs must be saved.
    if aggregate_trajectory and save_data == None:
        save_data = "production"

    fw_list = []

    # Setup initial Run and convergence of structure
    run_args = {"md_params": {"start_temp": temperature, "end_temp": temperature, "nsteps": 2000},
                "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"},
                "optional_fw_params": {
                    "override_default_vasp_params": {'user_incar_settings': {'ISIF': 1, 'LWAVE': False,
                                                                             'PREC': 'Normal'}},
                    "spec": {'_priority': priority}
                }
                }
    run_args = recursive_update(run_args, kwargs.get('converge_args', {}))

    # Setup Dictionary specifying parameters of the spawner for convergence tasks
    _spawner_args = {
        "converge_params": {"max_rescales": 15, "density_spawn_count": 1, "energy_spawn_count": 0,
                            'converge_type': kwargs.get('convergence_criteria',
                                                        [("density", 5), ('ionic', 0.001)])},
        "rescale_params": {"beta": 5e-7},
        "run_specs": run_args["run_specs"],
        "md_params": run_args["md_params"],
        "optional_fw_params": run_args["optional_fw_params"],
        "tag_id": tag_id
    }
    _spawner_args["md_params"].update({"start_temp": run_args["md_params"]["end_temp"]})
    _spawner_args = recursive_update(_spawner_args, kwargs.get('spawner_args', {}))

    # Converge the pressure (volume) of the system
    if not preconverged:
        insert_converge_data = True if save_data == "all" else False
        if converge_scheme == 'EOS':
            # Create structures for varying volumes
            images = kwargs.get('image_scale', [0.8, 1, 1.2])
            structures = [structure.copy() for i in images]
            for i, factor in enumerate(images):
                structures[i].scale_lattice(structure.volume * factor)

            # Create firework for each structure
            EOS_run_args = deepcopy(run_args)
            EOS_run_args['md_params']['nsteps'] = 1000  # Use less steps... Pressure usually converges rapidly
            EOS_run_args = recursive_update(EOS_run_args, kwargs.get('converge_args', {}))
            volume_fws = []
            for n, (i, vol_structure) in enumerate(zip(images, structures)):
                save_structure = True if n == len(images) - 1 else False
                _fw = MDFW(structure=vol_structure, name=f'volume_{str(i)}-{str(tag_id)}',
                           previous_structure=False, insert_db=insert_converge_data, save_structure=save_structure,
                           **EOS_run_args["md_params"], **EOS_run_args["run_specs"],
                           **EOS_run_args["optional_fw_params"])

                _fw = powerups.add_pass_pv(_fw)
                volume_fws.append(_fw)
            fw_list.extend(volume_fws)

            # Create firework to converge pressure/volume
            spawner_fw = MDFW(structure=structure, name="run1" + "-" + str(tag_id),
                              previous_structure=True, insert_db=insert_converge_data, parents=volume_fws,
                              **run_args["md_params"], **run_args["run_specs"], **run_args["optional_fw_params"])

            spawner_fw = powerups.add_pv_volume_rescale(spawner_fw)
            spawner_fw = powerups.add_pass_pv(spawner_fw)
            spawner_fw = powerups.add_converge_task(spawner_fw, **_spawner_args)
            fw_list.append(spawner_fw)
        else:
            fw1 = MDFW(structure=structure, name="run0" + "-" + str(tag_id),
                       previous_structure=False, insert_db=insert_converge_data,
                       **run_args["md_params"], **run_args["run_specs"], **run_args["optional_fw_params"])
            fw1 = powerups.add_converge_task(fw1, **_spawner_args)
            fw_list.append(fw1)

    # Production length MD runs
    insert_prod_data = True if save_data == "all" or save_data == "production" else False
    prod_steps = 0
    while prod_steps <= target_steps - max_steps:
        # Create Dictionary with production run parameters
        run_args = {"md_params": {"start_temp": run_args["md_params"]["end_temp"],
                                  "end_temp": run_args["md_params"]["end_temp"],
                                  "nsteps": max_steps},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"},
                    "optional_fw_params": {"override_default_vasp_params":
                                               {'user_incar_settings': {'ISIF': 1, 'LWAVE': False,
                                                                        'PREC': 'Normal'}},
                                           "spec": {'_priority': priority}}}
        run_args = recursive_update(run_args, kwargs.get('prod_args', {}))

        parents = fw_list[-1] if len(fw_list) > 0 else []
        previous_structure = False if preconverged and prod_steps == 0 else True
        fw = MDFW(structure=structure, name=f'{str(temperature)}_prod_run_{str(prod_count)}-{str(tag_id)}',
                  previous_structure=previous_structure, insert_db=insert_prod_data, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"], parents=parents)
        fw_list.append(fw)

        prod_steps += max_steps
        prod_count += 1

    pretty_name = structure.composition.reduced_formula

    if aggregate_trajectory:
        fw_list[-1] = powerups.aggregate_trajectory(fw_list[-1], tag_id=tag_id, notes=notes,
                                                    db_file=run_args["run_specs"]["db_file"])

    wf = Workflow(fireworks=fw_list, name=f'{pretty_name}_{temperature}_diffusion')
    return wf
