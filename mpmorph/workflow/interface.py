from atomate.vasp.workflows.base.adsorption import MPSurfaceSet
from fireworks import Workflow
import uuid
from mpmorph.fireworks.core import MDFW, OptimizeFW, StaticFW
from custodian.vasp.handlers import VaspErrorHandler, MeshSymmetryErrorHandler, UnconvergedErrorHandler, \
    PotimErrorHandler, FrozenJobErrorHandler, NonConvergingErrorHandler, PositiveEnergyErrorHandler, \
    StdErrHandler

handler_group = [VaspErrorHandler(), MeshSymmetryErrorHandler(), UnconvergedErrorHandler(),
                 NonConvergingErrorHandler(), PotimErrorHandler(),
                 PositiveEnergyErrorHandler(), FrozenJobErrorHandler(), StdErrHandler()]


def get_interface_wf(interfaces, match, slab_substrates=None, slab_films=None,
                     bulk_substrate=None, bulk_film=None, strained_substrate=None, strained_film=None,
                     md_prerelax=True, skip_bulk=False, substrate_layers=None, film_layers=None,
                     interface_labels=None, md_temp=500, h_label='bare_surface', **kwargs):
    """
    This workflow currently does not allow for calculation of interfacial energy.

    Notes:
        1) All calculations with vacuum are preceded with a static calculation to pre-converge the orbitals.
           Due to this, all jobs must be run on the same machine, so fireworks can copy the chgcar over.
        2) Substrate surface energy as calculated with the outputs of this workflow do not include optimization.
        3) MD is used to pre-relax structures before input into the ISIF=2 relaxation
    """
    tag_id = uuid.uuid4()
    sub_label = f'{substrate_layers}_{bulk_substrate.composition.reduced_formula}({"".join([str(i) for i in match["sub_miller"]])})'
    film_label = f'{film_layers}_{bulk_film.composition.reduced_formula}({"".join([str(i) for i in match["film_miller"]])})'
    label = f'{film_label}/{sub_label}'

    fw_list = []

    # Bulk fireworks
    if not skip_bulk:
        if bulk_substrate is not None:
            sub_static_fw = StaticFW(structure=bulk_substrate, name=f'bulk_substrate_static-{str(tag_id)}',
                                     vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<')
            fw_list.append(sub_static_fw)

        if bulk_film is not None:
            film_static_fw = StaticFW(structure=bulk_film, name=f'bulk_film_static-{str(tag_id)}',
                                      vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<')
            fw_list.append(film_static_fw)

    # Strained bulk substrate calculation
    if strained_substrate is not None:
        # Optimize ion positions only
        strained_sub_opt_fw = OptimizeFW(structure=strained_substrate,
                                         name=f'{label}-strained_bulk_substrate_optimize-{str(tag_id)}',
                                         vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<',
                                         override_default_vasp_params={'user_incar_settings': {'ISIF': 2, 'KPAR': 4}},
                                         handler_group=handler_group, job_type='normal')
        strained_sub_static_fw = StaticFW(structure=strained_substrate,
                                          name=f'{label}-strained_bulk_substrate_static-{str(tag_id)}',
                                          vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<', previous_structure=True,
                                          parents=[strained_sub_opt_fw])
        fw_list.extend([strained_sub_opt_fw, strained_sub_static_fw])

    # Strained bulk film calculation
    if strained_film is not None:
        strained_film_opt_fw = OptimizeFW(structure=strained_film,
                                          name=f'{label}-strained_bulk_film_optimize-{str(tag_id)}',
                                          vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<',
                                          override_default_vasp_params={'user_incar_settings': {'ISIF': 2, 'KPAR': 4}},
                                          handler_group=handler_group, job_type='normal')
        strained_film_static_fw = StaticFW(structure=strained_film,
                                           name=f'{label}-strained_bulk_film_static-{str(tag_id)}',
                                           vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<', previous_structure=True,
                                           parents=[strained_film_opt_fw])
        fw_list.extend([strained_film_opt_fw, strained_film_static_fw])

    # Substrate slab calculation
    if slab_substrates is not None:
        for i, slab_substrate in enumerate(slab_substrates):
            run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": None, "spec": {}},
                        "optional_fw_params": {"override_default_vasp_params": {}}}
            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {"IDIPOL": 3, 'KPAR': 4, 'AMIN': 0.01, 'LWAVE': True}})
            sub_slab_orbital_fw = StaticFW(structure=slab_substrate,
                                           name=f'{sub_label}_{i}-substrate_slab_orbital-{str(tag_id)}',
                                           **run_args["run_specs"], **run_args["optional_fw_params"],
                                           previous_structure=False)

            run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                      "spec": {}},
                        "optional_fw_params": {"override_default_vasp_params": {}}}
            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {'ISTART': 1, 'LVTOT': True, "LDIPOL": True, "IDIPOL": 3, 'KPAR': 4,
                                         'AMIN': 0.01}})
            sub_slab_static_fw = StaticFW(structure=slab_substrate,
                                          name=f'{sub_label}_{i}-substrate_slab_static-{str(tag_id)}',
                                          **run_args["run_specs"], **run_args["optional_fw_params"],
                                          previous_structure=False,
                                          prev_calc_loc=True, parents=[sub_slab_orbital_fw])

            fw_list.extend([sub_slab_orbital_fw, sub_slab_static_fw])

    # Film slab calculation
    if slab_films is not None:
        for i, slab_film in enumerate(slab_films):
            if md_prerelax:
                run_args = {"md_params": {"start_temp": md_temp, "end_temp": md_temp, "nsteps": 200},
                            "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"},
                            "optional_fw_params": {"override_default_vasp_params": {}, "spec": {}},
                            "label": "md_prerelax_"}

                run_args["optional_fw_params"]["override_default_vasp_params"].update(
                    {'user_incar_settings': {'ISIF': 1, 'LWAVE': False, 'AMIN': 0.01, 'ALGO': "Normal"}})
                fw = MDFW(structure=slab_film,
                          name=f'{sub_label}_{i}-{h_label}_film_slab_md_prerelax_0-{str(tag_id)}',
                          previous_structure=False, insert_db=True, **run_args["md_params"],
                          **run_args["run_specs"], **run_args["optional_fw_params"], parents=[])
                fw_list.append(fw)

                for q in range(1, 5):
                    fw = MDFW(structure=slab_film,
                              name=f'{sub_label}_{i}-{h_label}_film_slab_md_prerelax_{q}-{str(tag_id)}',
                              previous_structure=True, insert_db=True, **run_args["md_params"],
                              **run_args["run_specs"], **run_args["optional_fw_params"], parents=[fw_list[-1]])
                    fw_list.append(fw)

                run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": None,
                                          "spec": {}},
                            "optional_fw_params": {"override_default_vasp_params": {}}}
                run_args["optional_fw_params"]["override_default_vasp_params"].update(
                    {'user_incar_settings': {"IDIPOL": 3, 'KPAR': 4, 'AMIN': 0.01, 'LWAVE': True}})
                film_slab_orbital_fw = StaticFW(structure=slab_film,
                                                name=f'{sub_label}_{i}-{h_label}_film_slab_orbital-{str(tag_id)}',
                                                **run_args["run_specs"], **run_args["optional_fw_params"],
                                                previous_structure=True, parents=[fw_list[-1]])
            else:
                run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": None,
                                          "spec": {}},
                            "optional_fw_params": {"override_default_vasp_params": {}}}
                run_args["optional_fw_params"]["override_default_vasp_params"].update(
                    {'user_incar_settings': {"IDIPOL": 3, 'KPAR': 4, 'AMIN': 0.01, 'LWAVE': True}})
                film_slab_orbital_fw = StaticFW(structure=slab_film,
                                                name=f'{sub_label}_{i}-{h_label}_film_slab_orbital-{str(tag_id)}',
                                                **run_args["run_specs"], **run_args["optional_fw_params"],
                                                previous_structure=True)

            run_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                      "spec": {}},
                        "optional_fw_params": {"override_default_vasp_params": {}}}
            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {'ISTART': 1, 'LVTOT': True, 'LDIPOL': True, 'IDIPOL': 3, 'AMIN': 0.01,
                                         'EDIFF': 0.000001 * slab_film.num_sites, 'KPAR': 4,
                                         'ICHARG': 0, 'POTIM': 0.5, 'NSW': 200, 'IBRION': 2}})
            vasp_input_set = MPSurfaceSet(slab_film, auto_dipole=True,
                                          **run_args['optional_fw_params']["override_default_vasp_params"])
            run_args["run_specs"]["vasp_input_set"] = vasp_input_set
            film_slab_opt_fw = OptimizeFW(structure=slab_film,
                                          name=f'{sub_label}_{i}-{h_label}_film_slab_optimize-{str(tag_id)}',
                                          previous_structure=True,
                                          parents=[film_slab_orbital_fw], **run_args["run_specs"],
                                          **run_args["optional_fw_params"],
                                          max_force_threshold=None, handler_group=handler_group,
                                          prev_calc_loc=True, job_type='normal')

            run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                      "spec": {}},
                        "optional_fw_params": {"override_default_vasp_params": {}}}
            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {'LVTOT': True, "LDIPOL": True, "IDIPOL": 3, 'KPAR': 4,
                                         'AMIN': 0.01}})

            film_slab_static_fw = StaticFW(structure=slab_film,
                                           name=f'{sub_label}_{i}-{h_label}_film_slab_static-{str(tag_id)}',
                                           parents=[film_slab_opt_fw],
                                           **run_args["run_specs"], **run_args["optional_fw_params"],
                                           previous_structure=True,
                                           prev_calc_loc=True)
            fw_list.extend([film_slab_orbital_fw, film_slab_opt_fw, film_slab_static_fw])

    # Interface calculations:
    for i, interface_structure in enumerate(interfaces):
        interface_label = interface_labels[i]

        if md_prerelax:
            run_args = {"md_params": {"start_temp": md_temp, "end_temp": md_temp, "nsteps": 200},
                        "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"},
                        "optional_fw_params": {"override_default_vasp_params": {}, "spec": {}},
                        "label": "md_prerelax_"}

            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {'ISIF': 1, 'LWAVE': False, 'AMIN': 0.01, 'ALGO': "Normal"}})
            fw = MDFW(structure=interface_structure,
                      name=f'{sub_label}_{interface_label}-{h_label}_interface_md_prerelax_0-{str(tag_id)}',
                      previous_structure=False, insert_db=True, **run_args["md_params"],
                      **run_args["run_specs"], **run_args["optional_fw_params"], parents=[])
            fw_list.append(fw)

            for q in range(1, 5):
                fw = MDFW(structure=interface_structure,
                          name=f'{sub_label}_{interface_label}-{h_label}_interface_md_prerelax_{q}-{str(tag_id)}',
                          previous_structure=True, insert_db=True, **run_args["md_params"],
                          **run_args["run_specs"], **run_args["optional_fw_params"], parents=[fw_list[-1]])
                fw_list.append(fw)

            run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": None,
                                      "spec": {}},
                        "optional_fw_params": {"override_default_vasp_params": {}}}
            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {"IDIPOL": 3, 'KPAR': 4, 'AMIN': 0.01, 'LWAVE': True}})
            interface_orbital_fw = StaticFW(structure=interface_structure,
                                            name=f'{label}_{interface_label}-{h_label}_interface_orbital-{str(tag_id)}',
                                            **run_args["run_specs"], **run_args["optional_fw_params"],
                                            previous_structure=True, parents=[fw_list[-1]])
        else:

            run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": None,
                                      "spec": {}},
                        "optional_fw_params": {"override_default_vasp_params": {}}}
            run_args["optional_fw_params"]["override_default_vasp_params"].update(
                {'user_incar_settings': {"IDIPOL": 3, 'KPAR': 4, 'AMIN': 0.01, 'LWAVE': True}})
            interface_orbital_fw = StaticFW(structure=interface_structure,
                                            name=f'{label}_{interface_label}-{h_label}_interface_orbital-{str(tag_id)}',
                                            **run_args["run_specs"], **run_args["optional_fw_params"],
                                            previous_structure=True)

        run_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                  "spec": {}},
                    "optional_fw_params": {"override_default_vasp_params": {}}}
        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISTART': 1, 'LVTOT': True, 'LDIPOL': True, 'IDIPOL': 3, 'AMIN': 0.01,
                                     'EDIFF': 0.000001 * interface_structure.num_sites,
                                     'EDIFFG': 0.00001 * interface_structure.num_sites, 'KPAR': 4,
                                     'ICHARG': 0, 'POTIM': 0.5, 'NSW': 200, 'IBRION': 2}})
        vasp_input_set = MPSurfaceSet(interface_structure, auto_dipole=True,
                                      **run_args['optional_fw_params']["override_default_vasp_params"])
        run_args["run_specs"]["vasp_input_set"] = vasp_input_set
        interface_opt_fw = OptimizeFW(structure=interface_structure,
                                      name=f'{label}_{interface_label}-{h_label}_interface_optimize-{str(tag_id)}',
                                      previous_structure=True,
                                      parents=[interface_orbital_fw], **run_args["run_specs"],
                                      **run_args["optional_fw_params"],
                                      max_force_threshold=None, handler_group=handler_group,
                                      prev_calc_loc=True, job_type='normal')

        run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                                  "spec": {}},
                    "optional_fw_params": {"override_default_vasp_params": {}}}
        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'LVTOT': True, "LDIPOL": True, "IDIPOL": 3, 'KPAR': 4,
                                     'AMIN': 0.01}})

        interface_static_fw = StaticFW(structure=interface_structure,
                                       name=f'{label}_{interface_label}-{h_label}_interface_static-{str(tag_id)}',
                                       parents=[interface_opt_fw],
                                       **run_args["run_specs"], **run_args["optional_fw_params"],
                                       previous_structure=True,
                                       prev_calc_loc=True)
        fw_list.extend([interface_orbital_fw, interface_opt_fw, interface_static_fw])

    return Workflow(fw_list, name=f'{label}-{h_label}_interface_wf')
