from fireworks import Workflow, Firework
from atomate.vasp.fireworks.core import MDFW, OptimizeFW, StaticFW
from pymatgen.io.vasp.sets import MITMDSet
from mpmorph.runners.amorphous_maker import AmorphousMaker
from mpmorph.analysis.structural_analysis import get_sample_structures
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.write_inputs import ModifyIncar
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp import powerups
from pymatgen.core.structure import Structure
import os


def get_wf_density(structure, temperature, pressure_threshold=5.0, max_rescales=6, nsteps=2000, wall_time=19200,
                   vasp_input_set=None, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", name="density_finder",
                   optional_MDWF_params=None, override_default_vasp_params=None,
                   amorphous_maker_params=None, copy_calcs=False, calc_home="~/wflows",
                   cool=False, final_run = True, diffusion = True, priority_spec={}):
    """

    :param structure: Starting structure or dict of composition
    :param temperature: Temperature to run MD at
    :param pressure_threshold: Mean pressure for convergence
    :param max_rescales: Max iterations to converge pressure
    :param nsteps: number of steps for MD runs
    :param wall_time: Approximate time for vasp run
    :param vasp_input_set:
    :param vasp_cmd:
    :param db_file:
    :param name:
    :param optional_MDWF_params:
    :param override_default_vasp_params:
    :param amorphous_maker_params:
    :param copy_calcs:
    :param calc_home:
    :param cool:
    :param final_run:
    :param diffusion:
    :param priority_spec:
    :return:
    """
    # Ensure calculation storage folder exists and is not a duplicate
    if copy_calcs:
        if not os.path.exists(calc_home):
            raise ValueError("calc_home must be an existing folder.")
        elif os.path.exists(calc_home + "/" + name):
            raise ValueError("WF name already exists, choose a different name.")
        else:
            calc_home = os.path.join(calc_home, name)
            os.mkdir(calc_home)

    # If structure is in fact just a composition, create a random packed Structure!
    if not isinstance(structure, Structure) and isinstance(structure, dict):
        if not amorphous_maker_params:
            raise ValueError("amorphous_maker_params must be defined!")
        glass = AmorphousMaker(structure, **amorphous_maker_params)
        structure = glass.random_packed_structure

    # Create FW for first AIMD simulation
    fw1 = MDFW(structure=structure, start_temp=temperature, end_temp=temperature, nsteps=nsteps,
               name=name + "run0", vasp_input_set=vasp_input_set, db_file=db_file,
               vasp_cmd=vasp_cmd, wall_time=wall_time, override_default_vasp_params=override_default_vasp_params,
               **optional_MDWF_params)

    #Create FW that converges desired property
    t = []
    t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"]))
    if copy_calcs:
        t.append(CopyCalsHome(calc_home=calc_home, run_name="run0"))
    t.append(SpawnMDFWTask(pressure_threshold=pressure_threshold, max_rescales=max_rescales,
                           wall_time=wall_time, vasp_cmd=vasp_cmd, db_file=db_file,
                           copy_calcs=copy_calcs, calc_home=calc_home,
                           spawn_count=0, cool=cool, final_run=final_run, diffusion=diffusion,
                           temperature=temperature, priority_spec=priority_spec))

    fw2 = Firework(t, parents=[fw1], name=name + "_initial_spawn", spec=priority_spec)
    return Workflow([fw1, fw2], name=name + "_WF")


def get_wf_structure_sampler(xdatcar_file, n=10, steps_skip_first=1000, vasp_cmd=">>vasp_cmd<<",
                             db_file=">>db_file<<", name="structure_sampler", sim_anneal=False, copy_calcs=False,
                             calc_home="~/wflows", priority_spec={}, **kwargs):
    """
    :param xdatcar_file:
    :param n:
    :param steps_skip_first:
    :param vasp_cmd:
    :param db_file:
    :param name:
    :param sim_anneal:
    :param copy_calc:
    :param copy_home:
    :param kwargs:
    :return:
    """
    wfs = []
    structures = get_sample_structures(xdatcar_path=xdatcar_file, n=n, steps_skip_first=steps_skip_first)
    if sim_anneal:
        for i in range(len(structures)):
            s = structures[i]
            diffusion = True if i == 0 else False
            wflow_name=s.composition.reduced_formula
            _wf = get_simulated_anneal_wf(s, start_temp=2500, name='snap_' + str(i), diffusion=diffusion,
                                          wflow_name=wflow_name, calc_home=calc_home, copy_calcs=copy_calcs,
                                          db_file=db_file, snap_num=i, priority_spec=priority_spec)
            wfs.append(_wf)
    else:
        for s in structures:
            fw1 = OptimizeFW(s, vasp_cmd=vasp_cmd, db_file=db_file, parents=[], spec=priority_spec, **kwargs)
            fw2 = StaticFW(s, vasp_cmd=vasp_cmd, db_file=db_file, parents=[fw1], spec=priority_spec)
            wfs.append(Workflow([fw1, fw2], name=name + str(s.composition.reduced_formula)))
    return wfs


def get_relax_static_wf(structures, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<",
                        name="regular_relax", copy_calcs=False, calc_home="~/wflows", snap=0, priority_spec={}, **kwargs):
    """
    :param structures:
    :param vasp_cmd:
    :param db_file:
    :param name:
    :param kwargs:
    :return:
    """
    wfs = []
    calc_home = os.path.join(calc_home, "snap_" + str(snap))
    for s in structures:
        fw1 = OptimizeFW(s, vasp_cmd=vasp_cmd, db_file=db_file, parents=[], spec=priority_spec, **kwargs)
        fw2 = StaticFW(s, vasp_cmd=vasp_cmd, db_file=db_file, parents=[fw1], spec=priority_spec)
        t = [CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"])]
        if copy_calcs:
            t.append(CopyCalsHome(calc_home=calc_home, run_name=name))
        fw3 = Firework(t, name="relax_copy_calcs", parents=[fw2], spec=priority_spec)
        _wf = Workflow([fw1, fw2, fw3], name=str(s.composition.reduced_formula) + "_" + str(snap) + "_" + name)
        wfs.append(_wf)
    return wfs


def get_simulated_anneal_wf(structure, start_temp, snap_num, end_temp=500, temp_decrement=500, nsteps_cool=200, nsteps_hold=500,
                            wall_time=19200,
                            vasp_input_set=None, vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<", name="anneal",
                            optional_MDWF_params=None, override_default_vasp_params=None,
                            copy_calcs=False, calc_home="~/wflows", diffusion=False, wflow_name="", priority_spec={}):

    temperature = start_temp

    optional_MDWF_params = optional_MDWF_params or {}
    optional_MDWF_params['spec'] = priority_spec
    override_default_vasp_params = override_default_vasp_params or {}
    override_default_vasp_params['user_incar_settings'] = override_default_vasp_params.get('user_incar_settings') or {}
    override_default_vasp_params['user_incar_settings'].update({"ISIF": 1, "LWAVE": False})

    fw_list = []

    # Run first cool step
    fw1 = MDFW(structure=structure, start_temp=start_temp, end_temp=start_temp - temp_decrement, nsteps=nsteps_cool,
               name=name + "_cool_" + str(start_temp - temp_decrement), vasp_input_set=vasp_input_set, db_file=db_file,
               vasp_cmd=vasp_cmd, wall_time=wall_time, override_default_vasp_params=override_default_vasp_params,
               **optional_MDWF_params)
    fw_list.append(fw1)

    t = []
    t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"]))

    if copy_calcs:
        t.append(
            CopyCalsHome(calc_home=os.path.join(calc_home, name), run_name="cool_" + str(start_temp - temp_decrement)))

    # Run first hold step
    t.append(WriteSetTask(start_temp= start_temp - temp_decrement, end_temp = start_temp - temp_decrement, nsteps= nsteps_hold))
    t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                              handler_group="md", wall_time=wall_time, gzip_output=False))
    t.append(PassCalcLocs(name=str(name) + "_hold_" + str(start_temp - temp_decrement)))

    if copy_calcs:
        t.append(
            CopyCalsHome(calc_home=os.path.join(calc_home, name), run_name="hold_" + str(start_temp - temp_decrement)))

    fw_list.append(Firework(t, parents=[fw1], name=name + "_hold_" + str(start_temp - temp_decrement), spec=priority_spec))

    temperature -= temp_decrement

    while temperature > end_temp:
        # Cool Step
        cool_fw = get_mdfw(start_temp=temperature, end_temp=temperature-temp_decrement, nsteps=nsteps_cool, vasp_cmd=vasp_cmd, wall_time=wall_time,
                 name=name, copy_calcs=copy_calcs, calc_home=calc_home, parents=fw_list[len(fw_list)-1], priority_spec=priority_spec,
                           db_file=db_file, snap_num=snap_num)

        # Hold Step
        relax = temperature == end_temp+temp_decrement
        hold_fw = get_mdfw(start_temp=temperature-temp_decrement, end_temp=temperature-temp_decrement, nsteps=nsteps_hold,
                           vasp_cmd=vasp_cmd, wall_time=wall_time,name=name, copy_calcs=copy_calcs, calc_home=calc_home,
                           parents=cool_fw, priority_spec=priority_spec, relax=relax, diffusion=diffusion)


        if temperature == end_temp+temp_decrement:
            t.append(RelaxStaticTask(copy_calcs=copy_calcs, calc_home=calc_home, db_file=db_file, snap_num=snap_num, priority_spec=priority_spec))
            if diffusion:
                t.append(DiffusionTask(copy_calcs=copy_calcs, calc_home=calc_home, db_file=db_file, snap_num=snap_num, priority_spec=priority_spec))
        fw_list.append(Firework(t, name=name+"_hold_"+str(temperature-temp_decrement), parents=[fw_list[len(fw_list)-1]], spec=priority_spec))
        temperature -= temp_decrement

    wf = Workflow(fw_list, name=wflow_name + "_" + name + "simulated_anneal_WF")
    return wf


def get_mdfw(start_temp, end_temp, nsteps, vasp_cmd, wall_time, name, copy_calcs, calc_home, db_file, snap_num, parents=[],
             priority_spec={}, start_copy=False, relax=False, diffusion=False):
    t = []
    if start_copy:
        t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"]))
        if copy_calcs:
            t.append(
                CopyCalsHome(calc_home=os.path.join(calc_home, name), run_name="cool_" + str(end_temp)))

    t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True, additional_files=["XDATCAR", "OSZICAR", "DOSCAR"]))
    t.append(WriteSetTask(start_temp=start_temp, end_temp=end_temp, nsteps=nsteps))
    t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                              handler_group="md", wall_time=wall_time, gzip_output=False))
    t.append(PassCalcLocs(name=name + "_cool_" + str(end_temp)))
    if copy_calcs:
        t.append(
            CopyCalsHome(calc_home=os.path.join(calc_home, name), run_name="cool_" + str(end_temp)))

    if relax:
        t.append(RelaxStaticTask(copy_calcs=copy_calcs, calc_home=calc_home, db_file=db_file, snap_num=snap_num,
                                 priority_spec=priority_spec))
        if diffusion:
            t.append(DiffusionTask(copy_calcs=copy_calcs, calc_home=calc_home, db_file=db_file, snap_num=snap_num,
                                   priority_spec=priority_spec))
    return Firework(t, name=name + "_cool_" + str(end_temp), parents=parents, spec=priority_spec)


from mpmorph.workflow.old_temp_mdtasks import SpawnMDFWTask, CopyCalsHome, RelaxStaticTask, DiffusionTask, WriteSetTask
