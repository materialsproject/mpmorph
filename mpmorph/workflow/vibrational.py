from fireworks import Workflow
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from pymatgen.core.surface import SlabGenerator
from mpmorph.firetasks.vibtasks import AdsorbateGeneratorTask


def get_phonon_frequency_wf(structure, molecule, slab_gen_args, adsorbate_gen_args,
                            name = "phonon_frequency", vasp_cmd=">>vasp_cmd<<",
                            db_file=">>db_file<<", spec={}):

    #Generate Slab
    slab_gen = SlabGenerator(initial_structure=structure, **slab_gen_args)
    slab = slab_gen.get_slab()

    #Create fireworks to relax slab
    run_specs = {"vasp_cmd":vasp_cmd, "db_file":db_file, "spec":spec}

    fw1 = OptimizeFW(structure=slab, name=structure.composition.reduced_formula + "_slab_optimize",
                     **run_specs)
    fw2 = StaticFW(structure=slab, name=structure.composition.reduced_formula + "_slab_static",
                   prev_calc_loc=True, parents=[fw1], **run_specs)
    fw2.tasks.append(AdsorbateGeneratorTask(molecule=molecule, adsorbate_gen_args=adsorbate_gen_args, run_specs=run_specs))
    fws = [fw1, fw2]

    wf = Workflow(fws, name=name)
    return wf