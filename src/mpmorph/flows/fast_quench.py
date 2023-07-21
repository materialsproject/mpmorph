from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker
from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core.structure import Structure

from jobflow.core.flow import Flow

relax_settings = {
            "LREAL": "Auto",
            "NSW": 1000,
            "LDAUPRINT": 0,
            "SIGMA": 0.05, # Per VASP docs for ISMEAR, which should be 0 for large cells
            "LAECHG": False,
            "EDIFFG": None, # Per old MPRelax
            "EDIFF": 0.005, # From Hui
            "ENCUT": 520,
            "NELM": 100,
            "ENAUG": None, # Per old MPRelax
            "LVTOT": False, # per old MPRelax
        }

def get_original_relax_flow(struct: Structure):
    relax_input_generator = RelaxSetGenerator(
        user_incar_settings=relax_settings,
        user_kpoints_settings=Kpoints.monkhorst_automatic(2, 2, 2)
    )
    
    relax_maker = RelaxMaker(input_set_generator=relax_input_generator, name="Original Structure Relaxation")
    relax_job = relax_maker.make(struct)

    static_input_generator = StaticSetGenerator()

    static_job = StaticMaker(
        input_set_generator=static_input_generator,
        name=f'Original Structure Static: {struct.composition.to_pretty_string()}',
    ).make(relax_job.output.structure)

    flow = Flow([relax_job, static_job], output=static_job.output, name=f'Orig. Structure Relax-Static: {struct.composition.to_pretty_string()}')
    return flow  

def get_fast_quench_flow(struct: Structure):
    relax_input_generator = RelaxSetGenerator(
        user_incar_settings=relax_settings,
        user_kpoints_settings=Kpoints.monkhorst_automatic(2, 2, 2)
    )
    
    relax_maker_1 = RelaxMaker(input_set_generator=relax_input_generator)
    relax_maker_2 = RelaxMaker(input_set_generator=relax_input_generator)    
    
    relax_job = DoubleRelaxMaker(
        relax_maker1=relax_maker_1,
        relax_maker2=relax_maker_2,
        name=f'Double Relax: {struct.composition.to_pretty_string()}'
    ).make(struct)
    
    static_input_generator = StaticSetGenerator()

    static_job = StaticMaker(
        input_set_generator=static_input_generator,
        name=f'Static: {struct.composition.to_pretty_string()}',
    ).make(relax_job.output.structure)

    flow = Flow([relax_job, static_job], output=static_job.output, name=f'Fast Quench: {struct.composition.to_pretty_string()}')
    return flow