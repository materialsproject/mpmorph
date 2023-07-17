from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker
from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core.structure import Structure

from jobflow.core.flow import Flow

def get_fast_quench_flow(struct: Structure):
    relax_input_generator = RelaxSetGenerator(
        user_kpoints_settings=Kpoints(),
        user_incar_settings={
            "LREAL": "Auto",
            "NSW": 500,
            "LDAUPRINT": 0
        }
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