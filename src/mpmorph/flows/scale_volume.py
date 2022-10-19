from jobflow import Flow, Maker
from pymatgen.core.structure import Structure

def md_to_volume_flow(structure: Structure,
                      scale_factor: float,
                      md_maker: Maker,
                      pv_maker: Maker): 
    struct = structure.copy()
    struct.scale_lattice(struct.volume * scale_factor)    
    md_job = md_maker.make(struct)
    pv_job = pv_maker.make(md_job.output)

    return Flow([
        md_job,
        pv_job
    ], output=pv_job.output)