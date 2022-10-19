from jobflow import Flow

def md_to_volume_flow(structure, scale_factor, md_maker, pv_maker): 
    struct = structure.copy()
    struct.scale_lattice(struct.volume * scale_factor)    
    md_job = md_maker.make(struct)
    pv_job = pv_maker.make(md_job.output)

    return Flow([
        md_job,
        pv_job
    ], output=pv_job.output)