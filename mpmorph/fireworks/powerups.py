from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
def add_converge_task(firework, **kwargs):
    spawner_task = SpawnMDFW(**kwargs)
    firework.tasks.append(spawner_task)
    return firework

def add_cont_structure(firework, vasp_input_set=None):
    structure = #Get structure from database
    t = WriteVaspFromIOSet(structure, vasp_input_set)
    firework.tasks.insert(1, t)
    return firework