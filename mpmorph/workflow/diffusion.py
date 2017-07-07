from fireworks import explicit_serialize, Firework, Workflow,

def get_diffusion(structure, temperatures, unconverged=True):
    wf = []
    if structure == None:
        #Run full amorphous workflow
        AmorphousFW


    #converge pressure on structure
    # "density", "total_energy", "kinetic_energy"

    if unconverged:
        ConvergeTask("density")

    #
    structure