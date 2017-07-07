from fireworks import explicit_serialize, Firework, Workflow


def get_diffusion(structure, temperatures, unconverged=True):
    wf = []
    if structure == None:
        # Run full amorphous workflow with only one snap
        get_amorphous()

    # TODO: Converge
    # TODO: Production Run
    # TODO: Sampling & Sim Anneal/MP Relax
    # TODO: Diffusion Converge
    # TODO: Production Diffusion Run
    # converge pressure on structure
    # "density", "total_energy", "kinetic_energy"

    if unconverged:
        ConvergeTask("density")

    #
    structure