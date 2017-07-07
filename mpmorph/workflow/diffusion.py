from fireworks import explicit_serialize, Firework, Workflow
from pymatgen import Structure, Composition
from mpmorph.fireworks.core import ConvergeFW, DiffusionFW

def get_diffusion(structure, temperatures, unconverged=True, copy_calcs=False, calc_home="./wflows"):
    fw_list = []
    #If structure is composition, get amorphous structure first
    if not isinstance(structure, Structure):
        #Run full amorphous workflow with only one snap
        get_amorphous()

    #TODO: Diffusion Converge
    #converge pressure on structure
    # "density", "total_energy", "kinetic_energy"
    if unconverged:
        fw_list.append(ConvergeFW("density"))

    #TODO: Production Diffusion Run
    fw_list.append(DiffusionFW(temps))

    pretty_name=structure.composition.reduced_formula
    wf = Workflow(fireworks=fw_list, name = pretty_name + "_diffusion")

