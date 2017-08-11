from mpmorph.workflow.converge import get_converge

def get_diffusion(structure, temperatures=[500, 1000, 1500], prod_quants={"nsteps": 5000, "target": 40000}, converge_args={}, prod_args={}):
    wfs = []
    for temperature in temperatures:
        converge_args["md_params"]["start_temp"] = temperature
        converge_args["md_params"]["end_temp"] = temperature

        wf = get_converge(structure, prod_quants=prod_quants, converge_args=converge_args, prod_args=prod_args)
        wfs.append(wf)

    # Combine workflows
