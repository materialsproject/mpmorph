from jobflow import Flow, Maker

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure

from mpmorph.jobs.pv_from_calc import PVFromCalc


EQUILIBRATE_VOLUME_FLOW = "EQUILIBRATE_VOLUME_FLOW"
M3GNET_MD_FLOW = "M3GNET_MD_FLOW"
M3GNET_MD_CONVERGED_VOL_FLOW = "M3GNET_MD_CONVERGED_VOL_FLOW"
LAMMPS_VOL_FLOW = "LAMMPS_VOL_FLOW"


def get_md_flow(
    pv_md_maker,
    production_md_maker,
    structure,
    converge_first,
    initial_vol_scale,
    flow_name = "MD Flow"
):
    struct = structure.copy()
    if initial_vol_scale is not None:
        struct.scale_lattice(struct.volume * initial_vol_scale)

    if converge_first:
        return get_converge_flow(struct, pv_md_maker, production_md_maker, flow_name=flow_name)
    else:
        return Flow([production_md_maker.make(struct)], name=flow_name)


def get_converge_flow(
    structure: Structure,
    pv_md_maker: PVFromCalc,
    production_run_maker: Maker,
    flow_name: str
):
    eq_vol_maker = EquilibriumVolumeSearchMaker(pv_md_maker=pv_md_maker)

    equil_vol_job = eq_vol_maker.make(structure)

    final_md_job = production_run_maker.make(equil_vol_job.output)

    flow = Flow(
        [equil_vol_job, final_md_job],
        output=final_md_job.output,
        name=flow_name,
    )

    return flow


