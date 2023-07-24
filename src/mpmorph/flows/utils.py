import json
import uuid
import pandas as pd

from jobflow import Flow, Maker, job

from mpmorph.jobs.equilibrate_volume import EquilibriumVolumeSearchMaker
from pymatgen.core.structure import Structure
from pymatgen.core.trajectory import Trajectory

from mpmorph.jobs.pv_from_calc import PVExtractor
from ..jobs.equilibrate_volume import PVFromMDFlowMaker, GetPVDocFromMDMaker


EQUILIBRATE_VOLUME_FLOW = "EQUILIBRATE_VOLUME_FLOW"
M3GNET_MD_FLOW = "M3GNET_MD_FLOW"
M3GNET_MD_CONVERGED_VOL_FLOW = "M3GNET_MD_CONVERGED_VOL_FLOW"
LAMMPS_VOL_FLOW = "LAMMPS_VOL_FLOW"
VOLUME_TEMPERATURE_SWEEP = "VOLUME_TEMPERATURE_SWEEP"

def get_frames_from_trajectory(trajectory: Trajectory, step_separation = 300, num_steps = 5):
    earliest_step = step_separation * num_steps
    steps = trajectory[-1:-earliest_step:-step_separation]
    return steps

def get_md_flow(
    pv_md_maker: Maker,
    pv_extractor: PVExtractor,
    production_md_maker: Maker,
    structure,
    converge_first,
    initial_vol_scale,
    flow_name: str = "MD_FLOW",
):
    struct = structure.copy()
    if initial_vol_scale is not None:
        struct.scale_lattice(struct.volume * initial_vol_scale)

    if converge_first:
        return get_converge_flow(
            structure = struct,
            pv_md_maker = pv_md_maker,
            pv_extractor = pv_extractor,
            production_run_maker = production_md_maker,
            flow_name=flow_name
        )
    else:
        return Flow([production_md_maker.make(struct)], name=flow_name)


def get_converge_flow(
    structure: Structure,
    pv_md_maker: Maker,
    pv_extractor: PVExtractor,
    production_run_maker: Maker,
    flow_name: str = M3GNET_MD_CONVERGED_VOL_FLOW
):
    eq_vol_maker = EquilibriumVolumeSearchMaker(
        pv_from_md_maker=PVFromMDFlowMaker(
            md_maker=pv_md_maker,
            extract_maker=GetPVDocFromMDMaker(
                pv_extractor=pv_extractor
            )
        )
    )

    equil_vol_job = eq_vol_maker.make(structure)

    final_md_job = production_run_maker.make(equil_vol_job.output)

    flow = Flow(
        [equil_vol_job, final_md_job],
        output=final_md_job.output,
        name=flow_name,
    )

    return flow

@job
def collect_vt_results(v_outputs, ts, structure, output_fn, mp_id):
    result = {
        "structure": structure.as_dict(),
        "volumes": [get_converged_vol(v) for v in v_outputs],
        "temps": ts,
        "mp_id": mp_id,
        "reduced_formula": structure.composition.reduced_formula,
        "formula": structure.composition.formula,
        "uuid": str(uuid.uuid4()),
    }

    with open(output_fn, "+w") as f:
        f.write(json.dumps(result))
    return result


def get_converged_vol(v_output):
    df = pd.DataFrame.from_dict(v_output)
    total_steps = (len(df) - 1) * 10
    avging_window = int(total_steps / 30)
    vols = df.iloc[-avging_window::]["vol"]
    eq_vol = vols.values.mean()
    return float(eq_vol)


