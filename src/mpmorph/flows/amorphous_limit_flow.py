from pymatgen.core.trajectory import Trajectory

from jobflow import job, Response, Flow
from emmet.core.tasks import TaskDoc

from .vasp import get_md_flow_vasp
from .fast_quench import get_fast_quench_flow
from .utils import get_frames_from_trajectory
from ..runners.amorphous_maker import get_random_packed
from monty.json import MontyDecoder


def get_amorphous_limit_flow(comp, aimd_temp = 5000, aimd_steps_pv = 3000, aimd_steps_prod = 3000):
    original_amorphous_structure = get_random_packed(comp)

    md_flow = get_md_flow_vasp(
        original_amorphous_structure,
        temperature = aimd_temp,
        steps_prod=aimd_steps_prod,
        steps_pv=aimd_steps_pv,
    )


    fast_quenches =  make_quench_flow(md_flow.output)

    return Flow(
        [get_md_flow_vasp, fast_quenches],
        name=f'Amorphous Limit Flow: {comp}'
    )


@job(name="AMORPHOUS_LIMIT_MAKE_QUENCH_FLOWS")
def make_quench_flow(md_flow_output: TaskDoc):
    flow = get_quench_flow_from_aimd(md_flow_output)
    return Response(replace=flow)

def get_quench_flow_from_aimd(md_flow_output: TaskDoc):
    traj_dict = md_flow_output.vasp_objects["trajectory"]
    trajectory: Trajectory = MontyDecoder().process_decoded(traj_dict)

    frames = get_frames_from_trajectory(trajectory)
    flows = [get_fast_quench_flow(frame) for frame in frames]
    return Flow(flows, name=f'Fast Quench Frames: {trajectory[0].composition.to_pretty_string()}')



