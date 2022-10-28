from jobflow import Flow, Maker, Response, job
import json

from mpmorph.jobs.tasks.m3gnet_input import M3GNetMDInputs
from .helpers import get_volume_at_temp_m3gnet_job


class VolumeTemperatureSweepMaker(Maker):

    name: str = "VOLUME_TEMPERATURE_SWEEP"
    md_parameters: M3GNetMDInputs = None
    
    @job
    def make(
        self,
        structure,
        lower_bound=100,
        upper_bound=1100,
        temp_step=100,
        output_name="vt.out",
        steps=2000,
    ):
        if self.md_parameters is None:
            self.md_parameters = M3GNetMDInputs()

        vs = []
        volume_jobs = []
        temps = list(range(lower_bound, upper_bound, temp_step))

        

        for temp in temps:
            params = self.md_parameters.copy()
            params.temperature = temp
            params.steps = steps
            job = get_volume_at_temp_m3gnet_job(structure, params)
            volume_jobs.append(job)
            vs.append(job.output.volume)

        collect_job = _collect_vt_results(vs, temps, structure, output_name)

        new_flow = Flow([*volume_jobs, collect_job], output=collect_job.output)
        return Response(replace=new_flow)


@job
def _collect_vt_results(vs, ts, structure, output_fn):
    result = {"structure": structure.as_dict(), "volumes": vs, "temps": ts}

    with open(output_fn, "+w") as f:
        f.write(json.dumps(result))
    return result
