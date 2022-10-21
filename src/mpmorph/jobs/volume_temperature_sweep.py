from jobflow import Flow, Maker, Response, job
import json
from .helpers import get_volume_at_temp_m3gnet_job

class VolumeTemperatureSweepMaker(Maker):

    name: str = "VOLUME_TEMPERATURE_SWEEP"

    @job
    def make(self, 
            structure,
            lower_bound = 100,
            upper_bound = 1100,
            temp_step = 100,
            output_name = 'vt.out', 
            steps = 2000):
        vs = []
        volume_jobs = []
        temps = list(range(lower_bound,upper_bound,temp_step))

        for temp in temps:
            job = get_volume_at_temp_m3gnet_job(structure, temp, steps = steps)
            volume_jobs.append(job)
            vs.append(job.output.volume)
        
        collect_job = _collect_vt_results(vs, temps, structure, output_name)

        new_flow = Flow([*volume_jobs, collect_job], output = collect_job.output)
        return Response(replace=new_flow)


@job
def _collect_vt_results(vs, ts, structure, output_fn):
    result = {
        "structure": structure.as_dict(),
        "volumes": vs,
        "temps": ts
    }
    print(result)
    with open(output_fn, '+w') as f:
        f.write(json.dumps(result))
    return result