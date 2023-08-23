from jobflow import Flow, Maker, Response, job
from mpmorph.jobs.tasks.m3gnet_input import M3GNetMDInputs
import dataclasses


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
            params = dataclasses.replace(self.md_parameters)
            params.temperature = temp
            params.steps = steps(structure, params)
            volume_jobs.append(job)
            vs.append(job.output.volume)

        collect_job = _collect_vt_results(vs, temps, structure, output_name)

        new_flow = Flow([*volume_jobs, collect_job], output=collect_job.output)
        return Response(replace=new_flow)
