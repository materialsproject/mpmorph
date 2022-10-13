from dataclasses import dataclass, field

from jobflow import SETTINGS, Maker, job

from mpmorph.jobs.schema import M3GNetCalculation


@dataclass
class M3GNetMaker(Maker):

    name: str = "m3gnet_run"
    trajectory_db_name: str = "trajectory_db"

    @job(trajectory="trajectory", output_schema=M3GNetCalculation)
    def make(self, structure):
        trajectory_db = SETTINGS.JOB_STORE.additional_stores.get(self.entry_db_name)
