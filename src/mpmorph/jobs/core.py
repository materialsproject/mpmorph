from dataclasses import dataclass, field
from pathlib import Path

from jobflow import Maker, job
from m3gnet.models import MolecularDynamics
from pymatgen.core import Structure

from mpmorph.jobs.schema import M3GNetCalculation


@dataclass
class M3GNetMaker(Maker):
    """
    Maker to create M3GNet MD calculation jobs.

    Args:
        name: The name of the job
        temperature: Temperature of the MD run

    """

    name: str = "m3gnet_run"
    temperature: float = 2000
    ensemble: str = "nvt"
    timestep: int = 1
    trajectory: str = "out.traj"
    logfile: str = "out.log"
    loginterval: int = 100

    @job(trajectory="trajectory", output_schema=M3GNetCalculation)
    def make(self, structure: Structure, steps=1000, **kwargs):
        """
        Run MD using M3GNet

        Args:
            structure: the input structure
            steps: number of steps to run
        """
        md = MolecularDynamics(
            atoms=structure,
            temperature=self.temperature,
            ensemble=self.ensemble,
            timestep=self.timestep,
            trajectory=self.trajectory,
            logfile=self.logfile,
            loginterval=self.loginterval,
            **kwargs
        )

        md.run(steps=steps)

        doc = M3GNetCalculation.from_directory(Path.cwd())
        doc.task_label = self.name

        return doc
