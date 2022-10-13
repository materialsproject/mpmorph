from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ase import units
from jobflow import Maker, job
from m3gnet.models import MolecularDynamics
from pymatgen.core import Structure

from mpmorph.jobs.schema import M3GNetMDCalculation


@dataclass
class M3GNetMaker(Maker):
    """
    Maker to create M3GNet MD calculation jobs.

    Args:
        name: The name of the job. Defaults to "m3gnet_run".
        ensemble: The ensemble type to use for the MD calculation
            (options: ["nvt", "npt", "npt_berendsen"]). Defaults to "nvt".
            Use npt at your own risk.
        temperature (Kelvin): The absolute temperature to use for the MD calculation
        pressure (eV/A^3): The modeled pressure in the NPT ensemble. Defaults to 1 atm,
            which is 1.01325 bar = 6.42151e-07 eV/A^3.
        timestep (fs): The time step to use for the MD calculation, in femtoseconds.
            Defaults to 1 fs.
        potential: The potential to use for the MD calculation. Uses default model from
            m3gnet.models.MolecularDynamics if None.
        taut: Optional time constant for Berendsen temperature coupling
        taup: Optional time constant for pressure coupling
        compressibility_au: Optional compressibility of the material in A^3/eV
        trajectory_fn The name of the trajectory file to write. Defaults to "out.traj".
        logfile_fn: The name of the logfile to write. Defaults to "out.log".
        loginterval: The interval number of steps at which to write to the logfile.
            Defaults to 10.
        append_trajectory: Whether to append to previous trajectory file if it already
            exists. Defaults to False.
    """

    name: str = "m3gnet_run"
    ensemble: str = "nvt"
    temperature: float = 2000.0
    pressure: float = 1.01325 * units.bar
    timestep: float = 1.0
    potential: Optional[str] = None
    taut: Optional[float] = None
    taup: Optional[float] = None
    compressibility_au: Optional[float] = None
    trajectory_fn: str = "out.traj"
    logfile_fn: str = "out.log"
    loginterval: int = 10
    append_trajectory: bool = False

    @job(trajectory="trajectory", output_schema=M3GNetMDCalculation)
    def make(self, structure: Structure, steps: int = 1000, **kwargs):
        """
        Run MD using the M3GNet Molecular Dynamics interface. This runs molecular
        dynamics through the ASE interface with the default M3GNet potential.

        Args:
            structure: the input structure
            steps: number of steps to run
        """

        if self.potential is not None:
            kwargs["potential"] = self.potential

        md = MolecularDynamics(
            atoms=structure,
            ensemble=self.ensemble,
            temperature=self.temperature,
            timestep=self.timestep,
            pressure=self.pressure,
            taut=self.taut,
            taup=self.taup,
            compressibility_au=self.compressibility_au,
            trajectory=self.trajectory_fn,
            logfile=self.logfile_fn,
            loginterval=self.loginterval,
            append_trajectory=self.append_trajectory,
            **kwargs
        )

        md.run(steps=steps)

        d = M3GNetMDCalculation.from_directory(Path.cwd())
        d.task_label = self.name
        d.metadata = self.as_dict()

        return d
