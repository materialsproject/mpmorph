import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ase import units
from jobflow import Maker, job
from m3gnet.models import MolecularDynamics
from pymatgen.core import Structure

from ..schemas.m3gnet_md_calc import M3GNetMDCalculation


@dataclass
class M3GNetMDMaker(Maker):
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
        steps: The number of MD steps to run. Defaults to 1000.
        save_files: Whether to save the output and trajectory files. Defaults to False
            (i.e., files are removed from the directory where M3GNet was run).
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
    loginterval: int = 1
    append_trajectory: bool = False
    steps: int = 1000
    save_files: bool = False

    @job(trajectory="trajectory", output_schema=M3GNetMDCalculation)
    def make(self, structure: Structure, **kwargs):
        """
        Run MD using the M3GNet Molecular Dynamics interface. This runs molecular
        dynamics through the ASE interface with the default M3GNet potential.

        Args:
            structure: the input structure
        """

        if self.potential is not None:
            kwargs["potential"] = self.potential

        outfile_name = (
            f"{structure.composition.to_pretty_string()}-{self.temperature}K-{round(structure.volume, 3)}"
        )
        traj_fn = f"{outfile_name}.traj"
        log_fn = f"{outfile_name}.log"

        taut = 10 * units.fs

        md = MolecularDynamics(
            atoms=structure,
            ensemble=self.ensemble,
            temperature=self.temperature,
            timestep=self.timestep,
            pressure=self.pressure,
            taut=taut,
            taup=self.taup,
            compressibility_au=self.compressibility_au,
            trajectory=traj_fn,
            logfile=log_fn,
            loginterval=self.loginterval,
            append_trajectory=self.append_trajectory,
            **kwargs,
        )

        md.run(
            steps=self.steps,
        )

        d = M3GNetMDCalculation.from_directory(Path.cwd(), trajectory_fn=traj_fn)

        if not self.save_files:
            os.remove(traj_fn)
            os.remove(log_fn)

        d.task_label = self.name
        d.metadata = self.as_dict()

        return d
