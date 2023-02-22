from typing import Optional

from dataclasses import dataclass

from ase import units


def one_atmosphere():
    return 1.01325 * units.bar


@dataclass
class M3GNetMDInputs:
    ensemble: str = "nvt"
    temperature: float = 2000.0
    pressure: float = 1.01325 * units.bar
    timestep: float = 1.0
    potential: str = None
    taut: Optional[float] = None
    taup: Optional[float] = None
    compressibility_au: Optional[float] = None
    trajectory_fn: str = "out.traj"
    logfile_fn: str = "out.log"
    loginterval: int = 1
    append_trajectory: bool = False
    steps: int = 1000
    save_files: bool = True
