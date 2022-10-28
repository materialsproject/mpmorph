from typing import Optional

from pydantic import BaseModel, Field
from ase import units

def one_atmosphere():
    return 1.01325 * units.bar

class M3GNetMDInputs(BaseModel):

    ensemble: str = Field("nvt", description="The ensemble to use during the MD run.")

    temperature: float = Field(2000.0, description="The temperature of the simulation")

    pressure: float = Field(
        default_factory=one_atmosphere,
        description="The pressure at which this simulation should be run."
    )

    timestep: float = Field(
        1.0, 
        description="The timestep for this simulation in femtoseconds"
    )

    potential: Optional[str] = Field(
        None,
        description="The potential for use in this simulation"
    )
    taut: Optional[float] = Field(
        None, 
        description="The time constant for Berendsen temperature coupling"
    )
    taup: Optional[float] = Field(
        None,
        description="The time constant for pressure coupling"
    )
    compressibility_au: Optional[float] = Field(
        None,
        description="The compressibility of the material in A^3/eV"
    )
    trajectory_fn: str = Field(
        "out.traj", 
        description="The name of file which should be used to store the trajectory"
    )
    logfile_fn: str = Field(
        "out.log", 
        description="The name of the file for the simulation log"
    )
    loginterval: int = Field(
        1,
        description="Write to the logfile every {this number} steps"
    )
    append_trajectory: bool = Field(
        False, 
        description="Whether or not the current trajectory should be appended to the exitisting traj file"
    )
    steps: int = Field(
        1000,
        description="The number of steps for which this simulation should run"
    )
    save_files: bool = Field(
        True,
        description="Whether or not the files from this simulation should be saved"
    )





