from pydantic import BaseModel, Field
from pymatgen.core.composition import Composition
from pymatgen.core.trajectory import Trajectory as PmgTrajectory

from mpmorph.utils import datetime_str


class LammpsCalc(BaseModel):
    task_label: str = Field(None, description="The name of the task.")
    dir_name: str = Field(
        None, description="The directory where the LAMMPS calculation was run"
    )
    last_updated: str = Field(
        default_factory=datetime_str,
        description="Timestamp of when the document was last updated.",
    )
    trajectory: PmgTrajectory = Field(
        None, description="The pymatgen Trajectory object stored ad dictionary"
    )
    composition: Composition = Field(description="The composition of the structure.")
    reduced_formula: str = Field(
        description="The reduced formula of the structure's composition."
    )
    dump_data: dict = Field(
        None, description="Any additional data collected via LAMMPS dump files"
    )
    metadata: dict = Field(
        None,
        description=(
            "Important info about the calculation, including ensemble type,"
            " temperature, etc."
        ),
    )
