from pydantic import BaseModel, Field
from pymatgen.core.structure import Structure


class VTSweepDoc(BaseModel):
    task_label: str = Field(None, description="The name of the task.")
    volumes: float = Field(description="The volume at each temperature")
    temps: float = Field(
        description="The temperatures at which the volume was equilibrated"
    )
    structure: Structure = Field(
        description="The original structure for which this sweep was performed"
    )
