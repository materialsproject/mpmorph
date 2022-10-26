from pydantic import BaseModel, Field

from mpmorph.schemas.m3gnet_md_calc import M3GNetMDCalculation
from atomate2.vasp.schemas.task import TaskDocument

class MDPVDataDoc(BaseModel):

    task_label: str = Field(None, description="The name of the task.")
    volume: float = Field(
        description = "The volume data from the MD run"
    )
    pressure: float = Field(
        description = "The volume data from the MD run"
    )

class M3GnetPVDataDoc(MDPVDataDoc):

    m3gnet_calc: M3GNetMDCalculation = Field(
        description = "The m3gnet calculation from which this PV data was extracted"
    )

class VaspPVDataDoc(MDPVDataDoc):

    task_document: TaskDocument = Field(
        description = "The Vasp task document from which this PV data was extracted"
    )
    