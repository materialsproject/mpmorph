from pydantic import BaseModel, Field

class MDPVDataDoc(BaseModel):

    task_label: str = Field(None, description="The name of the task.")
    volume: float = Field(
        description = "The volume data from the MD run"
    )
    pressure: float = Field(
        description = "The volume data from the MD run"
    )
    