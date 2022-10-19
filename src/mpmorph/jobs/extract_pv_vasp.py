from jobflow import Maker
import numpy as np

from atomate2.vasp.schemas.task import TaskDocument

from mpmorph.schemas.pv_data_doc import MDPVDataDoc

class ExtractPVDataFromVASPMDMaker(Maker):
    """Simple job for extracting pressure-volume data from a VASP AIMD run.
    """

    def make(self, task_doc: TaskDocument):
        volume_data = task_doc_to_volume(task_doc)
        pressure_data = task_doc_to_pressure(task_doc)
        return MDPVDataDoc(
            volume = volume_data,
            pressure = pressure_data
        )

def task_doc_to_volume(task_doc: TaskDocument) -> float:
    volume = task_doc.calcs_reversed[-1].output.ionic_steps[-1].structure.lattice.volume
    return volume

def task_doc_to_pressure(task_doc: TaskDocument) -> float: #TODO
    stress_tensor = task_doc.calcs_reversed[-1].output.ionic_steps[-1].stress
    pressure = 1/3 * np.trace(stress_tensor)
    return pressure
