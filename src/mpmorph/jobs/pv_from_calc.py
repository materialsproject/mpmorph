from dataclasses import dataclass

from emmet.core.tasks import TaskDoc

from ..schemas.m3gnet_md_calc import M3GNetMDCalculation
from ..schemas.chgnet_md_calc import CHGNetMDCalculation

import numpy as np

from abc import ABC, abstractmethod
from monty.json import MSONable

class PVExtractor(ABC, MSONable):

    @abstractmethod
    def get_pressure(self, md_output):
        pass

    @abstractmethod
    def get_volume(self, md_output):
        pass


class PVFromM3GNet(PVExtractor):

    def get_pressure(self, m3gnet_calc: M3GNetMDCalculation):
        return m3gnet_calc_to_pressure(m3gnet_calc)

    def get_volume(self, m3gnet_calc: M3GNetMDCalculation):
        return m3gnet_calc_to_vol(m3gnet_calc)

def m3gnet_calc_to_vol(m3gnet_calc: M3GNetMDCalculation):
    volume = m3gnet_calc.trajectory[-1].lattice.volume
    return volume

def m3gnet_calc_to_pressure(m3gnet_calc: M3GNetMDCalculation):
    # retrieve stresses in Voigt order (i.e. xx, yy, zz, yz, xz, xy)
    stresses = m3gnet_calc.trajectory.frame_properties[-1]["stress"]
    trace = sum(stresses[0:3])
    pressure = 1 / 3 * trace
    return pressure

@dataclass
class PVFromCHGNet(PVExtractor):

    def get_pressure(self, chgnet_calc: CHGNetMDCalculation):
        return chgnet_calc_to_pressure(chgnet_calc)

    def get_volume(self, chgnet_calc: CHGNetMDCalculation):
        return chgnet_calc_to_vol(chgnet_calc)


def chgnet_calc_to_vol(chgnet_calc: CHGNetMDCalculation):
    volume = chgnet_calc.trajectory[-1].lattice.volume
    return volume

def chgnet_calc_to_pressure(chgnet_calc: CHGNetMDCalculation):
    # retrieve stresses in Voigt order (i.e. xx, yy, zz, yz, xz, xy)
    stresses = chgnet_calc.trajectory.frame_properties[-1]["stress"]
    trace = sum(stresses[0:3])
    pressure = 1 / 3 * trace
    return pressure
    
class PVFromVasp(PVExtractor):

    def get_volume(self, task_document: TaskDoc):
        return task_doc_to_volume(task_document)

    def get_pressure(self, task_document: TaskDoc):
        return task_doc_to_pressure(task_document)


def task_doc_to_volume(task_doc: TaskDoc) -> float:
    volume = task_doc.calcs_reversed[-1].output.ionic_steps[-1].structure.lattice.volume
    return volume


def task_doc_to_pressure(task_doc: TaskDoc) -> float:  # TODO
    stress_tensor = task_doc.calcs_reversed[-1].output.ionic_steps[-1].stress
    pressure = 1 / 3 * np.trace(stress_tensor)
    return pressure
