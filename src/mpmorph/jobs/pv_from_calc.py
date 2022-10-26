from dataclasses import dataclass
from jobflow import Maker, job, Flow, Response


from atomate2.vasp.schemas.task import TaskDocument
from atomate2.vasp.jobs.core import MDMaker

from mpmorph.schemas.pv_data_doc import VaspPVDataDoc
from mpmorph.jobs.core import M3GNetMDMaker
from ..schemas.m3gnet_md_calc import M3GNetMDCalculation

from mpmorph.schemas.pv_data_doc import M3GnetPVDataDoc


@dataclass
class PVFromCalc(Maker):
    @job
    def make(self, structure, scale_factor=None):
        struct = structure.copy()
        if scale_factor is not None:
            struct.scale_lattice(struct.volume * scale_factor)

        calc_job = self.md_maker.make(struct)

        doc_job = self.build_doc(
            calc_job.output,
        )

        return Response(replace=Flow([calc_job, doc_job], output=doc_job.output))


@dataclass
class PVFromM3GNet(PVFromCalc):

    name: str = "PV_FROM_M3GNET"
    md_maker: Maker = M3GNetMDMaker()

    @job
    def build_doc(self, m3gnet_calc: M3GNetMDCalculation):
        v_data = m3gnet_calc_to_vol(m3gnet_calc)
        p_data = m3gnet_calc_to_pressure(m3gnet_calc)
        return M3GnetPVDataDoc(volume=v_data, pressure=p_data, m3gnet_calc=m3gnet_calc)


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
class PVFromVasp(PVFromCalc):

    name: str = "PV_FROM_VASP"
    md_maker: Maker = MDMaker()

    @job
    def build_doc(self, task_document: M3GNetMDCalculation):
        v_data = task_doc_to_volume(task_document)
        p_data = task_doc_to_pressure(task_document)
        return VaspPVDataDoc(
            volume=v_data, pressure=p_data, task_document=task_document
        )


@job
def task_doc_to_volume(task_doc: TaskDocument) -> float:
    volume = task_doc.calcs_reversed[-1].output.ionic_steps[-1].structure.lattice.volume
    return volume


@job
def task_doc_to_pressure(task_doc: TaskDocument) -> float:  # TODO
    stress_tensor = task_doc.calcs_reversed[-1].output.ionic_steps[-1].stress
    pressure = 1 / 3 * np.trace(stress_tensor)
    return pressure
