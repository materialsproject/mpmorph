from dataclasses import Field
from jobflow import Maker, job
from ..schemas.m3gnet_md_calc import M3GNetMDCalculation

from mpmorph.schemas.pv_data_doc import MDPVDataDoc

class ExtractPVDataM3gnet(Maker):
    """Simple job for extracting pressure-volume data from a VASP AIMD run.
    """

    name: str = "EXTRACT_PV_M3GNET"

    @job
    def make(self, m3gnet_calc: M3GNetMDCalculation):
        volume_data = m3gnet_calc_to_vol(m3gnet_calc)
        pressure_data = m3gnet_calc_to_pressure(m3gnet_calc)
        return MDPVDataDoc(
            volume = volume_data,
            pressure = pressure_data
        )

def m3gnet_calc_to_vol(m3gnet_calc: M3GNetMDCalculation):
    volume = m3gnet_calc.trajectory[-1].lattice.volume
    return volume

def m3gnet_calc_to_pressure(m3gnet_calc: M3GNetMDCalculation):
    # retrieve stresses in Voigt order (i.e. xx, yy, zz, yz, xz, xy)
    stresses = m3gnet_calc.trajectory.frame_properties[-1]["stress"]
    trace = sum(stresses[0:3])
    pressure = 1/3 * trace
    return pressure
