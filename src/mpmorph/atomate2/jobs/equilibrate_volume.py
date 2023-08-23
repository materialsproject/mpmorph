from dataclasses import dataclass
from typing import List, Tuple
import numpy as np

from pymatgen.core.structure import Structure
from jobflow import Maker, job, Flow, Response

from .pv_from_calc import PVExtractor
from ..schemas.pv_data_doc import MDPVDataDoc
from ...runners import rescale_volume


MAX_MD_JOBS = 9  # if you can't converge with 6 additional calcs you're doing something wrong...
BUFFER = 0.1  # gives it enough room to slosh back

def get_scaled_structure(struct: Structure, scale_factor: float):
    copy: Structure = struct.copy()
    copy.scale_lattice(copy.volume * scale_factor)
    return copy

@dataclass
class GetPVDocFromMDMaker(Maker):

    name: str = "PV_DOC_FROM_MD_DOC"
    pv_extractor: PVExtractor = None

    @job
    def make(self, md_calc_output):
        return MDPVDataDoc(
            pressure=self.pv_extractor.get_pressure(md_calc_output),
            volume=self.pv_extractor.get_volume(md_calc_output)
        )

@dataclass
class PVFromMDFlowMaker(Maker):

    name: str = "MD_TO_PV"
    md_maker: Maker = None
    extract_maker: GetPVDocFromMDMaker = None

    def make(self, structure):
        md_job = self.md_maker.make(structure)
        extract_job = self.extract_maker.make(md_job.output)

        return Flow([md_job, extract_job], output=extract_job.output)


@dataclass
class EquilibriumVolumeSearchMaker(Maker):
    """Iteratively identifies the equilibrium volume
    of a structure at a particular temperature by fitting a growing
    P-V dataset using the Birch Murnaghan equation of state method.
    """

    name: str = "EQUIL_VOL_SEARCH"
    pv_from_md_maker: PVFromMDFlowMaker = None
    scale_factor_increment: float = 0.2

    @job
    def make(
        self,
        original_structure: Structure,
        pv_data_docs: List[MDPVDataDoc] = None
    ):
        if pv_data_docs is not None and len(pv_data_docs) > MAX_MD_JOBS:
            raise RuntimeError(
                "Maximum number of MD runs for equilibrium volume search exceeded"
            )

        if pv_data_docs is None:
            initial_scale_factors = [1 - self.scale_factor_increment, 1, 1 + self.scale_factor_increment]

            scaled_structs = [
                get_scaled_structure(original_structure, factor)
                for factor in initial_scale_factors
            ]

            new_jobs = [
                self.pv_from_md_maker.make(struct)
                for struct in scaled_structs
            ]

            pv_data_docs = [job.output for job in new_jobs]
        else:
            volumes = [doc.volume for doc in pv_data_docs]
            pressures = [doc.pressure for doc in pv_data_docs]
            pv_pairs = np.array(list(zip(pressures, volumes)))

            max_explored_volume = max(volumes)
            min_explored_volume = min(volumes)

            new_job_vol_scales = []
            try:
                params = rescale_volume.fit_BirchMurnaghanPV_EOS(pv_pairs)
                equil_volume = params[0]
                if (
                    equil_volume < max_explored_volume
                    and equil_volume > min_explored_volume
                ):
                    final_structure = original_structure.copy()
                    final_structure.scale_lattice(equil_volume)
                    return final_structure
                elif equil_volume > max_explored_volume:
                    new_job_vol_scales.append(get_new_max_volume(equil_volume, original_structure))
                elif equil_volume < min_explored_volume:
                    new_job_vol_scales.append(get_new_min_volume(equil_volume, original_structure))
            except ValueError:
                print("Unable to converge EoS fit for volume optimization, expanding search range")
                new_job_max = expand_upper_bound(max_explored_volume, original_structure)
                new_job_min = expand_lower_bound(max_explored_volume, original_structure)
                new_job_vol_scales.append(new_job_max)
                new_job_vol_scales.append(new_job_min)
        
            # This is specific to the type of MD run you're doing
            scaled_structs = [
                get_scaled_structure(original_structure, factor)
                for factor in new_job_vol_scales
            ]

            new_jobs =  [self.pv_from_md_maker.make(struct) for struct in scaled_structs]

            for new_job in new_jobs:
                pv_data_docs.append(new_job.output)

        expanded_search_job = self.make(original_structure, pv_data_docs)

        flow = Flow([*new_jobs, expanded_search_job])

        return Response(replace=flow, output=expanded_search_job.output)


def get_new_max_volume(equil_guess, original_structure):
    return equil_guess / original_structure.volume + BUFFER

def expand_upper_bound(old_max_vol, original_structure):
    return old_max_vol / original_structure.volume + 0.2

def get_new_min_volume(equil_guess, original_structure):
    return equil_guess / original_structure.volume - BUFFER

def expand_lower_bound(old_max_vol, original_structure):
    return old_max_vol / original_structure.volume - 0.2
