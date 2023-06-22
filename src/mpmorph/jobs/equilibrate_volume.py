from dataclasses import dataclass
from typing import List, Tuple
import numpy as np

from atomate2.vasp.jobs.core import MDMaker

from pymatgen.core.structure import Structure
from jobflow import Maker, job, Flow, Response

from .pv_from_calc import PVFromVasp, PVExtractor
from ..schemas.pv_data_doc import MDPVDataDoc
from ..runners import rescale_volume


MAX_MD_JOBS = 5  # if you can't converge with five additional calcs you're doing something wrong...
OFFSET = 0.1  # gives it enough room to slosh back

def get_scaled_structure(struct: Structure, scale_factor: float):
    copy: Structure = struct.copy()
    copy.scale_lattice(copy.volume * scale_factor)
    return copy

@dataclass
class EquilibriumVolumeSearchMaker(Maker):
    """Iteratively identifies the equilibrium volume
    of a structure at a particular temperature by fitting a growing
    P-V dataset using the Birch Murnaghan equation of state method.
    """

    name: str = "EQUIL_VOL_SEARCH"
    md_maker: Maker = None
    pv_extractor: PVExtractor = None
    initial_scale_factors: Tuple[float] = (0.8, 1, 1.2)

    @job
    def make(
        self,
        original_structure: Structure,
        md_calc_outputs: List[MDPVDataDoc] = None
    ):
        if md_calc_outputs is not None and len(md_calc_outputs) > MAX_MD_JOBS:
            raise RuntimeError(
                "Maximum number of jobs for equilibrium volume search exceeded"
            )

        if md_calc_outputs is None:
            scaled_structs = [
                get_scaled_structure(original_structure, factor)
                for factor in self.initial_scale_factors
            ]

            new_jobs = [
                self.md_maker.make(struct)
                for struct in scaled_structs
            ]

            md_calc_outputs = [job.output for job in new_jobs]
        else:
            volumes = [self.pv_extractor.get_volume(doc) for doc in md_calc_outputs]
            pressures = [self.pv_extractor.get_pressure(doc) for doc in md_calc_outputs]
            pv_pairs = np.array(list(zip(pressures, volumes)))

            max_explored_volume = max(volumes)
            min_explored_volume = min(volumes)

            try:
                params = rescale_volume.fit_BirchMurnaghanPV_EOS(pv_pairs)
                equil_volume = params[0]
            except ValueError:
                return MDPVDataDoc()
            

            if (
                equil_volume < max_explored_volume
                and equil_volume > min_explored_volume
            ):
                final_structure = original_structure.copy()
                final_structure.scale_lattice(equil_volume)
                return final_structure

            elif equil_volume > max_explored_volume:
                new_vol_scale = get_new_max_volume(equil_volume, original_structure)

            elif equil_volume < min_explored_volume:
                new_vol_scale = get_new_min_volume(equil_volume, original_structure)

            # This is specific to the type of MD run you're doing
            new_job = self.md_maker.make(original_structure, new_vol_scale)
            new_jobs = [new_job]
            md_calc_outputs.append(new_job.output)

        expanded_search_job = EquilibriumVolumeSearchMaker(
            md_maker=self.md_maker,
            pv_extractor=self.pv_extractor
        ).make(original_structure, md_calc_outputs)

        flow = Flow([*new_jobs, expanded_search_job])

        return Response(replace=flow, output=expanded_search_job.output)


def get_new_max_volume(equil_guess, original_structure):
    return equil_guess / original_structure.volume + OFFSET


def get_new_min_volume(equil_guess, original_structure):
    return equil_guess / original_structure.volume - OFFSET
