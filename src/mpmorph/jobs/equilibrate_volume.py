from dataclasses import dataclass
from typing import List
import numpy as np

from pymatgen.core.structure import Structure
from jobflow import Maker, job, Flow, Response
from atomate2.vasp.jobs.core import MDMaker

from ..flows.scale_volume import md_to_volume_flow
from .extract_pv_vasp import ExtractPVDataFromVASPMDMaker
from ..schemas.pv_data_doc import MDPVDataDoc
from ..runners import rescale_volume


MAX_MD_JOBS = 5 # if you can't converge with five additional calcs you're doing something wrong...
OFFSET = 0.1 # gives it enough room to slosh back

@dataclass
class EquilibriumVolumeSearchMaker(Maker):
    """Iteratively identifies the equilibrium volume
    of a structure at a particular temperature by fitting a growing
    P-V dataset using the Birch Murnaghan equation of state method.
    """

    name: str = "EQUIL_VOL_SEARCH"
    md_maker: Maker = MDMaker()
    pv_maker: Maker = ExtractPVDataFromVASPMDMaker()

    @job
    def make(self, original_structure: Structure, md_pv_data_docs: List[MDPVDataDoc]):
        if len(md_pv_data_docs) > MAX_MD_JOBS:
            raise RuntimeError("Maximum number of jobs for equilibrium volume search exceeded")


        volumes = [doc.volume for doc in md_pv_data_docs]
        pressures = [doc.pressure for doc in md_pv_data_docs]
        vp_pairs = np.array(list(zip(volumes, pressures)))

        max_explored_volume = max(volumes)
        min_explored_volume = min(volumes)

        params = rescale_volume.fit_BirchMurnaghanPV_EOS(vp_pairs)
        equil_volume = params[0]
        
        if equil_volume < max_explored_volume and equil_volume > min_explored_volume:
            final_structure = original_structure.copy()
            return final_structure.scale_lattice(equil_volume)

        elif equil_volume > max_explored_volume: 
            new_vol_scale = get_new_max_volume(equil_volume, original_structure)

        elif equil_volume < min_explored_volume: 
            new_vol_scale = get_new_min_volume(equil_volume, original_structure)

        # This is specific to the type of MD run you're doing
        new_job = md_to_volume_flow(
            original_structure,
            new_vol_scale,
            self.md_maker,
            self.pv_maker
        )

        md_pv_data_docs.append(new_job.output)
        expanded_search_job = EquilibriumVolumeSearchMaker(
            md_maker = self.md_maker,
            pv_maker = self.pv_maker,
        ).make(original_structure, md_pv_data_docs)

        flow = Flow([new_job, expanded_search_job])

        return Response(replace = flow)

def get_new_max_volume(equil_guess, original_structure):
    return equil_guess / original_structure.volume + OFFSET

def get_new_min_volume(equil_guess, original_structure):
    return equil_guess / original_structure.volume - OFFSET


