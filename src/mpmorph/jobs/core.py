from dataclasses import dataclass
from jobflow import Maker, job
from pymatgen.core import Structure

from .tasks.m3gnet_input import M3GNetMDInputs
from .tasks.m3gnet_md_task import run_m3gnet
from ..schemas.m3gnet_md_calc import M3GNetMDCalculation

@dataclass
class M3GNetMDMaker(Maker):
    """
    Maker to create M3GNet MD calculation jobs.

    Args:
        name: The name of the job. Defaults to "m3gnet_run".
        parameters: An M3GNetMDInputs object containing the parameters for
            the MD run.
    """

    name: str = "m3gnet_run"
    parameters: M3GNetMDInputs = None

    @job(trajectory="trajectory", output_schema=M3GNetMDCalculation)
    def make(self, structure: Structure, **kwargs):
        """
        Run MD using the M3GNet Molecular Dynamics interface. This runs molecular
        dynamics through the ASE interface with the default M3GNet potential.

        Args:
            structure: the input structure
        """

        calc_doc = run_m3gnet(
            structure,
            self.parameters,
            self.name,
            **kwargs
        )

        return calc_doc
