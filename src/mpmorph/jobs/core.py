from dataclasses import dataclass, field
from jobflow import Maker, job
from pymatgen.core import Structure

from .tasks.m3gnet_input import M3GNetMDInputs
from .tasks.chgnet_input import CHGNetMDInputs
from .tasks.m3gnet_md_task import run_m3gnet
from .tasks.chgnet_md_task import run_chgnet
from ..schemas.m3gnet_md_calc import M3GNetMDCalculation
from ..schemas.chgnet_md_calc import CHGNetMDCalculation


def empty_inputs():
    return M3GNetMDInputs()

def empty_inputs_chgnet():
    return CHGNetMDInputs()


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
    parameters: M3GNetMDInputs = field(
        default_factory=empty_inputs,
    )

    @job(trajectory="trajectory", output_schema=M3GNetMDCalculation)
    def make(self, structure: Structure, **kwargs):
        """
        Run MD using the M3GNet Molecular Dynamics interface. This runs molecular
        dynamics through the ASE interface with the default M3GNet potential.

        Args:
            structure: the input structure
        """

        calc_doc = run_m3gnet(structure, self.parameters, self.name, **kwargs)

        return calc_doc


@dataclass
class CHGNetMDMaker(Maker):
    """
    Maker to create CHGNet MD calculation jobs.

    Args:
        name: The name of the job. Defaults to "chgnet_run".
        parameters: An CHGNetMDInputs object containing the parameters for
            the MD run.
    """

    name: str = "chgnet_run"
    parameters: CHGNetMDInputs = field(
        default_factory=empty_inputs_chgnet,
    )

    @job(trajectory="trajectory", output_schema=CHGNetMDCalculation)
    def make(self, structure: Structure, **kwargs):
        """
        Run MD using the M3GNet Molecular Dynamics interface. This runs molecular
        dynamics through the ASE interface with the default M3GNet potential.

        Args:
            structure: the input structure
        """

        calc_doc = run_chgnet(structure, self.parameters, self.name, **kwargs)

        return calc_doc