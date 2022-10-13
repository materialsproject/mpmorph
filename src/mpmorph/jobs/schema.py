from datetime import datetime
from pathlib import Path
from typing import Union

from ase.io.trajectory import Trajectory as AseTrajectory
from pydantic import BaseModel, Field
from pymatgen.core.trajectory import Trajectory as PmgTrajectory
from pymatgen.io.ase import AseAtomsAdaptor

from mpmorph.utils import datetime_str


class M3GNetCalculation(BaseModel):
    task_label: str = Field(None, description="The name of the task.")
    last_updated: datetime = Field(
        default_factory=datetime_str,
        description="Timestamp of when the document was last updated.",
    )
    trajectory: PmgTrajectory = Field(
        None, description="The pymatgen Trajectory object stored ad dictionary"
    )
    metadata: dict = Field(
        None,
        description="The key info of the MD run like the temperature and time step",
    )

    @classmethod
    def from_directory(
        cls,
        dir_name: Union[Path, str],
        trajectory_fn="out.traj",
    ):

        """
        Create a M3GnetCalculation document from a directory containing output files of
        a M3GNet MD run.
        """
        dir_name = Path(dir_name)

        trajectory = AseTrajectory(filename=dir_name / trajectory_fn)

        return cls.from_trajectory(trajectory)

    @classmethod
    def from_trajectory(cls, trajectory: AseTrajectory, **kwargs):
        """
        Args:
            trajectory: the trajectory file saved from  MolecularDynamics to out.traj
            metadata: The metadata dictionary. like the temperature and timestep of the MD run.
            **kwargs: Additional keyword arguments to pass to the M3GNetCalculation.
        """
        traj = AseTrajectory(trajectory)
        structure_list = []
        for _, atoms in enumerate(traj[:]):  # first MD run
            structure = AseAtomsAdaptor.get_structure(atoms)
            structure_list.append(structure)

        traj_pmg = PmgTrajectory.from_structures(structure_list)
        traj_dict = traj_pmg.as_dict()

        metadata = {}

        d = {"trajectory": traj_dict, "metadata": metadata}

        return cls(**d, **kwargs)
