from pydantic import BaseModel, Field
from datetime import datetime
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.trajectory import Trajectory as AseTrajectory
from pymatgen.core.trajectory import Trajectory as PmgTrajectory


class M3GNetMDCalculation(BaseModel):
    task_label: str = Field(None, description="The name of the task.")
    last_updated: datetime = Field(
        description="Timestamp of when the document was last updated.",
    )
    trajectory: PmgTrajectory = Field(
        None, description="The pymatgen Trajectory object stored ad dictionary")
    metadata: dict = Field(
        None, description="The key info of the MD run like the temperature and time step")

    @classmethod
    def from_directory(cls, trajectory, metadata: dict, **kwargs):

        """
        Create a NetworkTask from the reaction network workflow outputs.
        Args:
            trajectory: the trajectory file saved from  MolecularDynamics to out.traj
            metadata: The metadata dictionary. like the temperature and timestep of the MD run.
            **kwargs: Additional keyword arguments to pass to the M3GNetMDCalculation.
        """

        traj = AseTrajectory(trajectory)
        structure_list = []
        for i, atoms in enumerate(traj[:]):  # first MD run
            structure = AseAtomsAdaptor.get_structure(atoms)
            structure_list.append(structure)

        traj_pmg = PmgTrajectory.from_structures(structure_list)
        traj_dict = traj_pmg.as_dict()

        d = {"trajectory": traj_dict, 'metadata': metadata}

        return cls(**d, **kwargs)
