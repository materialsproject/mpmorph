from pathlib import Path
from typing import Any, Dict, Iterable, List, Union

from ase.io.trajectory import Trajectory as AseTrajectory
from pydantic import BaseModel, Field
from pymatgen.core.composition import Composition
from pymatgen.core.trajectory import Trajectory as PmgTrajectory
from pymatgen.io.ase import AseAtomsAdaptor

from mpmorph.utils import datetime_str


class M3GNetMDCalculation(BaseModel):
    task_label: str = Field(None, description="The name of the task.")
    dir_name: str = Field(
        None, description="The directory where the M3GNet calculation was run"
    )
    last_updated: str = Field(
        default_factory=datetime_str,
        description="Timestamp of when the document was last updated.",
    )
    trajectory: PmgTrajectory = Field(
        None, description="The pymatgen Trajectory object stored ad dictionary"
    )
    composition: Composition = Field(description="The composition of the structure.")
    reduced_formula: str = Field(
        description="The reduced formula of the structure's composition."
    )
    metadata: dict = Field(
        None,
        description=(
            "Important info about the calculation, including ensemble type,"
            " temperature, etc."
        ),
    )

    @classmethod
    def from_directory(
        cls,
        dir_name: Union[Path, str],
        trajectory_fn: str = "out.traj",
        frame_properties: Iterable[str] = (
            "total_energy",
            "potential_energy",
            "kinetic_energy",
            "stress",
        ),
        **kwargs,
    ):
        """
        Create a M3GnetCalculation document from a directory containing output files of
        a M3GNet MD run.
        """
        dir_name = Path(dir_name)

        trajectory = AseTrajectory(filename=str(dir_name / trajectory_fn))

        return cls.from_trajectory(
            trajectory,
            frame_properties=frame_properties,
            dir_name=str(dir_name),
            **kwargs,
        )

    @classmethod
    def from_trajectory(
        cls,
        trajectory: AseTrajectory,
        frame_properties: Iterable[str] = (
            "total_energy",
            "potential_energy",
            "kinetic_energy",
            "stress",
            "temperature",
        ),
        **kwargs,
    ):
        """
        Create a M3GnetCalculation document from an ASE trajectory object.

        Args:
            trajectory: the ASE trajectory file loaded from the out.traj file
            **kwargs: Additional keyword arguments to pass to the M3GNetCalculation
                constructor.
        """
        structures = []
        frame_properties_list: List[Dict[str, Any]] = []

        initial_structure = AseAtomsAdaptor.get_structure(trajectory[0])

        for atoms in trajectory:
            struct = AseAtomsAdaptor.get_structure(atoms)
            frame_props = {k: getattr(atoms, f"get_{k}")() for k in frame_properties}

            structures.append(struct)
            frame_properties_list.append(frame_props)

        traj_pmg = PmgTrajectory.from_structures(
            structures,
            frame_properties=frame_properties_list,
            time_step=trajectory.description["timestep"],
        )

        d = {
            "trajectory": traj_pmg,
            "composition": initial_structure.composition,
            "reduced_formula": initial_structure.composition.reduced_formula,
        }

        return cls(**d, **kwargs)
