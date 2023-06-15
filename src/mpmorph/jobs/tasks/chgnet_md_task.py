from pathlib import Path

import os

from ase import units
from chgnet.model.dynamics import MolecularDynamics
from pymatgen.core import Structure
import dataclasses

from mpmorph.jobs.tasks.chgnet_input import CHGNetMDInputs
from ...schemas.chgnet_md_calc import CHGNetMDCalculation

def run_chgnet(
    structure: Structure, inputs: CHGNetMDInputs, name: str = "chgnet_run", **kwargs
):
    """
    Run MD using the CHGNet Molecular Dynamics interface. This runs molecular
    dynamics through the ASE interface with the default CHGNet potential.

    Args:
        structure: the input structure
        inputs: The parameters for the MD run
        name: The name of the calculation document returned by this simulation
    """

    outfile_name = f"{structure.composition.to_pretty_string()}-{inputs.temperature}K-{round(structure.volume, 3)}"
    traj_fn = f"{outfile_name}.traj"
    log_fn = f"{outfile_name}.log"

    taut = 10 * units.fs

    if inputs.potential is not None:
        kwargs["potential"] = inputs.potential

    md = MolecularDynamics(
        atoms=structure,
        ensemble=inputs.ensemble,
        temperature=inputs.temperature,
        timestep=inputs.timestep,
        pressure=inputs.pressure,
        use_device=inputs.use_device,
        taut=taut,
        taup=inputs.taup,
        compressibility_au=inputs.compressibility_au,
        trajectory=traj_fn,
        logfile=log_fn,
        loginterval=inputs.loginterval,
        append_trajectory=inputs.append_trajectory,
        **kwargs,
    )

    md.run(
        steps=inputs.steps,
    )

    d = CHGNetMDCalculation.from_directory(Path.cwd(), trajectory_fn=traj_fn)

    if not inputs.save_files:
        os.remove(traj_fn)
        os.remove(log_fn)

    d.task_label = name
    d.metadata = dataclasses.asdict(inputs)

    return d
