import json
import os
import re
import zlib

import gridfs
import numpy as np
from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk, get_logger
from atomate.vasp.drones import VaspDrone
from bson import ObjectId
from mpmorph.database import VaspMDCalcDb
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from monty.json import MontyEncoder
from pymatgen import Structure
from pymatgen.core.trajectory import Trajectory
from collections import defaultdict

__author__ = 'Eric Sivonxay and Jianli Cheng'

logger = get_logger(__name__)


@explicit_serialize
class VaspMDToDb(FiretaskBase):
    """
    Enter a VASP run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.
    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains VASP
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        parse_dos (bool): whether to parse the DOS and store in GridFS.
            Defaults to False.
        bandstructure_mode (str): Set to "uniform" for uniform band structure.
            Set to "line" for line mode. If not set, band structure will not
            be parsed.
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        defuse_unsuccessful (bool): Defuses children fireworks if VASP run state
            is not "successful"; i.e. both electronic and ionic convergence are reached.
            Defaults to True.
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos", "bandstructure_mode",
                       "additional_fields", "db_file", "fw_spec_field",
                       "md_structures", "defuse_unsuccessful"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", False),
                          parse_bader=False,
                          bandstructure_mode=self.get("bandstructure_mode", False),
                          store_volumetric_data=[])

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = VaspMDCalcDb.from_db_file(db_file, admin=True)

            # prevent duplicate insertion
            mmdb.db.tasks.find_one_and_delete({'formula_pretty': task_doc['formula_pretty'],
                                               'task_label': task_doc['task_label']})

            t_id = mmdb.insert_task(task_doc,
                                    parse_dos=self.get("parse_dos", False),
                                    parse_bs=bool(self.get("bandstructure_mode", False)),
                                    md_structures=self.get("md_structures", True))

            logger.info("Finished parsing with task_id: {}".format(t_id))

        if self.get("defuse_unsuccessful", True):
            defuse_children = (task_doc["state"] != "successful")
        else:
            defuse_children = False

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children)


@explicit_serialize
class TrajectoryDBTask(FiretaskBase):
    """
    Obtain all production runs and insert them into the db. This is done by
    searching for a unique tag
    """
    required_params = ["tag_id", "db_file"]
    optional_params = ['notes']

    def run_task(self, fw_spec):
        notes = self.get('notes', None)
        tag_id = self['tag_id']

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)
        mmdb = VaspMDCalcDb.from_db_file(db_file, admin=True)
        # mmdb.db.trajectories.find_one_and_delete({"runs_label": tag_id})
        runs = mmdb.db['tasks'].find(
            {"task_label": re.compile(f'.*prod_run.*{tag_id}.*')})

        runs_sorted = sorted(runs, key=lambda x: int(re.findall('run[_-](\d+)', x['task_label'])[0]))

        # Remove duplicates of the same run (if they exist)
        labels = [result['task_label'] for result in runs_sorted]
        nums = [int(re.findall('run[_-](\d+)', label)[0]) for label in labels]
        duplicates = np.where((nums - np.roll(nums, 1)) == 0)[0]
        runs_sorted = [runs_sorted[i] for i in range(len(runs_sorted)) if i not in duplicates]

        trajectory_doc = runs_to_trajectory_doc(runs_sorted, db_file, tag_id, notes)

        mmdb.db.trajectories.insert_one(trajectory_doc)


def runs_to_trajectory_doc(runs, db_file, runs_label, notes=None):
    """
    Takes a list of task_documents, aggregates the trajectories from the ionics_steps gridfs storage, then dumps
    the pymatgen.core.Trajectory object into the 'trajectories_fs' collection and makes a dictionary doc to track
    the entry.

    :param runs: list of MD runs
    :param db_file:
    :param runs_label: unique identifier to the runs
    :param notes: (optional) any notes or comments on the specific run
    :return:
    """
    mmdb = VaspMDCalcDb.from_db_file(db_file, admin=True)

    trajectory = load_trajectories_from_gfs(runs, mmdb)

    traj_dict = json.dumps(trajectory.as_dict(), cls=MontyEncoder)
    gfs_id, compression_type = insert_gridfs(traj_dict, mmdb.db, "trajectories_fs")

    traj_doc = {
        'formula_pretty': trajectory[0].composition.reduced_formula,
        'formula': trajectory[0].composition.formula.replace(' ', ''),
        'temperature': int(runs[0]["input"]["incar"]["TEBEG"]),
        'runs_label': runs_label,
        'compression': compression_type,
        'fs_id': gfs_id,
        'step_fs_ids': [i["calcs_reversed"][0]["output"]['ionic_steps_fs_id']
                        for i in runs],
        'structure': trajectory[0].as_dict(),
        'dimension': list(np.shape(trajectory.frac_coords)),
        'time_step': runs[0]["input"]["incar"]["POTIM"] * 1e-3,
        'notes': notes
    }
    return traj_doc


def load_trajectories_from_gfs(runs, mmdb, gfs_keys=None):
    if gfs_keys is None:
        # Attempt to automatically determine where the trajectory is stored
        gfs_keys = []
        for run in runs:
            # 3 cases to deal with: 1) Trajectory 2) previous_runs (old mpmorph) 3) structures_fs
            if 'trajectory' in run.keys():
                gfs_keys.append((run['trajectory']['fs_id'], 'trajectories_fs'))
            elif "INCAR" in run.keys():
                # for backwards compatibility with older version of mpmorph
                gfs_keys.append((run["ionic_steps_fs_id"], 'previous_runs_gfs'))
            elif "input" in run.keys():
                gfs_keys.append((run["calcs_reversed"][0]["output"]["ionic_steps_fs_id"], 'structures_fs'))

    trajectory = None
    for i, (fs_id, fs) in enumerate(gfs_keys):

        if fs == 'trajectories_fs' or fs == 'rebuild_trajectories_fs':
            # Load stored Trajectory
            print(fs_id, 'is stored in trajectories_fs')
            _trajectory = load_trajectory(fs_id=fs_id, db=mmdb.db, fs=fs)
        else:
            print(fs_id)
            # This is compatibility code for when
            # Load Ionic steps from gfs, then convert to trajectory before extending (compatibility code)
            ionic_steps_dict = load_ionic_steps(fs_id=fs_id, db=mmdb.db, fs=fs)
            ionic_steps_defaultdict = defaultdict(list)
            for d in ionic_steps_dict:
                for key, val in d.items():
                    ionic_steps_defaultdict[key].append(val)
            ionic_steps = dict(ionic_steps_defaultdict.items())

            # extract structures from dictionary
            structures = [Structure.from_dict(struct) for struct in ionic_steps['structure']]
            del ionic_steps['structure']

            frame_properties = {}
            for key in ['e_fr_energy', 'e_wo_entrp', 'e_0_energy', 'kinetic', 'lattice kinetic', 'nosepot',
                        'nosekinetic', 'total']:
                frame_properties[key] = ionic_steps[key]

            # Create trajectory
            _trajectory = Trajectory.from_structures(structures, constant_lattice=True,
                                                     frame_properties=frame_properties,
                                                     time_step=run[0]['input']['incar']['POTIM'])
        if trajectory is None:
            trajectory = _trajectory
        else:
            # Eliminate duplicate structure at the start of each trajectory
            # (since vasp will output the input structure)
            trajectory.extend(_trajectory[1:])
    return trajectory


def process_traj(data):
    i, fs_id, fs, db_file = data[0], data[1], data[2], data[3]
    mmdb = VaspMDCalcDb.from_db_file(db_file, admin=True)
    ionic_steps_dict = load_ionic_steps(fs_id, mmdb.db, fs)

    structure = Structure.from_dict(ionic_steps_dict[0]['structure'])
    positions = [0] * len(ionic_steps_dict)
    for i, step in enumerate(ionic_steps_dict):
        _step = [atom['abc'] for atom in step["structure"]["sites"]]
        positions[i] = _step

    traj = Trajectory(structure.lattice.matrix, structure.species, positions, 0.002)
    return i, traj.as_dict()


def load_trajectory(fs_id, db, fs=None):
    if not fs:
        # Default to trajectories_fs
        fs = gridfs.GridFS(db, 'trajectories_fs')
    elif not isinstance(fs, gridfs.GridFS):
        # Handle fs supplied as str
        fs = gridfs.GridFS(db, fs)

    trajectories_json = zlib.decompress(fs.get(fs_id).read())
    trajectories_dict = json.loads(trajectories_json.decode())
    try:
        trajectory = Trajectory.from_dict(trajectories_dict)
    except AttributeError:
        trajectories_dict = json.loads(trajectories_dict)
        trajectory = Trajectory.from_dict(trajectories_dict)
    return trajectory


def load_ionic_steps(fs_id, db, fs):
    if not isinstance(fs, gridfs.GridFS):
        fs = gridfs.GridFS(db, fs)
    ionic_steps_json = zlib.decompress(fs.get(fs_id).read())
    ionic_steps_dict = json.loads(ionic_steps_json.decode())
    del ionic_steps_json
    return ionic_steps_dict


def insert_gridfs(d, db, collection="fs", compress=True, oid=None, task_id=None):
    """
    Insert the given document into GridFS.
    Args:
        d (dict): the document
        collection (string): the GridFS collection name
        compress (bool): Whether to compress the data or not
        oid (ObjectId()): the _id of the file; if specified, it must not already exist in GridFS
        task_id(int or str): the task_id to store into the gridfs metadata
    Returns:
        file id, the type of compression used.
    """
    oid = oid or ObjectId()
    compression_type = None

    if compress:
        d = zlib.compress(d.encode(), compress)
        compression_type = "zlib"

    fs = gridfs.GridFS(db, collection)
    if task_id:
        # Putting task id in the metadata subdocument as per mongo specs:
        # https://github.com/mongodb/specifications/blob/master/source/gridfs/gridfs-spec.rst#terms
        fs_id = fs.put(d, _id=oid, metadata={"task_id": task_id, "compression": compression_type})
    else:
        fs_id = fs.put(d, _id=oid, metadata={"compression": compression_type})

    return fs_id, compression_type
