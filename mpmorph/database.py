from __future__ import division, print_function, unicode_literals, absolute_import

import json
import zlib

import gridfs
from atomate.utils.utils import get_logger
from atomate.vasp.database import VaspCalcDb
from monty.json import MontyEncoder
from pymatgen.core import Structure
from pymatgen.core.trajectory import Trajectory
from collections import defaultdict
import numpy as np

logger = get_logger(__name__)

__author__ = 'Eric Sivonxay <esivonxay@lbl.gov>'


class VaspMDCalcDb(VaspCalcDb):
    """
    Adapted from atomate.vasp.database

    Class to help manage database insertions of Vasp drones
    """

    def __init__(self, host="localhost", port=27017, database="vasp", collection="tasks", user=None,
                 password=None, **kwargs):
        super(VaspMDCalcDb, self).__init__(host, port, database, collection, user, password, **kwargs)

    def insert_task(self, task_doc, parse_dos=False, parse_bs=False, md_structures=False):
        """
        Inserts a task document (e.g., as returned by Drone.assimilate()) into the database.
        Handles putting DOS and band structure into GridFS as needed.
        Args:
            task_doc: (dict) the task document
            parse_dos: (bool) attempt to parse dos in task_doc and insert into Gridfs
            parse_bs: (bool) attempt to parse bandstructure in task_doc and insert into Gridfs
        Returns:
            (int) - task_id of inserted document
        """

        # insert dos into GridFS
        if parse_dos and "calcs_reversed" in task_doc:
            if "dos" in task_doc["calcs_reversed"][0]:  # only store idx=0 DOS
                dos = json.dumps(task_doc["calcs_reversed"][0]["dos"], cls=MontyEncoder)
                gfs_id, compression_type = self.insert_gridfs(dos, "dos_fs")
                task_doc["calcs_reversed"][0]["dos_compression"] = compression_type
                task_doc["calcs_reversed"][0]["dos_fs_id"] = gfs_id
                del task_doc["calcs_reversed"][0]["dos"]

        # insert band structure into GridFS
        if parse_bs and "calcs_reversed" in task_doc:
            if "bandstructure" in task_doc["calcs_reversed"][0]:  # only store idx=0 BS
                bs = json.dumps(task_doc["calcs_reversed"][0]["bandstructure"], cls=MontyEncoder)
                gfs_id, compression_type = self.insert_gridfs(bs, "bandstructure_fs")
                task_doc["calcs_reversed"][0]["bandstructure_compression"] = compression_type
                task_doc["calcs_reversed"][0]["bandstructure_fs_id"] = gfs_id
                del task_doc["calcs_reversed"][0]["bandstructure"]

        # insert structures  at each ionic step into GridFS
        if md_structures and "calcs_reversed" in task_doc:
            # insert ionic steps into gridfs
            # #TODO: Deprecate this and move to only storing trajectory
            # ionic_steps_json = json.dumps(task_doc["calcs_reversed"][0]['output']['ionic_steps'], cls=MontyEncoder)
            # gfs_id, compression_type = self.insert_gridfs(ionic_steps_json, "structures_fs")
            # task_doc["calcs_reversed"][0]['output']['ionic_steps_compression'] = compression_type
            # task_doc["calcs_reversed"][0]['output']['ionic_steps_fs_id'] = gfs_id

            # Aggregate a trajectory
            ## Convert from a list of dictionaries to a dictionary of lists
            ionic_steps_dict = task_doc["calcs_reversed"][0]['output']['ionic_steps']
            del task_doc["calcs_reversed"][0]['output']['ionic_steps']
            ionic_steps_defaultdict = defaultdict(list)
            for d in ionic_steps_dict:
                for key, val in d.items():
                    ionic_steps_defaultdict[key].append(val)
            ionic_steps = dict(ionic_steps_defaultdict.items())

            ## extract structures from dictionary
            structures = [Structure.from_dict(struct) for struct in ionic_steps['structure']]
            del ionic_steps['structure']

            frame_properties = {}
            for key in ['e_fr_energy', 'e_wo_entrp', 'e_0_energy', 'kinetic', 'lattice kinetic', 'nosepot',
                        'nosekinetic', 'total']:
                frame_properties[key] = ionic_steps[key]

            # Create trajectory
            trajectory = Trajectory.from_structures(structures, constant_lattice=True,
                                                    frame_properties=frame_properties,
                                                    time_step=task_doc['input']['incar']['POTIM'])
            traj_dict = json.dumps(trajectory, cls=MontyEncoder)
            gfs_id, compression_type = self.insert_gridfs(traj_dict, "trajectories_fs")

            task_doc['trajectory'] = {
                'formula': trajectory[0].composition.formula.replace(' ', ''),
                'temperature': int(task_doc["input"]["incar"]["TEBEG"]),
                'compression': compression_type,
                'fs_id': gfs_id,
                'dimension': list(np.shape(trajectory.frac_coords)),
                'time_step': task_doc["input"]["incar"]["POTIM"],
            }

        # insert the task document and return task_id
        return self.insert(task_doc)

    # def get_ionic_steps(self, task_id):
    #     m_task = self.collection.find_one({"task_id": task_id}, {"calcs_reversed": 1})
    #     fs_id = m_task["calcs_reversed"][0]["output"]["ionic_steps_fs_id"]
    #     fs = gridfs.GridFS(self.db, "structures_fs")
    #     ionic_steps_json = zlib.decompress(fs.get(fs_id).read())
    #     ionic_steps_dict = json.loads(ionic_steps_json.decode())
    #     return ionic_steps_dict
    #
    # def get_structures(self, task_id):
    #     ionic_steps_dict = self.get_ionic_steps(task_id)
    #     structures = [Structure.from_dict(step["structure"]) for step in ionic_steps_dict]
    #     return structures
