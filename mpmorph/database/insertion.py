import os
import yaml
import pymongo

from fireworks import LaunchPad
from pymongo import MongoClient
from atomate.utils.utils import get_logger
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone

logger =get_logger(__name__)

class PostRunInsertion():
    """
    Inserts existing amorphous calculations into MongoDB without using fireworks
    """
    def __init__(self):
        with open("/global/homes/s/sivonxay/.conda/envs/knl_env/config/my_launchpad.yaml", "r") as f:
            t = yaml.load(f)
        client = MongoClient(t['host'])
        db = client.fw_es_vasp
        db.authenticate(t['username'], t['password'])
        self.tasks = db.tasks

        db_file = "/global/homes/s/sivonxay/.conda/envs/knl_env/config/db.json"
        self.mmdb = VaspCalcDb.from_db_file(db_file, admin=True)


    def run(self, directory):

        directories = os.walk(directory)
        dirs = []
        for dir in directories:
            if os.path.isfile(os.path.join(dir[0], "vasprun.xml")):
                dirs.append(dir[0])

        for calc_dir in dirs:
            #TODO Add functionality for md runs
            name = "vasp_run"
            if "relax" in calc_dir:
                task_doc = self.get_task_doc(calc_dir, name="static")
                self.insert_task_doc(task_doc=task_doc)
        return

    def get_task_doc(self, calc_dir, name):
        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))
        drone = VaspDrone(additional_fields={"task_label": name}, parse_dos=False, compress_dos=1,
                          bandstructure_mode=False, compress_bs=1)

        task_doc = drone.assimilate(calc_dir)
        return task_doc

    def insert_task_doc(self, task_doc):
        matches = list(self.tasks.find({'completed_at': task_doc['completed_at'], 'formula_pretty': task_doc['formula_pretty']}))
        if len(matches) == 0:
            t_id = self.mmdb.insert_task(task_doc, parse_dos=False, parse_bs=False)
            logger.info("Finished parsing with task_id: {}".format(t_id))
        else:
            print(len(matches))
            logger.info("Entry exists in db")
        return