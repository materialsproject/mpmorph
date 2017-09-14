from atomate.utils.database import CalcDb
import yaml
from pymongo import MongoClient
from bson import ObjectId
import gridfs
import zlib
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.outputs import Vasprun
from monty.json import MontyEncoder
import json
import yaml
from pymongo import MongoClient
import shutil
import os


def archive_system(directory):
    osw = os.walk(directory)
    dirs = []
    i = 0
    for _dir in osw:
        if len(_dir[1]) == 0 and "POSCAR" in _dir[2]:
            dirs.append(_dir[0])
            i += 1

    with open("/global/homes/s/sivonxay/.conda/envs/knl_env/config/my_launchpad.yaml", "r") as f:
        t = yaml.load(f)
    client = MongoClient(t['host'])
    db = client.fw_es_vasp
    db.authenticate(t['username'], t['password'])

    for _dir in dirs:
        task_doc = {}
        if os.path.exists(_dir + "/POSCAR"):
            _poscar = Poscar.from_file(_dir + "/POSCAR")
        elif os.path.exists(_dir + "/POSCAR.gz"):
            _poscar = Poscar.from_file(_dir + "/POSCAR")
        else:
            continue
        task_doc["POSCAR"] = _poscar.as_dict()

        if os.path.exists(_dir + "/INCAR"):
            _incar = Incar.from_file(_dir + "/INCAR")
        elif os.path.exists(_dir + "/INCAR.gz"):
            _incar = Incar.from_file(_dir + "/INCAR")
        else:
            continue
        task_doc["INCAR"] = _incar.as_dict()
        task_doc["task_label"] = get_label(_dir)
        if get_snap(_dir) != None:
            task_doc["snap"] = get_snap(_dir)
        task_doc["chemsys"] = [str(el) for el in _poscar.structure.composition.elements]
        task_doc["nsites"] = _poscar.structure.num_sites
        task_doc["formula_pretty"] = _poscar.structure.composition.reduced_formula
        task_doc["nelements"] = len(_poscar.structure.composition.elements)
        if doc_exists(task_doc, db):
            print(get_label(_dir) + " exists")
            print("deleting " + get_label(_dir))
            shutil.rmtree(_dir)
            continue
        if "diffusion" in _dir:
            task_doc["diffusion"] = True
        if "longrun" in _dir or "diffusion_run" in _dir:
            print(_dir)
            # load vasprun.xml and save the structures
            if os.path.exists(_dir + "/vasprun.xml"):
                _vr = Vasprun(_dir + "/vasprun.xml")
            elif os.path.exists(_dir + "/vasprun.xml.gz"):
                _vr = Vasprun(_dir + "/vasprun.xml.gz")
            ionic_steps = json.dumps(_vr.as_dict()["output"]["ionic_steps"], cls=MontyEncoder)
            gfs_id, compression_type = insert_gridfs(ionic_steps, "previous_runs_gfs")
            task_doc["ionic_steps_fs_id"] = gfs_id
            task_doc["ionic_steps_compression"] = compression_type
        print("inserting " + get_label(_dir))
        insert(task_doc, db)
        print("deleting " + get_label(_dir))
        shutil.rmtree(_dir)
    print("-" * 100)
    prune_dirs(directory)

def get_label(directory):
    index = directory.rfind("/")
    label = directory[index+1:]
    return label

def get_snap(directory):
    index = directory.rfind("snap")
    if index != -1:
        snap = directory[index+5]
    else:
        return None
    return snap

def doc_exists(task_doc, db):
    collection = db["previous_runs"]
    if not collection.find_one({"POSCAR": task_doc["POSCAR"], \
                                "INCAR": task_doc["INCAR"], "task_label": task_doc["task_label"]}):
        return False
    return True

def insert(task_doc, db):
    collection = db["previous_runs"]
    if not collection.find_one({"POSCAR": task_doc["POSCAR"], "INCAR": task_doc["INCAR"]}):
        collection.insert_one(task_doc)

def insert_gridfs(d, collection="fs", compress=True, oid=None):
        """
        Insert the given document into GridFS.
        Args:
            d (dict): the document
            collection (string): the GridFS collection name
            compress (bool): Whether to compress the data or not
            oid (ObjectId()): the _id of the file; if specified, it must not already exist in GridFS
        Returns:
            file id, the type of compression used.
        """
        with open("/global/homes/s/sivonxay/.conda/envs/knl_env/config/my_launchpad.yaml", "r") as f:
            t = yaml.load(f)
        client = MongoClient(t['host'])
        db = client.fw_es_vasp
        db.authenticate(t['username'], t['password'])

        oid = oid or ObjectId()
        if compress:
            d = zlib.compress(d.encode(), compress)
        fs = gridfs.GridFS(db, collection)
        fs_id = fs.put(d, _id=oid)
        return fs_id, "zlib"

def prune_dirs(directory, change=True):
    _change = False
    if change == False:
        return
    else:
        for _dir in os.walk(directory):
            if len(_dir[1]) == 0 and len(_dir[2]) == 0:
                os.rmdir(_dir[0])
                _change=True
    prune_dirs(directory, change=_change)
    return