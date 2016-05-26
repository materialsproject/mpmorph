from fireworks import explicit_serialize, FireTaskBase, FWAction
from mpmorph.runners.amorphous_maker import AmorphousMaker

__author__ = 'Muratahan Aykol <maykol@lbl.gov>'

@explicit_serialize
class AmorphousMakerTask(FireTaskBase):
    """
    Create a constraint random packed structure from composition and box dimensions.
    Required params:
        composition: (dict) a dict of target composition with integer atom numbers
                        e.g. {"V":22, "Li":10, "O":75, "B":10}
        box_scale: (float) all lattice vectors are multiplied with this scalar value.
                        e.g. edge length of a cubic simulation box.
    Optional params:
        tol (float): tolerance factor for how close the atoms can get (angstroms).
                        e.g. tol = 2.0 angstroms
        packmol_path (str): path to the packmol executable. Defaults to "packmol"
        clean (bool): whether the intermedite files generated are deleted. Defaults to True.
    """

    required_params = ["composition", "box_scale"]
    optional_params = ["packmol_path", "clean", "tol"]

    def run_task(self, fw_spec):
        glass = AmorphousMaker(self.get("composition"), self.get("box_scale"), self.get("tol", 2.0),
                               packmol_path=self.get("packmol_path", "packmol"),
                               clean=self.get("clean", True))
        structure = glass.random_packed_structure
        return FWAction(stored_data=structure)

@explicit_serialize
class VaspMdToDbTask(FireTaskBase):
    pass

@explicit_serialize
class VaspMdToDiffusion(FireTaskBase):
    pass

@explicit_serialize
class VaspMdToStructuralAnalysis(FireTaskBase):
    pass