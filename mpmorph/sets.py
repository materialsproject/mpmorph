from pymatgen.io.vasp.sets import MPRelaxSet

class FrozenPhononSet(MPRelaxSet):
    """
    """

    def __init__(self, structure, incar_changes={}, **kwargs):
        super(FrozenPhononSet, self).__init__(structure, **kwargs)
        self.incar_changes = incar_changes
        self.kwargs = kwargs

        incar_updates = {"IBRION": 5, "NFREE": 4, "NSW": 100, "EDIFF": 1e-6,
                         "ISMEAR": 2, "SIGMA": 0.2, "ALGO": "FAST", "MAXMIN": 60}
        self._config_dict["INCAR"].update(incar_updates)
        self._config_dict["INCAR"].update(incar_changes)