from fireworks import explicit_serialize, FireTaskBase, FWAction, Firework

__author__ = "Eric Sivonxay <esivonxay@lbl.gov>"

@explicit_serialize
class DiffusionFW(FireTaskBase):
