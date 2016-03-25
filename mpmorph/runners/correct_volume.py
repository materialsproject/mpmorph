import numpy as np
from pymatgen.io.vasp import Poscar

class CorrectVolume(object):
    """
    TODO : be able to get from Poscar, so that predictor corrector
    parameters can be returned correctly.

    Class for adjusting the volume of an input simulation box
    to zero pressure.
    Args:
        - Structure
        - initial pressure
        - initial temperature

    """
    def __init__(self, structure, pressure=0.0, temperature=1000.0, poscar = None):

        self.structure = structure
        self.initial_pressure = pressure  # in bars
        self.initial_temperature = temperature  # in K
        self.target_pressure = 0.0  # default target pressure is zero
        self.alpha = 10e-5  # default thermal exp coeff. is 10e-5 /K
        self.beta = 10e-7  # default compressibility is 10e-7 /bar
        self.poscar = poscar

    @classmethod
    def from_poscar(cls, poscar_path, pressure=0.0, temperature=1000.0):
        '''
        TODO :Get the Poscar object rather than
        Returns:

        '''
        poscar = Poscar.from_file(poscar_path)
        cls.poscar = poscar
        return cls(poscar.structure, pressure, temperature=1000.0, poscar = poscar)


    def rescale_vol(self, scale='pressure', target_pressure=None, target_temperature=None, alpha=None, beta=None):

        if scale == 'pressure':
            if target_pressure:
                self.target_pressure = target_pressure
            if beta:
                self.beta = beta
            v2_v1 = np.exp(-self.beta * (self.target_pressure - self.initial_pressure))
            self.initial_pressure = self.target_pressure  # Ensures multiple calls don't explode the structure

        elif scale == 'temperature':
            if target_temperature:
                self.target_temperature = target_temperature
            if alpha:
                self.alpha = alpha
            v2_v1 = np.exp(self.alpha * (self.target_temperature - self.initial_temperature))
            self.initial_temperature = self.target_temperature  # Ensures multiple calls don't explode the structure

        else:
            print("Unknown rescale option! Exiting.")
            return None

        self.structure.scale_lattice(self.structure.volume * v2_v1)

        return self.structure


def poly_rescale(p_v):
    """
    :param p_v: a numpy array of pressure-volume pairs p_v = [[p1,v1],[p2,v2],...]
    if p_v has 2 elements, fit a line and return volume for zero pressure
    if p_v has 3 or more elements, fit a second order polynomial and return
    :return: the vovlume fore zero pressure
    """
    if len(p_v)<2:
        raise ValueError("At least 2 p-v points required to estimate")
    if len(p_v)==2:
        eqs = np.poly1d(np.polyfit(p_v[:,0], p_v[:,1],1))
    else:
        eqs = np.poly1d(np.polyfit(p_v[:,0], p_v[:,1],2))
    #Return volume at zero pressure:
    return eqs(0)

if __name__ == '__main__':
    p_v = np.array([[-0.1,100.0],[0.5,90],[2.0,70]])
    print CorrectVolume.poly_rescale(p_v)