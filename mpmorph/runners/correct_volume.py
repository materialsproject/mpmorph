import numpy as np
from pymatgen.io.vasp import Poscar

class RescaleVolume(object):
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
    def __init__(self, structure, initial_pressure=0.0, initial_temperature=1000.0,
                         target_pressure=0.0, target_temperature=1000.0,
                         alpha=10e-5, beta=10e-7, poscar=None):
        self.structure = structure
        self.initial_pressure = initial_pressure  # in bars
        self.initial_temperature = initial_temperature  # in K
        self.target_pressure = target_pressure  # in bars
        self.target_temperature = target_temperature  # in K
        self.alpha = alpha  #  /K
        self.beta = beta  # /bar
        self.poscar = poscar

    def rescale_structure_volume(self, v2_v1):
        """
        Scales the volume of a structure by the given factor v2_v1.
        Returns the same structure with the new volume.
        """
        return self.structure.scale_lattice(self.structure.volume * v2_v1)

    def by_thermo(self, scale='pressure'):
        """
        Scales the volume of structure using thermodynamic functions
        Args:
            scale (str): thermodynamic function used to scale; 'temperature' or 'pressure'.
        Returns:

        """
        if scale == 'pressure':
            v2_v1 = np.exp(-self.beta * (self.target_pressure - self.initial_pressure))
            self.rescale_structure_volume(v2_v1)
            self.initial_pressure = self.target_pressure
        elif scale == 'temperature':
            v2_v1 = np.exp(self.alpha * (self.target_temperature - self.initial_temperature))
            self.rescale_structure_volume(v2_v1)
            self.initial_temperature = self.target_temperature
        else:
            raise ValueError("scale function must be specified as temperature or pressure.")

        return self.structure

    def by_eos(self, p_v, eos='polynomial'):
        """
        Args:
            p_v (numpy array): an array of pressure-volume pairs; e.g. p_v = [[p1,v1],[p2,v2],...]
            eos (str): type of equation of state to fit to. Currently only polynomial is supported.
                - 'polynomial': if p_v has two elements, fits line, otherwise; quadratic polynomial
                - 'bm': Birch-Murnaghan EOS. Not implemented yet.
        Returns:

        """
        if eos=='polynomial':
            v1 = self.structure.volume
            v2_v1 = poly_rescale(p_v, target_pressure=self.target_pressure)/v1
            self.rescale_structure_volume(v2_v1)
            self.initial_pressure = self.target_pressure
        elif eos=='bm':

            raise ValueError("BM EOS is not implemented yet")

        return self.structure



    @classmethod
    def of_poscar(cls, poscar_path, initial_pressure=0.0, initial_temperature=1000.0,
                         target_pressure=0.0, target_temperature=1000.0,
                         alpha=10e-5, beta=10e-7):
        """
        Convenience constructor that accepts a poscar file as input

        """
        poscar = Poscar.from_file(poscar_path)
        return cls(poscar.structure, initial_pressure=initial_pressure, initial_temperature=initial_temperature,
                         target_pressure=target_pressure, target_temperature=target_temperature,
                         alpha=alpha, beta=beta, poscar=poscar)

def poly_rescale(p_v, target_pressure=0.0):
    """
    :param p_v: a numpy array of pressure-volume pairs p_v = [[p1,v1],[p2,v2],...]
    if p_v has 2 elements, fit a line and return volume for zero pressure
    if p_v has 3 or more elements, fit a second order polynomial and return
    :return: the volume at target_pressure
    """
    if len(p_v)<2:
        raise ValueError("At least 2 p-v points required to estimate")
    if len(p_v)==2:
        eqs = np.poly1d(np.polyfit(p_v[:,0], p_v[:,1],1))
    else:
        eqs = np.poly1d(np.polyfit(p_v[:,0], p_v[:,1],2))
    #Return volume at zero pressure:
    return eqs(target_pressure)

def BirchMurnaghanEOS(p_v, target_pressure=0.0):



    return eqs(target_pressure)

def energy_bm(V,E0,B0,V0,B0p):
    """
    Args:
        V: volume
        E0: equilibrium energy
        B0: bulk modulus
        V0: equilibrium volume
        B0p: pressure derivative of B0

    Returns:
        Energy of Birch-Murnaghan EOS at V with given parameters E0, B0, V0 and B0p

    """
    n = (V/V0)**(1./3)
    return E0 + 9.*B0*V0/16*(n**2-1)**2 *(6.+B0p*(n**2-1)-4*n**2)

def pressure_bm(V,E0,B0,V0,B0p):
    """
    Args:
        V: volume
        E0: equilibrium energy
        B0: bulk modulus
        V0: equilibrium volume
        B0p: pressure derivative of B0

    Returns:
        Pressure of Birch-Murnaghan EOS at V with given parameters E0, B0, V0 and B0p

    """
    n = (V/V0)**(1./3)
    return


if __name__ == '__main__':
    p_v = np.array([[-0.1,100.0],[0.5,90],[2.0,70]])
