import numpy as np
from pymatgen.io.vasp import Poscar

class RescaleVolume(object):
    """
    Class for adjusting the volume of an input simulation box based on conditions.
    """

    def __init__(self, structure, initial_pressure=0.0, initial_temperature=1000.0,
                         target_pressure=0.0, target_temperature=1000.0,
                         alpha=10e-5, beta=10e-7, poscar=None):
        """
        Args:
            structure:
            initial_pressure: in bars
            initial_temperature: in Kelvins
            target_pressure: in bars
            target_temperature: n Kelvins
            alpha:
            beta:
            poscar:

        Returns:

        """
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
        Scales the volume of structure using thermodynamic functions, which basically give linear
        Equations of State. For more advanced EOS, one should use the by_EOS method.
        Args:
            scale (str): thermodynamic function used to scale; 'temperature' or 'pressure'.
        Returns:
            rescaled structure

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

    def by_EOS(self, p_v, eos='polynomial'):
        """
        Args:
            p_v (numpy array): an array of pressure-volume pairs; e.g. p_v = [[p1,v1],[p2,v2],...]
            eos (str): type of equation of state to fit to. Currently only polynomial is supported.
                - 'polynomial': if p_v has two elements, fits line, otherwise; quadratic polynomial
                - 'bm': Birch-Murnaghan EOS. Not implemented yet.
        Returns:

        """
        v1 = self.structure.volume
        if eos=='polynomial':
            v2_v1 = poly_rescale(p_v, target_pressure=self.target_pressure)/v1
            self.rescale_structure_volume(v2_v1)
            self.initial_pressure = self.target_pressure
        elif eos=='Murnaghan':
            raise ValueError("not implemented yet")
        elif eos=='BirchMurnaghan':
            v2_v1 = fit_BirchMurnaghanPV_EOS(p_v, target_pressure=self.target_pressure)/v1
            self.rescale_structure_volume(v2_v1)
            self.initial_pressure = self.target_pressure
        else:
            raise ValueError("Unknown EOS. Volume not rescaled.")
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
    return eqs(target_pressure)


def BirchMurnaghanPV_EOS(V,params):
    """
    Args:
        V: volume
        params: tuple of B0,V0,B0p
    Returns:
        Pressure of Birch-Murnaghan EOS at V with given parameters E0, B0, V0 and B0p
    """
    V0, B0, B0p = params[0], params[1], params[2]
    n = (V0/V)**(1./3) # not this definition is different from the Energy EOS
    p = 3./2*B0*(n**7-n**5)*(1.+3./4*(B0p-4)*(n**2-1))
    return p

def fit_BirchMurnaghanPV_EOS(p_v):
    # Borrows somewhat from pymatgen/io/abinitio/EOS
    # Initial guesses for the parameters
    from scipy.optimize import leastsq
    eqs = np.polyfit(p_v[:,1], p_v[:,0], 2)
    V0 = np.mean(p_v[:,1]) # still use mean to ensure we are at reasonable volumes
    B0 = -1*(2*eqs[0]*V0**2+eqs[1]*V0)
    B0p = -1*(2*eqs[0]*V0**2+eqs[1])
    initial_params = (V0,B0,B0p)
    print initial_params
    Error=lambda params,x,y: BirchMurnaghanPV_EOS(x,params) - y
    found_params, check = leastsq(Error,initial_params,args=(p_v[:,1],p_v[:,0]))
    print check
    if check not in [1,2,3,4]:
        raise ValueError("fitting not converged")
    else:
        return found_params

def BirchMurnaghan_rescale(p_v, target_pressure=0):
    """
    Calls fit_BirchMurnaghanPV_EOS to find params of EOS and returns V corresponding to target_pressure
    Args:
        p_v:
        target_pressure:

    Returns:

    """
    params = fit_BirchMurnaghanPV_EOS(p_v)
    if target_pressure==0:
        return params[0]
    else:
        # TODO: find volume corresponding to this target_pressure
        pass

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    p_v = np.array([[-0.1,100.0],[0.5,90],[3.0,70]])
    p = p_v[:,0]
    v = p_v[:,1]
    params= fit_BirchMurnaghanPV_EOS(p_v)
    print params
    xx = np.linspace(v.min(),v.max(),50)
    yy = BirchMurnaghanPV_EOS(xx,params)
    plt.plot(xx,yy,'r-',v,p,'bo')
    plt.show()
    print BirchMurnaghan_rescale(p_v)




