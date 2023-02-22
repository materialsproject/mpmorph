import itertools

from monty.io import zopen


class Xdatcar_Writer:
    def write_xdatcar(self, filename, **kwargs):
        """
        Writes Xdatcar to a file. The supported kwargs are the same as those for
        the Xdatcar_from_structs.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_string_from_struct(**kwargs))

    def get_string_from_struct(
        self, structures, system="unknown system", significant_figures=6
    ):
        format_str = "{{:.{0}f}}".format(significant_figures)

        for si, structure in enumerate(structures):
            lines = [system, "1.0", str(structure.lattice)]
            lines.append(" ".join(self.get_site_symbols(structure)))
            lines.append(" ".join([str(x) for x in self.get_natoms(structure)]))

            lines.append("Direct configuration=     " + str(si + 1))
            for i, site in enumerate(structure):
                coords = site.frac_coords
                line = " ".join([format_str.format(c) for c in coords])
                line += " " + site.species_string
                lines.append(line)

        return "\n".join(lines) + "\n"

    def get_site_symbols(self, structure):
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in
        vasp 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in structure]
        return [a[0] for a in itertools.groupby(syms)]

    def get_natoms(self, structure):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """
        syms = [site.specie.symbol for site in structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]


class Xdatcar_Writer_Trajectory:
    def __init__(self, trajectory):
        self.trajectory = trajectory

    def write_xdatcar(self, filename, **kwargs):
        """
        Writes Xdatcar to a file. The supported kwargs are the same as those for
        the Xdatcar_from_structs.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_string(**kwargs))

    def get_string(self, system="unknown system", significant_figures=6):
        lines = [system, "1.0"]
        for direction in self.trajectory.lattice:
            lines.append(" ".join([str(i) for i in direction]))
        lines.append(" ".join(self.get_site_symbols()))
        lines.append(" ".join([str(x) for x in self.get_natoms()]))

        format_str = "{{:.{0}f}}".format(significant_figures)
        positions = self.trajectory.frac_coords
        #  positions = np.add(self.trajectory[0].frac_coords, self.trajectory.displacements)
        atoms = [site.specie.symbol for site in self.trajectory[0]]

        for si, position_array in enumerate(positions):
            lines.append("Direct configuration=     " + str(si + 1))
            for i, coords in enumerate(position_array):
                line = " ".join([format_str.format(c) for c in coords])
                line += " " + atoms[i]
                lines.append(line)

        return "\n".join(lines) + "\n"

    def get_site_symbols(self):
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in
        vasp 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in self.trajectory[0]]
        return [a[0] for a in itertools.groupby(syms)]

    def get_natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """
        syms = [site.specie.symbol for site in self.trajectory[0]]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]
