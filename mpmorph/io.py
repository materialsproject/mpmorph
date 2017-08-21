from monty.io import zopen
import itertools

class Xdatcar_Writer():

    def write_xdatcar(self, filename, **kwargs):
        """
        Writes Xdatcar to a file. The supported kwargs are the same as those for
        the Xdatcar_from_structs.get_string method and are passed through directly.
        """
        with zopen(filename, "wt") as f:
            f.write(self.get_string_from_struct(**kwargs))

    def get_string_from_struct(self, structures, system="unknown system", significant_figures=6):
        lines = [system, "1.0", str(structures[0].lattice)]
        lines.append(" ".join(self.get_site_symbols(structures[0])))
        lines.append(" ".join([str(x) for x in self.get_natoms(structures[0])]))

        format_str = "{{:.{0}f}}".format(significant_figures)
        for (si, structure) in enumerate(structures):
            lines.append("Direct configuration=     " + str(si + 1))
            for (i, site) in enumerate(structure):
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