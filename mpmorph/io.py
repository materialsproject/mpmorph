from pymatgen.io.vasp.outputs import Xdatcar
from monty.io import zopen

class Xdatcar_Structures(Xdatcar):

    def write_xdatcar(self, **kwargs):
        """
        Writes Xdatcar to a file. The supported kwargs are the same as those for
        the Xdatcar_from_structs.get_string method and are passed through directly.
        """
        with zopen(self.filename, "wt") as f:
            f.write(self.get_string_from_struct(**kwargs))

    def get_string_from_struct(self, structures, system="unknown system", significant_figures=6):
        lines = [system, "1.0", str(structures[0].lattice)]
        lines.append(" ".join(self.site_symbols))
        lines.append(" ".join([str(x) for x in self.natoms]))

        format_str = "{{:.{0}f}}".format(significant_figures)
        for (si, structure) in enumerate(structures):
            lines.append("Direct configuration=     " + str(si + 1))
            for (i, site) in enumerate(structure):
                coords = site.frac_coords
                line = " ".join([format_str.format(c) for c in coords])
                line += " " + site.species_string
                lines.append(line)

        return "\n".join(lines) + "\n"