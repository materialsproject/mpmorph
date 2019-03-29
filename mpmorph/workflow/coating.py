import numpy as np
from pymatgen import Lattice
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, get_mi_vec
from pymatgen.core.surface import Slab
from itertools import product
from pymatgen.util.coord import pbc_shortest_vectors
import warnings

from atomate.vasp.workflows.base.adsorption import MPSurfaceSet
from mpmorph.fireworks import powerups
from fireworks import Workflow
import uuid
from mpmorph.fireworks.core import MDFW, OptimizeFW, StaticFW
from custodian.vasp.handlers import VaspErrorHandler, MeshSymmetryErrorHandler, UnconvergedErrorHandler, \
    PotimErrorHandler, FrozenJobErrorHandler, NonConvergingErrorHandler, PositiveEnergyErrorHandler, \
    StdErrHandler

handler_group = [VaspErrorHandler(), MeshSymmetryErrorHandler(), UnconvergedErrorHandler(),
                 NonConvergingErrorHandler(), PotimErrorHandler(),
                 PositiveEnergyErrorHandler(), FrozenJobErrorHandler(), StdErrHandler()]


def get_film_wf(bulk_substrate, bulk_film, slab_substrate, slab_film, descriptor, md_prerelax=True, skip_bulk=False,
                wf_descriptor="", offset=None, **kwargs):
    """
    This workflow currently does not allow for calculation of interfacial energy.
    """
    fw_list = []

    tag_id = uuid.uuid4()
    # Check if lattices are equal. If not, strain them to match
    latt_1 = slab_substrate.lattice.matrix.copy()
    latt_1[2, :] = [0, 0, 1]
    latt_2 = bulk_substrate.lattice.matrix.copy()
    latt_2[2, :] = [0, 0, 1]
    if not Lattice(latt_1) == Lattice(latt_2):
        # Calculate lattice strained to match:
        matched_slab_substrate, matched_slab_film = strain_slabs(slab_substrate, slab_film)
    else:
        matched_slab_substrate = slab_substrate
        matched_slab_film = slab_film

    # Ensure substrate has positive c-direction:
    if matched_slab_substrate.lattice.matrix[2, 2] < 0:
        latt = matched_slab_substrate.lattice.matrix.copy()
        latt[2, 2] *= -1
        new_struct = matched_slab_substrate.copy()
        new_struct.lattice = Lattice(latt)
        matched_slab_substrate = new_struct

    # Ensure film has positive c-direction:
    if matched_slab_film.lattice.matrix[2, 2] < 0:
        latt = matched_slab_film.lattice.matrix.copy()
        latt[2, 2] *= -1
        new_struct = matched_slab_film.copy()
        new_struct.lattice = Lattice(latt)
        matched_slab_film = new_struct

    if not skip_bulk:
        # Add optimize + static calculations for the bulk structures
        sub_static_fw = StaticFW(structure=bulk_substrate, name=descriptor + '-bulk_sub_static-' + str(tag_id),
                                 vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<')
        fw_list.extend([sub_static_fw])

        film_static_fw = StaticFW(structure=bulk_film, name=descriptor + '-bulk_film_static-' + str(tag_id),
                                  vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<')
        fw_list.extend([film_static_fw])

    # Add static calculations for the substrate slab
    # sub_slab_opt_fw = OptimizeFW(structure=matched_slab_substrate)
    slab_substrate_with_vacuum = add_vaccuum(slab_substrate)
    sub_slab_static_fw = StaticFW(structure=slab_substrate_with_vacuum, name=descriptor + '-slab_sub_static-' + str(tag_id),
                                  vasp_cmd='>>vasp_cmd<<', db_file='>>db_file<<')
    fw_list.extend([sub_slab_static_fw])

    # Add optimize + static calculations for the film slab
    slab_film_with_vacuum = add_vaccuum(slab_film)
    matched_slab_film_sd = assign_selective_dynamics_by_height(slab_film_with_vacuum, 3)
    # matched_slab_film_sd = fix_one_atom(slab_film_with_vacuum)

    if md_prerelax:
        run_args = {"md_params": {"start_temp": 500, "end_temp": 500, "nsteps": 200},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"},
                    "optional_fw_params": {"override_default_vasp_params": {}, "spec": {}},
                    "label": "md_prerelax_"}

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False, 'AMIN': 0.01, 'ALGO': "Normal"}})
        fw = MDFW(structure=matched_slab_film_sd, name=descriptor + "-" + run_args["label"] + "0" + "-" + str(tag_id),
                  previous_structure=False, insert_db=True, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"], parents=[])
        fw_list.append(fw)

        for i in range(1, 5):
            fw = MDFW(structure=matched_slab_film_sd,
                      name=descriptor + "-" + run_args["label"] + str(i) + "-" + str(tag_id),
                      previous_structure=True, insert_db=True, **run_args["md_params"],
                      **run_args["run_specs"], **run_args["optional_fw_params"], parents=[fw_list[-1]])
            fw_list.append(fw)

    parents = [fw_list[-1]] if md_prerelax else []
    run_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "spec": {}},
                "optional_fw_params": {"override_default_vasp_params": {}}}
    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'LVTOT': True, 'LDIPOL': True, 'IDIPOL': 3, 'AMIN': 0.01,
                                 'EDIFF': 0.000001 * matched_slab_film_sd.num_sites, 'ALGO': 'All', 'KPAR': 4,
                                 'ICHARG': 0, 'POTIM': 0.25, 'NSW': 200, 'IBRION': 2}})
    vasp_input_set = MPSurfaceSet(matched_slab_film_sd, auto_dipole=True,
                                  **run_args['optional_fw_params']["override_default_vasp_params"])
    run_args["run_specs"]["vasp_input_set"] = vasp_input_set
    film_slab_opt_fw = OptimizeFW(structure=matched_slab_film_sd,
                                  name=descriptor + '-slab_film_optimize-' + str(tag_id),
                                  previous_structure=md_prerelax,
                                  parents=parents, **run_args["run_specs"], **run_args["optional_fw_params"],
                                  max_force_threshold=None, handler_group=handler_group, job_type='normal')

    run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "spec": {}},
                "optional_fw_params": {"override_default_vasp_params": {}}}
    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'LVTOT': True, "LDIPOL": True, "IDIPOL": 3, 'AMIN': 0.01,
                                 'EDIFF': 0.000001 * matched_slab_film_sd.num_sites, 'ALGO': "All",
                                 'ICHARG': 0, 'KPAR': 4}})
    film_slab_static_fw = StaticFW(structure=matched_slab_film_sd, name=descriptor + '-slab_film_static-' + str(tag_id),
                                   parents=[film_slab_opt_fw],
                                   **run_args["run_specs"], **run_args["optional_fw_params"], previous_structure=True)
    fw_list.extend([film_slab_opt_fw, film_slab_static_fw])

    # Create interface structure and its corresponding optimization and static fireworks
    interface_structure = make_interface(matched_slab_substrate, matched_slab_film, offset)

    if md_prerelax:
        run_args = {"md_params": {"start_temp": 500, "end_temp": 500, "nsteps": 200},
                    "run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"},
                    "optional_fw_params": {"override_default_vasp_params": {}, "spec": {}},
                    "label": "md_prerelax_"}

        run_args["optional_fw_params"]["override_default_vasp_params"].update(
            {'user_incar_settings': {'ISIF': 1, 'LWAVE': False, 'AMIN': 0.01, 'ALGO': "Normal"}})
        fw = MDFW(structure=interface_structure, name=descriptor + "-" + run_args["label"] + "0" + "-" + str(tag_id),
                  previous_structure=False, insert_db=True, **run_args["md_params"],
                  **run_args["run_specs"], **run_args["optional_fw_params"], parents=[])
        fw_list.append(fw)

        for i in range(1, 5):
            fw = MDFW(structure=interface_structure,
                      name=descriptor + "-" + run_args["label"] + str(i) + "-" + str(tag_id),
                      previous_structure=True, insert_db=True, **run_args["md_params"],
                      **run_args["run_specs"], **run_args["optional_fw_params"], parents=[fw_list[-1]])
            fw_list.append(fw)

    parents = [fw_list[-1]] if md_prerelax else []
    run_args = {"run_specs": {"vasp_input_set": None, "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "spec": {}},
                "optional_fw_params": {"override_default_vasp_params": {}}}
    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'LVTOT': True, 'LDIPOL': True, 'IDIPOL': 3, 'AMIN': 0.01,
                                 'EDIFF': 0.000001 * interface_structure.num_sites, 'ALGO': 'All', 'KPAR': 4,
                                 'ICHARG': 0, 'POTIM': 0.25, 'NSW': 200, 'IBRION': 2}})
    vasp_input_set = MPSurfaceSet(interface_structure, auto_dipole=True, \
                                  **run_args['optional_fw_params']["override_default_vasp_params"])
    run_args["run_specs"]["vasp_input_set"] = vasp_input_set
    interface_opt_fw = OptimizeFW(structure=interface_structure, name=descriptor + '-interface_optimize-' + str(tag_id),
                                  previous_structure=md_prerelax,
                                  parents=parents, **run_args["run_specs"], **run_args["optional_fw_params"],
                                  max_force_threshold=None, handler_group=handler_group, job_type='normal')

    run_args = {"run_specs": {"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<",
                              "spec": {}},
                "optional_fw_params": {"override_default_vasp_params": {}}}
    run_args["optional_fw_params"]["override_default_vasp_params"].update(
        {'user_incar_settings': {'LVTOT': True, "LDIPOL": True, "IDIPOL": 3, 'AMIN': 0.01,
                                 'EDIFF': 0.000001 * interface_structure.num_sites, 'ALGO': "All",
                                 'ICHARG': 0, 'KPAR': 4}})
    interface_static_fw = StaticFW(structure=interface_structure, name=descriptor + '-interface_static-' + str(tag_id),
                                   parents=[interface_opt_fw],
                                   **run_args["run_specs"], **run_args["optional_fw_params"], previous_structure=True)
    fw_list.extend([interface_opt_fw, interface_static_fw])

    return Workflow(fw_list, name=wf_descriptor + '-coating_wf'), interface_structure


def make_interface(matched_slab_substrate, matched_slab_film, offset=None):
    matched_slab_substrate_sd = assign_selective_dynamics(matched_slab_substrate, sd_options='other')
    matched_slab_film_sd = assign_selective_dynamics(matched_slab_film, sd_options='all_true')

    if offset is None:
        offset = get_offset(matched_slab_substrate_sd, matched_slab_film_sd)
    print(offset)
    _structure = combine_slabs(matched_slab_substrate_sd, matched_slab_film_sd, *offset)
    orthogonal_structure = _structure.get_orthogonal_c_slab()
    orthogonal_structure.sort()

    if not orthogonal_structure.is_valid(tol=1):
        warnings.warn("Check generated structure, it may contain atoms too closely placed")

    return orthogonal_structure


def fix_one_atom(structure):
    layers = np.unique(structure.cart_coords[:, 2])
    index = np.where(structure.cart_coords[:, 2] == layers[int(len(layers) / 2)])[0][0]
    sd_list = [[True, True, True]] * len(structure.sites)
    sd_list[index] = [False, False, False]
    new_sp = structure.site_properties
    new_sp['selective_dynamics'] = sd_list
    return structure.copy(site_properties=new_sp)


def assign_selective_dynamics_by_height(structure, subsurface_thickness=3):
    layers = np.unique(np.round(structure.cart_coords[:, 2], decimals=3))

    upper_cutoff = max(layers) - subsurface_thickness
    lower_cutoff = min(layers) + subsurface_thickness
    sd_list = []
    for site in structure.sites:
        if site.coords[2] < upper_cutoff and site.coords[2] > lower_cutoff:
            sd_list.append([False, False, False])
        else:
            sd_list.append([True, True, True])
    new_sp = structure.site_properties
    new_sp['selective_dynamics'] = sd_list
    return structure.copy(site_properties=new_sp)


def assign_selective_dynamics(slab, sd_options='all_true'):
    vec = get_mi_vec(slab)
    mi_vec = vec * np.sign(vec[2])
    asf = AdsorbateSiteFinder(slab, mi_vec=mi_vec.tolist())
    if sd_options == 'vertical_only':
        # Only allow vertical movement (for film_subsurface)
        sd_list = [[False, False, True] if site.properties['surface_properties'] == 'subsurface'
                   else [True, True, True] for site in asf.slab.sites]
    elif sd_options == 'all_true':
        sd_list = [[True, True, True] for site in asf.slab.sites]
    elif sd_options == 'all_false':
        sd_list = [[False, False, False] for site in asf.slab.sites]
    else:
        sd_list = [[False, False, False] if site.properties['surface_properties'] == 'subsurface'
                   else [True, True, True] for site in asf.slab.sites]

    new_sp = slab.site_properties
    new_sp['selective_dynamics'] = sd_list
    new_sp['surface_properties'] = asf.slab.site_properties['surface_properties']
    return slab.copy(site_properties=new_sp)


def add_vaccuum(structure, vacuum=20):
    cell_z = np.abs(np.dot([0, 0, 1], structure.lattice.matrix[2, :]))
    multiplier = 1 + vacuum / cell_z

    struct_lattice = structure.lattice.matrix.copy()
    struct_lattice[2, :] *= multiplier

    new_structure = Slab(lattice=Lattice(struct_lattice), species=structure.species,
                         coords=structure.cart_coords,
                         miller_index=structure.miller_index,
                         oriented_unit_cell=structure.oriented_unit_cell,
                         shift=structure.shift,
                         scale_factor=structure.scale_factor,
                         coords_are_cartesian=True, energy=structure.energy,
                         reorient_lattice=False, to_unit_cell=True,
                         site_properties=structure.site_properties)
    return new_structure


def strain_slabs(sub_slab, film_slab):
    sub_struct = sub_slab.copy()
    latt_1 = sub_struct.lattice.matrix.copy()
    film_struct = align_x(film_slab, get_ortho_axes(sub_struct)).copy()
    latt_2 = film_struct.lattice.matrix.copy()

    # Rotate film so its diagonal matches with the sub's diagonal
    diag_vec = np.add(latt_1[0, :], latt_1[1, :])
    sub_norm_diag_vec = diag_vec / np.linalg.norm(diag_vec)
    sub_b = np.cross(sub_norm_diag_vec, [0, 0, 1])
    sub_matrix = np.vstack([sub_norm_diag_vec, sub_b, [0, 0, 1]])

    diag_vec = np.add(latt_2[0, :], latt_2[1, :])
    film_norm_diag_vec = diag_vec / np.linalg.norm(diag_vec)
    film_b = np.cross(film_norm_diag_vec, [0, 0, 1])
    film_matrix = np.vstack([film_norm_diag_vec, film_b, [0, 0, 1]])

    rotation = np.dot(np.linalg.inv(film_matrix), sub_matrix)
    new_latt = Lattice(np.dot(film_struct.lattice.matrix, rotation))
    film_struct.lattice = new_latt

    # Average the two lattices (Should get equal strain?)
    mean_a = np.mean([film_struct.lattice.matrix[0, :], sub_struct.lattice.matrix[0, :]], axis=0)
    mean_b = np.mean([film_struct.lattice.matrix[1, :], sub_struct.lattice.matrix[1, :]], axis=0)
    new_latt = np.vstack([mean_a, mean_b, sub_struct.lattice.matrix[2, :]])
    sub_struct.lattice = Lattice(new_latt)
    new_latt = np.vstack([mean_a, mean_b, film_struct.lattice.matrix[2, :]])
    film_struct.lattice = Lattice(new_latt)

    return sub_struct, film_struct


def transf_mat(A, B):
    return np.dot(np.linalg.inv(A), B)


def third_vect(a, b):
    c = np.cross(a, b);
    return c / np.linalg.norm(c)


def get_ortho_axes(structure):
    sub_a = structure.lattice.matrix[0, :] / np.linalg.norm(structure.lattice.matrix[0, :])
    sub_c = third_vect(sub_a, structure.lattice.matrix[1, :])

    sub_b = third_vect(sub_c, sub_a)
    sub_b = sub_b / np.linalg.norm(sub_b)

    return np.vstack((sub_a, sub_b, sub_c))


def align_x(slab, orthogonal_basis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
    sub_ortho_axes = get_ortho_axes(slab)
    rotation = transf_mat(sub_ortho_axes, orthogonal_basis)
    new_sub_lattice = Lattice(np.dot(slab.lattice.matrix[0:3], rotation))
    slab.lattice = new_sub_lattice
    return slab


def get_offset(substrate, film):
    # Coarse search for z offset

    last_offset = (0.2, -0.5, -0.5)
    for (slab_offset, x_offset, y_offset) in product(np.arange(0.2, 3, 0.1), [-0.5, 0, 0.5], [-0.5, 0, 0.5]):
        valid = is_valid_offset(substrate, film, 20, slab_offset, x_offset, y_offset)
        #         print((slab_offset, x_offset, y_offset))
        if valid:
            last_offset = (slab_offset, x_offset, y_offset)
            break

    # Finer search for offset
    for (slab_offset, x_offset, y_offset) in product(np.arange(last_offset[0] - 0.4, last_offset[0], 0.05),
                                                     np.arange(-1, 1, 0.1), np.arange(-1, 1, 0.1)):
        valid = is_valid_offset(substrate, film, 20, slab_offset, x_offset, y_offset)
        if valid:
            last_offset = (slab_offset, x_offset, y_offset)
            break

    return last_offset


# def assign_selective_dynamics(slab, sd_options='all_true'):
#     asf = AdsorbateSiteFinder(slab)
#     if sd_options == 'vertical_only':
#         # Only allow vertical movement (for film_subsurface)
#         sd_list = [[False, False, True] if site.properties['surface_properties'] == 'subsurface'
#                    else [True, True, True] for site in asf.slab.sites]
#     elif sd_options == 'all_true':
#         sd_list = [[True, True, True] for site in asf.slab.sites]
#     elif sd_options == 'all_false':
#         sd_list = [[False, False, False] for site in asf.slab.sites]
#     else:
#         sd_list = [[False, False, False] if site.properties['surface_properties'] == 'subsurface'
#                    else [True, True, True] for site in asf.slab.sites]
#
#     new_sp = slab.site_properties
#     new_sp['selective_dynamics'] = sd_list
#     new_sp['surface_properties'] = asf.slab.site_properties['surface_properties']
#     return slab.copy(site_properties=new_sp)


def combine_slabs(substrate, film, slab_offset, x_offset, y_offset, vacuum=20, **kwargs):
    # strain film to match substrate
    new_latt = film.lattice.matrix.copy()
    new_latt[:2, :2] = substrate.lattice.matrix[:2, :2]
    film.lattice = Lattice(new_latt)

    combined_species = [*substrate.species, *film.species]
    if kwargs.get('cell_height'):
        height = kwargs.get('cell_height')
    else:
        added_height = vacuum + slab_offset + film.lattice.c
        height = added_height + substrate.lattice.matrix[2, 2]
    combined_lattice = substrate.lattice.matrix.copy()
    combined_lattice[2, :] *= height / substrate.lattice.matrix[2, 2]

    max_substrate = np.max(substrate.cart_coords[:, 2])
    min_substrate = np.min(film.cart_coords[:, 2])
    offset = max_substrate - min_substrate + slab_offset
    offset_film_coords = [np.add(coord, [x_offset, y_offset, offset]) for coord in film.cart_coords]
    combined_coords = [*substrate.cart_coords, *offset_film_coords]
    combined_site_properties = {}
    for key, item in substrate.site_properties.items():
        combined_site_properties[key] = [*substrate.site_properties[key], *film.site_properties[key]]

    combined_structure = Slab(lattice=Lattice(combined_lattice), species=combined_species,
                              coords=combined_coords,
                              miller_index=substrate.miller_index,
                              oriented_unit_cell=substrate,
                              shift=substrate.shift,
                              scale_factor=substrate.scale_factor,
                              coords_are_cartesian=True, energy=substrate.energy,
                              reorient_lattice=False, to_unit_cell=True,
                              site_properties=combined_site_properties)
    return combined_structure


def is_valid_offset(substrate, film, vacuum, slab_offset, x_offset, y_offset):
    # Calculate new lattice
    max_substrate = np.max(substrate.cart_coords[:, 2])
    min_substrate = np.min(film.cart_coords[:, 2])

    added_height = vacuum + slab_offset + film.lattice.c
    offset = max_substrate - min_substrate + slab_offset
    offset_film_coords = film.cart_coords + [x_offset, y_offset, offset]

    height = added_height + substrate.lattice.matrix[2, 2]
    combined_lattice = substrate.lattice.matrix.copy()
    combined_lattice[2, :] *= height / substrate.lattice.matrix[2, 2]
    #             combined_lattice = np.add(substrate.lattice.matrix, np.diag([0, 0, added_height]))
    combined_lattice = Lattice(combined_lattice)

    combined_species = [*substrate.species, *film.species]
    combined_coords = np.vstack((substrate.cart_coords, offset_film_coords))
    combined_species, combined_coords = zip(*sorted(zip(combined_species, combined_coords), key=lambda x: x[0]))
    frac_coords = np.dot(combined_coords, combined_lattice.inv_matrix)

    v, d2 = pbc_shortest_vectors(combined_lattice, frac_coords, frac_coords, return_d2=True)
    combined_distance_matrix = np.sqrt(d2)

    # Ensure that bonds aren't shorter than bond_lens in i) substrate or ii) film
    combined_bond_lens = get_bond_lens(combined_species, combined_distance_matrix)
    film_bond_lens = get_bond_lens_from_structure(film)
    substrate_bond_lens = get_bond_lens_from_structure(substrate)

    struct_is_valid = True
    tolerance = -0.005
    for pair, item in combined_bond_lens.items():
        _pair_in_film = pair in film_bond_lens.keys()
        _pair_in_substrate = pair in substrate_bond_lens.keys()
        if _pair_in_film and _pair_in_substrate:
            min_bond = min([film_bond_lens[pair][1], substrate_bond_lens[pair][1]])
            if (item[1] - min_bond) / min_bond < tolerance:
                struct_is_valid = False
                break
        else:
            if _pair_in_film:
                if (item[1] - film_bond_lens[pair][1]) / film_bond_lens[pair][1] < tolerance:
                    struct_is_valid = False
                    break

            if _pair_in_substrate:
                if (item[1] - substrate_bond_lens[pair][1]) / substrate_bond_lens[pair][1] < tolerance:
                    struct_is_valid = False
                    break

            if not (_pair_in_film or _pair_in_substrate):
                # TODO: Query MP for structure with that bond
                raise Exception('A bonding pair does not exist in either film or substrate')
    return struct_is_valid


def get_bond_lens(species, distance_matrix):
    species = np.array([str(item) for item in species])
    species_indices = [(el, np.where(species == el)[0]) for el in np.unique(species)]
    species_indices = dict(species_indices)

    bond_lens = {}
    for i, j in product(species_indices.keys(), species_indices.keys()):
        i_indices = (np.min(species_indices[i]), np.max(species_indices[i]))
        j_indices = (np.min(species_indices[j]), np.max(species_indices[j]))

        dm = distance_matrix[i_indices[0]:i_indices[1], j_indices[0]:j_indices[1]]
        bond_lens[(i, j)] = np.unique(np.round(dm, decimals=3))
    return bond_lens


def get_bond_lens_from_structure(structure):
    return get_bond_lens(structure.species, structure.distance_matrix)
