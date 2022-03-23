from pyxmolpp2 import Molecule, aId
from scipy.spatial import cKDTree  # to compute distacnes


def extract_residues_on_interface(partner_A: Molecule,
                                  partner_B: Molecule,
                                  cutoff: float):
    """
    Extract residues on interface of complex between partner_A and partner_B
    :param partner_A:
    :param partner_B:
    :param cutoff: distance cutoff in angstroms
    :return: list of [residues of partner_A, residues of partner_B] - residues on interface
                  of interaction between partner_A and partner_B within the distance cutoff
    """
    # get coordinates and construct binary trees
    a_tree = cKDTree(partner_A.atoms.coords.values)
    b_tree = cKDTree(partner_B.atoms.coords.values)

    # compute distances with a cut-off
    dist = a_tree.sparse_distance_matrix(b_tree, max_distance=cutoff)
    dist = dist.toarray()

    # get atom indexes in cut-off range
    a_index, b_index = dist.nonzero()

    # get corresponding residues
    first_a_id = partner_A.atoms[0].id  # to correct indexes from matrix
    # atom_index + first_atom_index + (first_atom_id - first_atom_index) = atom_index + first_atom_id
    resdue_selection_A = partner_A.atoms.filter(aId.is_in(set(a_index + first_a_id))).residues

    first_b_id = partner_B.atoms[0].id  # to correct indexes from matrix
    resdue_selection_B = partner_B.atoms.filter(aId.is_in(set(b_index + first_b_id))).residues

    return [resdue_selection_A, resdue_selection_B]


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, mName

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    [interface_residues_partner_A, interface_residues_partner_B] = extract_residues_on_interface(
        partner_A=frame.molecules.filter(mName == "A"),
        partner_B=frame.molecules.filter(mName == "B"),
        cutoff=5)
