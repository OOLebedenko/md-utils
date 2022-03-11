from pyxmolpp2 import Molecule, Residue, AngleValue, ResidueSelection
from scipy.spatial import cKDTree # to compute distacnes
from typing import List
import numpy as np


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
    a_index, b_index = np.where(dist > 0)
    a_index, b_index = set(a_index), set(b_index)

    # collect corresponding residues into sets
    res_a, res_b = set(), set()

    a_atoms, b_atoms = partner_A.atoms, partner_B.atoms
    for index in a_index:
      res_a.add(a_atoms[index].residue.index)

    for index in b_index:
      res_b.add(b_atoms[index].residue.index)

    return [list(res_a), list(res_b)]
    

if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, mName

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    [interface_residues_partner_A, interface_residues_partner_B] = extract_residues_on_interface(
                                                                    partner_A=frame.molecules.filter(mName == "A"),
                                                                    partner_B=frame.molecules.filter(mName == "B"),
                                                                    cutoff=5)
