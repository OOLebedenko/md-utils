from pyxmolpp2 import Molecule, ResidueSelection
from typing import List


def extract_residues_on_interface(partner_A: Molecule,
                                  partner_B: Molecule,
                                  cutoff: float) -> List[ResidueSelection]:
    """
    Extract residues on interface of complex between partner_A and partner_B
    :param partner_A:
    :param partner_b:
    :param cutoff: distance cutoff in angstroms
    :return: list of [residues of partner_A, residues of partner_B] - residues on interface
                  of interaction between partner_A and partner_B within the distance cutoff
    """
    ...


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, mName

    path_to_pdb = "7jzu.pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    [interface_residues_partner_A, interface_residues_partner_B] = extract_residues_on_interface(
        partner_A=frame.molecules.filter(mName == "A"),
        partner_B=frame.molecules.filter(mName == "B"),
        cutoff=5)
