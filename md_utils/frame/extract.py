from pyxmolpp2 import Molecule, Residue, AngleValue, ResidueSelection
from typing import List


def extract_dihedral_angle_phi(residue: Residue) -> AngleValue:
    """
    see https://proteopedia.org/wiki/index.php/Dihedral/Dihedral_angles_in_proteins
    :param residue:
    :return: dihedral angle phi
    """
    ...


def extract_dihedral_angle_psi(residue: Residue) -> AngleValue:
    """
    see https://proteopedia.org/wiki/index.php/Dihedral/Dihedral_angles_in_proteins
    :param residue:
    :return: dihedral angle psi
    """
    ...


def extract_residues_on_interface(partner_A: Molecule,
                                  partner_B: Molecule,
                                  cutoff: float) -> List[ResidueSelection]:
    """
    Extract residues on interface of complex between partner_A and partner_B
    :param partner_A:
    :param partner_B:
    :param cutoff: distance cutoff in angstroms
    :return: list of [residues of partner_A, residues of partner_B] - residues on interface
                  of interaction between partner_A and partner_B within the distance cutoff
    """
    ...


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, mName

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    phi_angles = []
    psi_angles = []
    omega_angles = []

    for residue in frame.residues:
        phi_angles.append(extract_dihedral_angle_phi(residue))
        psi_angles.append(extract_dihedral_angle_psi(residue))

    [interface_residues_partner_A, interface_residues_partner_B] = extract_residues_on_interface(
                                                                    partner_A=frame.molecules.filter(mName == "A"),
                                                                    partner_B=frame.molecules.filter(mName == "B"),
                                                                    cutoff=5)
