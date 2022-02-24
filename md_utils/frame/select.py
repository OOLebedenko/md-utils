from typing import Dict
from pyxmolpp2 import Frame, AtomSelection


def get_ss_residies(path_to_pdb: str) -> Dict[str, str]:
    """
    to predict the secondary structure residues use DSSP https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
    :param path_to_pdb
    :return: ss_rids_rnames_dict: key - residue id, value - type of secondary structure in DSSP notation
                         helix - H (Alpha helix (4-12), G (3-10 helix) and I (Pi helix),
                         strand E (Strand) and B (Isolated beta-bridge residue),
                         loop - S (Bend), T (Turn)
                         None (-)
    """
    ...


def select_ss_ca_atoms(frame: Frame,
                       ) -> AtomSelection:
    """
    :param path_to_pdb:
    :param molnames:
    :return: CA atoms belongs to secondary structured regions
    """
    ...


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    ss_ca_atoms = select_ss_ca_atoms(frame)
