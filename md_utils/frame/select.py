from typing import List
from pyxmolpp2 import AtomSelection


def get_sec_str_residue_ids(path_to_pdb: str,
                            molname: str
                            ) -> List[int]:
    """
    to predict the secondary structure residues use DSSP https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
    :param path_to_pdb:
    :param molname:
    :return: list of residue id, which belongs to secondary-structure region
             helix - H (Alpha helix (4-12), G (3-10 helix) and I (Pi helix)
             or
             strand E (Strand) and B (Isolated beta-bridge residue),
    """
    ...


def select_sec_str_ca_atoms(path_to_pdb: str,
                            molname: str,
                            ) -> AtomSelection:
    """
    :param: frame: frame in MD trajectory
    :param path_to_pdb:
    :param molname: name of molecules.
                    For example, reference structure is complex so that it contains two molecules of interacting partner
    :return: CA atoms belongs to secondary structured regions
    """
    sec_str_rids = get_sec_str_residue_ids(path_to_pdb=path_to_pdb,
                                           molname=molname)
    ...


if __name__ == "__main__":
    path_to_pdb = "7jzu.pdb"

    ss_ca_atoms = select_sec_str_ca_atoms(path_to_pdb, molname="A")
