from typing import List
from pyxmolpp2 import AtomSelection, aName, rId, PdbFile, Frame
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import os


def get_sec_str_residue_ids(frame: Frame,
                            molname: str
                            ) -> List[int]:
    """
    to predict the secondary structure residues use DSSP https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
    :param frame:
    :param molname:
    :return: list of residue id, which belongs to secondary-structure region
             helix - H (Alpha helix (4-12), G (3-10 helix) and I (Pi helix)
             or
             strand E (Strand) and B (Isolated beta-bridge residue),
    """
    # predict secondary structure from input pdb
    frame.to_pdb("temp.pdb")
    dssp_tuple = dssp_dict_from_pdb_file("temp.pdb")
    dssp_dict = dssp_tuple[0]

    # get structured residues from target molecule
    sec_str_rids = []

    for key in dssp_dict.keys():
        if (key[0] == molname) and (dssp_dict[key][1] in ['H', 'G', 'I', 'E', 'B']):
            sec_str_rids.append(key[1][1])

    os.remove("temp.pdb")

    return sec_str_rids


def select_sec_str_ca_atoms(frame: Frame,
                            molname: str,
                            ) -> AtomSelection:
    """
    :param frame:
    :param molname: name of molecules.
                    For example, reference structure is complex so that it contains two molecules of interacting partner
    :return: CA atoms belongs to secondary structured regions
    """
    sec_str_rids = get_sec_str_residue_ids(frame=frame,
                                           molname=molname)
    # get frame
    frame = PdbFile(path_to_pdb).frames()[0]
    
    # return atoms filtered by residue id and CA type
    return frame[molname].atoms.filter((rId.is_in(set(sec_str_rids))) & (aName == "CA")))
                                           

if __name__ == "__main__":
    path_to_pdb = "7jzu.pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    sec_str_rids = get_sec_str_residue_ids(frame, molname="A")
    sec_str_CA_atoms = select_sec_str_ca_atoms(frame, molname="A")
