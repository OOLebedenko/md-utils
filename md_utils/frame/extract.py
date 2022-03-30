from pyxmolpp2 import (AtomPredicate, Frame, rName, aName, mName, Molecule,
                       MoleculeSelection, ResidueSelection, AtomSelection, Atom, aId)
from scipy.spatial import cKDTree
from typing import Dict, List, Union, Tuple
import numpy as np
import functools
import operator
import json


def extract_one_letter_amino_acid_seq(frame: Frame,
                                      molname: mName) -> list:
    """
    extract one letter amino asids sequence for the specific molecule in frame
    :param frame:
    :param molname: name of molecule in frame
    :return: list of one letter amino acids for the specific molecule in frame
    """
    # dictionary to convert canonical residue names
    three_to_one_leter_name = {"GLY": "G",
                               "ALA": "A",
                               "VAL": "V",
                               "LEU": "L",
                               "ILE": "I",
                               "SER": "S",
                               "THR": "T",
                               "MET": "M",
                               "CYS": "C",
                               "ASN": "N",
                               "GLN": "Q",
                               "ASP": "D",
                               "GLU": "E",
                               "LYS": "K",
                               "ARG": "R",
                               "HIS": "H",
                               "PHE": "F",
                               "TRP": "W",
                               "TYR": "Y",
                               "PRO": "P",
                               }

    return [three_to_one_leter_name[residue.name] for residue in frame.molecules.filter(mName == molname).residues]


def extract_residues_on_interface(partner_A: Union[Molecule, MoleculeSelection],
                                  partner_B: Union[Molecule, MoleculeSelection],
                                  cutoff: float
                                  ) -> List[ResidueSelection]:
    """
    Extract residues on interface of complex between partner_A and partner_B
    :param partner_A: first partner of intermolecular interaction
    :param partner_B: second partner of intermolecular interaction
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
    resdue_selection_A = partner_A.atoms.filter(aId.is_in(set(a_index + first_a_id))).residues

    first_b_id = partner_B.atoms[0].id  # to correct indexes from matrix
    resdue_selection_B = partner_B.atoms.filter(aId.is_in(set(b_index + first_b_id))).residues

    return [resdue_selection_A, resdue_selection_B]


def get_predicate_donors(donors_dict: Dict[str, list]
                         ) -> AtomPredicate:
    donors_predicate = []
    for rname in donors_dict:
        atoms_donor = []
        for atom_pair in donors_dict[rname]:
            atoms_donor.append(aName == atom_pair[0])
        donors_predicate.append((rName == rname) & (functools.reduce(operator.or_, atoms_donor)))

    return functools.reduce(operator.or_, donors_predicate)

def get_predicate_acceptors(acceptor_dict: Dict[str, list]
                            ) -> AtomPredicate:
    acceptor_predicate = []
    for ind, rname in enumerate(acceptor_dict):
        acceptor_predicate.append((rName == rname) & (aName.is_in(set(acceptor_dict[rname]))))

    return functools.reduce(operator.or_, acceptor_predicate)

def extract_hydrogen_bonds(partner_A: Union[Molecule, MoleculeSelection],
                           partner_B: Union[Molecule, MoleculeSelection],
                           distance_cutoff=3,
                           angle_cutoff=30,
                           ) -> List[Tuple[Union[Atom, AtomSelection]]]:
    """
    extract hydrogen bond is formed between an atom with a hydrogen bonded to it (the donor, D)
            and another atom (the acceptor, A)
    :param partner_A:
    :param partner_B:
    :param distance_cutoff: distance D-A is less than the cut-off distance (default 3.0 Angstroms)
    :param angle_cutoff: angle D-H-A is less than the cut-off angle (default 30 degrees)
    :return:
    """
    # load data about donors and acceptors
    with open('./supporting_data/HydrogenBondSites.json', 'r') as f:
      hydrogen_bonds_sites = json.load(f)

    donors_dict = hydrogen_bonds_sites["HydrogenBondDonors"]
    acceptors_dict = hydrogen_bonds_sites["HydrogenBondAcceptors"]

    # get donor and acceptors predicates
    donors_predicate = get_predicate_donors(donors_dict)
    acceptors_predicate = get_predicate_acceptors(acceptors_dict)

    # check for rosseta notation in atom names
    general_notation_h = set({'HD21', 'HD22', 'HH11', 'HH12', 'H21', 'H22', 'HE21', 'HE22',
                              'HZ1', 'HZ2', 'HZ3', 'H1', 'H2', 'H3'})
    if len(partner_A.atoms.filter(aName.is_in(general_notation_h))) == 0:
      rosetta_notation = True
    else:
      rosetta_notation = False
    
    def get_hydrogens(donor_atom): # retreive hydrogens attached to the donor atom
      hydrogens = []
      dh_pairs =  hydrogen_bonds_sites['HydrogenBondDonors'][donor_atom.residue.name]
      for pair in dh_pairs:
        if pair[0] == donor_atom.name:
          if rosetta_notation:
            h_rosetta_name = pair[1][-1] + pair[1][:-1]
            new_h = donor_atom.residue.atoms.filter(aName == h_rosetta_name)
          else:
            new_h = donor_atom.residue.atoms.filter(aName == pair[1])
          if len(new_h) > 0:
            hydrogens.append(new_h[0])
      return hydrogens

    def angle(a, b, h): # calculate criteria angle
      l1 = (h-a) / np.linalg.norm(h-a) # normalization of vectors
      l2 = (h-b) / np.linalg.norm(h-b) # normalization of vectors
      cos = np.dot(l1, l2)
      rad = np.arccos(np.clip(cos, -1.0, 1.0))
      deg = np.rad2deg(rad)
      return 180 - deg
    
    def get_donor_acceptor_pairs(donor_chain, acceptor_chain):
      # specify donor and acceptor atoms
      donor_atoms = donor_chain.atoms.filter(donors_predicate)
      acceptor_atoms = acceptor_chain.atoms.filter(acceptors_predicate)

      # get coordinates donor/acceptor atoms and construct binary trees
      donor_tree = cKDTree(donor_atoms.coords.values)
      acceptor_tree = cKDTree(acceptor_atoms.coords.values)
      
      # compute distances with a cut-off
      dist = donor_tree.sparse_distance_matrix(acceptor_tree, max_distance=distance_cutoff)
      dist = dist.toarray()
      donor_sites, acceptor_sites = dist.nonzero()
      
      # retrieve close atoms
      selected_donors, selected_acceptors = [], []
      for ind in a_sites:
        selected_donors.append(donor_atoms[ind])
      for ind in b_sites:
        selected_acceptors.append(acceptor_atoms[ind])

      # for each pair find donor's hydrogens and calculate angle
      # if angle does not exceed threshold - add pair to output
      h_bonded = []

      for don_atom, acc_atom in zip(selected_donors, selected_acceptors):
        hydrogens = get_hydrogens(don_atom)
        if len(hydrogens) > 0:
          for h in hydrogens:
            if angle(don_atom.r.values, acc_atom.r.values, h.r.values) <= angle_cutoff:
              h_bonded.append((don_atom, acc_atom))

      return h_bonded

    return get_donor_acceptor_pairs(partner_A, partner_B) + get_donor_acceptor_pairs(partner_B, partner_A)


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, mName

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]

    # example usage of extract_residues_on_interface
    interface_residues_partner_A, interface_residues_partner_B = extract_residues_on_interface(
        partner_A=frame.molecules.filter(mName == "A"),
        partner_B=frame.molecules.filter(mName == "B"),
        cutoff=5)
    
    # example usage of extract_hydrogen_bonds
    h_bonds = extract_hydrogen_bonds(partner_A=interface_residues_partner_A,
                           partner_B=interface_residues_partner_B,
                           distance_cutoff=4,
                           angle_cutoff=30)
    
    # example usage of extract_one_letter_amino_acid_seq
    amino_acid_seq = extract_one_letter_amino_acid_seq(frame, molname="A")
