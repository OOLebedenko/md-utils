from pyxmolpp2 import Frame


def pdb_to_fasta(frame: Frame) -> None:
    """
    writes one or more sequences from pdb file into a file in FASTA format
    :param path_to_pdb_file:
    :return: file in FASTA format
    """
    # dictionary to convert canonical residue names
    aa_dictionary = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
     'GLH': 'E', 'ASH': 'D', 'HID': 'H', 'HIP': 'H', 'HIE': 'H',
     'CYX': 'C', 'CYM': 'C'}

    sequences = [] # list to save sequences
    
    # filter out protein residues only
    molecules = frame.atoms.residues.filter(rName.is_in(set(aa_dictionary.keys()))).molecules
    
    # iterate over chains in frame
    for mol in molecules:
      name = mol.name # save sequence name
      seq = '' # string to append sequence
      # iterate over residues to save sequence
      for res in mol.residues:
          r = aa_dictionary[res.name]
          seq += r

    # save sequences in .fasta
    with open('new_file.fasta' , 'w') as f:
      for mol in sequences:
        if len(mol[1]) > 0: # check if valid sequence was appended
          f.write(f'>{mol[0]}\n')
          f.write(f'{mol[1]}\n')

    return


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile
    path_to_pdb = '7jzu.pdb'
    frame = PdbFile(path_to_pdb).frames()[0]
    sequence = pdb_to_fasta(frame=frame)
