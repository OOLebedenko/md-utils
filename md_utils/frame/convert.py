from pyxmolpp2 import Frame


def pdb_to_fasta(frame: Frame) -> None:
    """
    writes one or more sequences from pdb file into a file in FASTA format
    :param path_to_pdb_file:
    :return: file in FASTA format
    """
    # dictionary to convert canonical residue names
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    sequences = [] # list to save sequences
    
    # iterate over chains in frame
    for mol in frame.molecules:
      name = mol.name # save sequence name
      seq = '' # string to append sequence
      # iterate over residues to save sequence
      for res in mol.residues:
        try:
          r = d[res.name]
          seq += r
        except: # if there are unrecognized residues break and leave an empty record
          seq = ''
          break

      sequences.append([name, seq])

    # save sequences in .fasta
    with open('new_file.fasta' , 'w') as f:
      for mol in sequences:
        if len(mol[1]) > 0: # check if valid sequence was appended
          f.write(f'>{mol[0]}\n')
          f.write(f'{mol[1]}\n')

    return


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile

    # можно дополнительно усовершенствовать функцию
    # передавать в функцию имя файла, чтобы отразить его в фаста

    path_to_pdb = '7jzu.pdb'
    frame = PdbFile(path_to_pdb).frames()[0]
    sequence = pdb_to_fasta(frame=frame)
