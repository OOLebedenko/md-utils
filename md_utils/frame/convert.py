from pyxmolpp2 import Frame


def pdb_to_fasta(frame: Frame) -> None:
    """
    writes one or more sequences from pdb file into a file in FASTA format
    :param path_to_pdb_file:
    :return: file in FASTA format
    """  
    aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    output = []
    for i in range(len(frame.molecules)):
        seq = []
        for mol in frame.molecules[i].residues:
            if mol.name in aa_code:
                seq.append(aa_code[mol.name])

        output.append('>' + str(frame.molecules[i].name))
        output.append(str(''.join(seq)))
    
    print(*output, sep = '\n')

if __name__ == "__main__":
    from pyxmolpp2 import PdbFile

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]
    sequence = pdb_to_fasta(frame=frame)
