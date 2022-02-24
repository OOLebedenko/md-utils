from pyxmolpp2 import Frame


def pdb_to_fasta(frame: Frame) -> None:
    """
    writes one or more sequences from pdb file into a file in FASTA format
    :param path_to_pdb_file:
    :return: file in FASTA format
    """
    ...


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile

    path_to_pdb = "path to your pdb"
    frame = PdbFile(path_to_pdb).frames()[0]
    sequence = pdb_to_fasta(frame=frame)
