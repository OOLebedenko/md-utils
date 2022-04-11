from pyxmolpp2 import Frame, AtomPredicate


def calc_rmsd(reference_frame: Frame,
              probe_frame: Frame,
              by_atoms: AtomPredicate,
              ) -> float:
    """
    Comparison similarity in teo three-dimensional structures by the RMSD after optimal rigid body superposition
    https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions

    :param reference_frame: for example X-ray structure or any other case
    :param probe_frame: for example frame from MD trajectory or any other case
    :param by_atoms: rigid body superposition is carried out by this atoms. In case of proteins,
                     this often calculated for the secondary-structure Cα atoms.
    :return: rmsd value by_atoms between reference_frame and probe_frame
    """
