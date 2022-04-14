from pyxmolpp2 import Frame, AtomPredicate, calc_rmsd
from select import get_sec_str_residues_predicate
import numpy as np
import sys

def frame_get_rmsd(reference_frame: Frame,
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
    
    # Extract frames by predicates from reference and probe structures

    reference = reference_frame.atoms.frame(by_atoms)
    probe = probe_frame.atoms.frame(by_atoms)
    
    # Similarity estimation by residue - correlation matrix

    prob_ref_alignment = probe.alignment_to(reference)
    
    # Obtain new coordinates considering alignment

    imposition = np.dot(probe.coords.values, prob_ref_alignment.matrix3d().T) + prob_ref_alignment.vector3d().values
    
    return calc_rmsd(reference.coords.values, impositions)



if __name__ == "__main__":
    from pyxmolpp2 import PdbFile
    p1 = sys.argv[1]
    p2 = sys.argv[2]
    
    frame1 = PdbFile("path1").frames()[0]
    frame2 = PdbFile("path2").frames()[0]
    pred = get_sec_str_residues_predicate(frame=frame1, molnames=["B"])