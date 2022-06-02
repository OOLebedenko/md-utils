import time
from md_utils.frame.calc import get_rmsd
from pyxmolpp2.pipe import TrajectoryProcessor
from pyxmolpp2 import Frame, AtomPredicate


class CalcRmsd(TrajectoryProcessor):

    def __init__(self, reference: Frame, predicate: AtomPredicate, ns_stride: int = 1):
        self.reference = reference
        self.predicate = predicate
        self.ns_stride = ns_stride
        self.time = 1
        self.out_file = None

    def before_first_iteration(self, frame: Frame) -> None:
        # open file for rmsd writing
        timestr = time.strftime("%Y%m%d_%H%M%S")  # time stamp to save file
        self.out_file = open(f"RMSD_{timestr}.csv", "w")

        # create columns in out file
        self.out_file.write("time_ns,RMSD\n")  # write header

    def after_last_iteration(self, exc_type, exc_value, traceback) -> bool:
        # close file
        self.out_file.close()
        self.reference = None
        self.predicate = None
        self.ns_stride = None
        self.time = None
        self.out_file = None
        return False

    def __call__(self, frame: Frame) -> Frame:
        # calculate RMSD
        rmsd = get_rmsd(self.reference, frame, self.predicate)
        # write new line and update time
        self.out_file.write(f"{self.time},{rmsd}\n")
        self.time += self.ns_stride
        return frame


if __name__ == "__main__":
    from md_utils.frame.select import get_sec_str_residues_predicate
    from pyxmolpp2 import Trajectory, PdbFile, TrjtoolDatFile
    from pyxmolpp2.pipe import Run
    import os

    #  setup trajectory parameters
    path_to_trajectory = "path to  your reference"
    path_to_reference = "path to your trajectory"
    trajectory_start = 1
    trajectory_length = 10

    #  load trajectory
    reference = PdbFile(os.path.join(path_to_trajectory, "1_build/ref.pdb"))
    trajectory = Trajectory(reference.frames()[0])
    for ind in range(trajectory_start, trajectory_length + 1):
        fname = "{pattern}.{filetype}".format(pattern="run%05d", filetype="dat")
        trajectory.extend(TrjtoolDatFile(os.path.join(os.path.join(path_to_trajectory, "5_run"), fname % (ind))))

    # get predicate for alignment
    predicate = get_sec_str_residues_predicate(frame=reference, molnames=["A", "B"])

    trajectory | CalcRmsd(reference=reference, predicate=predicate, ns_stride=1) | Run()
