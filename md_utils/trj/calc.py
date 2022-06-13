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
    # test CalcRMSD using pytest
    import pytest
    CalcRMSD_retcode = pytest.main(["../../tests/test_trj_rmsd.py"])
    print(CalcRMSD_retcode)
