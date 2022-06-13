from md_utils.trj.calc import CalcRmsd
from pyxmolpp2 import mName, aName, Trajectory, PdbFile, AmberNetCDF
from pyxmolpp2.pipe import Run, AssembleQuaternaryStructure
import os
import glob


class TestClass:
    # load trajectory
    reference = PdbFile("../../tests/test_data/reference.pdb").frames()[0]
    trajectory = Trajectory(reference)
    trajectory.extend(AmberNetCDF("../../tests/test_data/sample_trj.nc"))

    # get predicate for Ca atom alignment
    predicate = (aName == "CA")

    # variable for RMSD data
    RMDS_data = []
    # this reference data is computed as Ca RMSD in VMD
    reference_RMDS_data = [1.719, 1.724, 1.767, 1.746, 1.525, 1.498, 1.399, 1.493, 1.647, 1.439]

    def test_one(self):
        # check if CalcRmsd() runs smoothly
        rmsd_file = 'rmds.csv'
        try:
            self.trajectory[::100] | AssembleQuaternaryStructure(of=(mName.is_in("A", "B")), by=self.predicate,
                                                                 reference=self.reference) | CalcRmsd(
                reference=self.reference, predicate=self.predicate) | Run()
            rmsd_file = glob.glob('*.csv')[0]  # save RMSD data for the next test
            with open(rmsd_file, "r") as f:
                for line in f:
                    self.RMDS_data.append(line.strip().split(','))
        except Exception as exc:
            assert False, f"Running CalcRmsd() raised an exception {exc}"
        finally:
            if os.path.exists(rmsd_file):  # delete RMSD.csv
                os.remove(rmsd_file)

    def test_two(self):
        # check if RMSD values in output file are valid
        assert len(self.RMDS_data[1:]) == len(self.reference_RMDS_data)

        for calculated, reference in zip(self.RMDS_data[1:], self.reference_RMDS_data):
            assert abs(float(calculated[1]) - reference) < 0.01
