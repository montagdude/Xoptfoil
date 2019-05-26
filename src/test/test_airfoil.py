#!/usr/bin/python3

import os
import sys
import unittest
import subprocess

cwd = os.getcwd()
sys.path.append(os.path.join(cwd, '..'))

from airfoil import Airfoil, SeedAirfoil
import data

class TestAirfoil(unittest.TestCase):

    def test_readFromFile(self):
        airfoil = Airfoil()

        # Try to read non-existent file
        retval, errmsg = airfoil.readFromFile("nonsense.dat")
        self.assertEqual(retval, False)
        self.assertIn("Unable to read", errmsg)

        # Try to read non-airfoil file
        retval, errmsg = airfoil.readFromFile("data/not_an_airfoil_file.txt")
        self.assertEqual(retval, False)
        self.assertIn("Unable to convert", errmsg)

        # Airfoil file with only 1 coordinate
        retval, errmsg = airfoil.readFromFile("data/one_coordinate.dat")
        self.assertEqual(retval, False)
        self.assertIn("Less than 3", errmsg)

        # Airfoil file with missing data
        retval, errmsg = airfoil.readFromFile("data/missing_data.dat")
        self.assertEqual(retval, False)
        self.assertIn("Not enough data", errmsg)

        # Read plain coordinate file
        retval, errmsg = airfoil.readFromFile("data/ag18_plain.dat")
        self.assertEqual(retval, True)
        self.assertEqual(airfoil.x.shape[0], 160)

        # Read labeled coordinate file
        retval, errmsg = airfoil.readFromFile("data/mh45_labeled.dat")
        self.assertEqual(retval, True)
        self.assertEqual(airfoil.x.shape[0], 67)

    def testFindLE(self):
        mydata = data.Data()

        # Run Fortran reference code and read expected results
        truth_file = "data/le_truthdata.dat"
        if os.path.isfile(truth_file):
            os.remove(truth_file)
        os.chdir("reference_code")
        subprocess.run(["./test_airfoil"])
        os.chdir("..")
        try:
            f = open(truth_file)
        except IOError:
            self.fail("Cannot open {:s}.".format(truth_file))
        xle = float(f.readline())
        yle = float(f.readline())
        leidx = int(f.readline()) - 1
        addpoint_loc = int(f.readline())
        xoffset = float(f.readline())
        yoffset = float(f.readline())
        foilscale = float(f.readline())
        f.close()

        # Read and process seed airfoil
        retval, errmsg = mydata.readSeedAirfoil("data/mh45_labeled.dat")
        if not retval:
            self.fail(errmsg)
        retval, errmsg = mydata.processSeedAirfoil()
        if not retval:
            self.fail(errmsg)

        # Check LE coordinates and split location
        self.assertAlmostEqual(mydata.seed_airfoil.xle, xle)
        self.assertAlmostEqual(mydata.seed_airfoil.yle, yle)
        self.assertEqual(mydata.seed_airfoil.leidx, leidx)
        self.assertEqual(mydata.seed_airfoil.addpoint_loc, addpoint_loc)

        # Check offset and scale
        self.assertAlmostEqual(mydata.xoffset, xoffset)
        self.assertAlmostEqual(mydata.yoffset, yoffset)
        self.assertAlmostEqual(mydata.foilscale, foilscale)


if __name__ == "__main__":
    unittest.main()
