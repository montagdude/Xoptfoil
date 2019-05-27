#!/usr/bin/python3

import os
import sys
import unittest
import subprocess
import numpy as np

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

    def runReference(self, airfoil_file):
        # Run Fortran reference code and read expected results
        truth_file = "data/le_truthdata.dat"
        if os.path.isfile(truth_file):
            os.remove(truth_file)
        os.chdir("reference_code")
        subprocess.run(["./test_airfoil", airfoil_file])
        os.chdir("..")

    def checkFindLE(self, airfoil_file):
        mydata = data.Data()

        truth_file = "data/le_truthdata.dat"
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
        retval, errmsg = mydata.readSeedAirfoil(airfoil_file)
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

    def testFindLE(self):
        '''Tests LE and split location for a few airfoils'''
        self.runReference(os.path.join(os.getcwd(), "data/ag18_plain.dat"))
        self.checkFindLE(os.path.join(os.getcwd(), "data/ag18_plain.dat"))
        self.runReference(os.path.join(os.getcwd(), "data/mh45_labeled.dat"))
        self.checkFindLE(os.path.join(os.getcwd(), "data/mh45_labeled.dat"))
        self.runReference(os.path.join(os.getcwd(), "data/naca_flyingwing.dat"))
        self.checkFindLE(os.path.join(os.getcwd(), "data/naca_flyingwing.dat"))

    def checkSplitAirfoil(self, airfoil_file):
        mydata = data.Data()

        truth_file = "data/splitfoil_truthdata.dat"
        try:
            f = open(truth_file)
        except IOError:
            self.fail("Cannot open {:s}.".format(truth_file))

        # Read top surface
        pointst = int(f.readline())
        xt = np.zeros((pointst))
        yt = np.zeros((pointst))
        for i in range(pointst):
            splitline = f.readline().split()
            xt[i] = float(splitline[0])
            yt[i] = float(splitline[1])

        # Read bottom surface
        pointsb = int(f.readline())
        xb = np.zeros((pointsb))
        yb = np.zeros((pointsb))
        for i in range(pointsb):
            splitline = f.readline().split()
            xb[i] = float(splitline[0])
            yb[i] = float(splitline[1])
        f.close()

        # Read and process seed airfoil
        retval, errmsg = mydata.readSeedAirfoil(airfoil_file)
        if not retval:
            self.fail(errmsg)
        retval, errmsg = mydata.processSeedAirfoil()
        if not retval:
            self.fail(errmsg)

        # Split seed airfoil and compare
        mydata.seed_airfoil.split(symmetrical=False)
        np.testing.assert_almost_equal(mydata.seed_airfoil.xt, xt)
        np.testing.assert_almost_equal(mydata.seed_airfoil.yt, yt)
        np.testing.assert_almost_equal(mydata.seed_airfoil.xb, xb)
        np.testing.assert_almost_equal(mydata.seed_airfoil.yb, yb)

    def testSplitAirfoil(self):
        '''Tests splitting a few airfoils into upper and lower surfaces'''
        self.runReference(os.path.join(os.getcwd(), "data/ag18_plain.dat"))
        self.checkSplitAirfoil(os.path.join(os.getcwd(), "data/ag18_plain.dat"))
        self.runReference(os.path.join(os.getcwd(), "data/mh45_labeled.dat"))
        self.checkSplitAirfoil(os.path.join(os.getcwd(), "data/mh45_labeled.dat"))
        self.runReference(os.path.join(os.getcwd(), "data/naca_flyingwing.dat"))
        self.checkSplitAirfoil(os.path.join(os.getcwd(), "data/naca_flyingwing.dat"))


if __name__ == "__main__":
    unittest.main()
