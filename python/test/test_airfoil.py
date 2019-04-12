#!/usr/bin/python3

import os
import sys
import unittest

cwd = os.getcwd()
sys.path.append(os.path.join(cwd, '..'))

from airfoil import Airfoil

class TestAirfoil(unittest.TestCase):

    def test_readFromFile(self):
        airfoil = Airfoil()

        # Try to read non-existent file
        retval = airfoil.readFromFile("nonsense.dat")
        self.assertEqual(retval, 1)

        # Try to read non-airfoil file
        retval = airfoil.readFromFile("data/not_an_airfoil_file.txt")
        self.assertEqual(retval, 2)

        # Airfoil file with only 1 coordinate
        retval = airfoil.readFromFile("data/one_coordinate.dat")
        self.assertEqual(retval, 2)

        # Airfoil file with missing data
        retval = airfoil.readFromFile("data/missing_data.dat")
        self.assertEqual(retval, 2)

        # Read plain coordinate file
        retval = airfoil.readFromFile("data/ag18_plain.dat")
        self.assertEqual(retval, 0)
        self.assertEqual(airfoil.x.shape[0], 160)

        # Read labeled coordinate file
        retval = airfoil.readFromFile("data/mh45_labeled.dat")
        self.assertEqual(retval, 0)
        self.assertEqual(airfoil.x.shape[0], 67)

if __name__ == "__main__":
    unittest.main()
