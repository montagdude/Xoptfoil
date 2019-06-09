#!/usr/bin/python3

import os
import sys
import subprocess
import numpy as np
import unittest
from matplotlib import pyplot as plt

cwd = os.getcwd()
sys.path.append(os.path.join(cwd, '..'))

from airfoil import SeedAirfoil
import data

class TestParametrization(unittest.TestCase):

    def readAirfoil(self, airfoil_file, domaint, domainb, nshapest, nshapesb, shapetype,
                    symmetrical=False):
        # Read and proess seed airfoil
        mydata = data.Data()
        retval, errmsg = mydata.readSeedAirfoil(airfoil_file)
        if not retval:
            self.fail(errmsg)
        retval, errmsg = mydata.processSeedAirfoil()
        if not retval:
            self.fail(errmsg)
        retval, errmsg = mydata.seed_airfoil.setOptimizationLimits(domaint, domainb)
        if not retval:
            self.fail(errmsg)
        retval, errmsg = mydata.seed_airfoil.split(symmetrical)
        if not retval:
            self.fail(errmsg)
        retval, errmsg = mydata.seed_airfoil.setupShapeFunctions(nshapest, nshapesb, shapetype)
        if not retval:
            self.fail(errmsg)

        return mydata.seed_airfoil, mydata.xoffset, mydata.yoffset, mydata.foilscale

    def readTruthAirfoil(self, truth_file):
        try:
            f = open(truth_file)
        except IOError:
            self.fail("Cannot open {:s}.".format(truth_file))

        x = []
        y = []
        for line in f:
            splitline = line.split()
            x.append(float(splitline[0]))
            y.append(float(splitline[1]))
        f.close()

        return np.array(x), np.array(y)

    def test_HicksHenne(self):
        # Shape function parameters
        st = [0.05, -0.03, 0.1, 0.25]
        t1 = [0.1, 0.3, 0.8, 0.5]
        t2 = [0.3, 0.5, 1.0, 2.0]
        nshapest = len(st)
        modest = []
        for i in range(nshapest):
            modest.append(st[i])
            modest.append(t1[i])
            modest.append(t2[i])

        st = [0.15, 0.04, -0.1, -0.03]
        t1 = [0.05, 0.33, 0.6, 0.9]
        t2 = [2.0, 0.3, 0.8, 3.0]
        nshapesb = len(st)
        modesb = []
        for i in range(nshapest):
            modesb.append(st[i])
            modesb.append(t1[i])
            modesb.append(t2[i])

        # Read and process seed airfoil
        domaint = [0.,1.]
        domainb = [0.,1.]
        seed, xoffset, yoffset, foilscale = self.readAirfoil("data/ls417.dat", domaint, 
                                                 domainb, nshapest, nshapesb, "Hicks-Henne")

        # Create a new airfoil
        foil = seed.createPerturbedAirfoil(modest, modesb, "Hicks-Henne", 1., 1., xoffset,
                                           yoffset, foilscale)

        # Run Fortran reference code and read expected results
        truth_file = "data/hickshenne_truthdata.dat"
        if os.path.isfile(truth_file):
            os.remove(truth_file)
        os.chdir("reference_code")
        subprocess.run(["./test_hickshenne", "../data/ls417.dat"])
        os.chdir("..")
        xtruth, ytruth = self.readTruthAirfoil(truth_file)

        # Plot
        fig, ax = plt.subplots()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Hicks-Henne parametrization")
        ax.grid()
        ax.plot(foil.x, foil.y)
        ax.plot(xtruth, ytruth)
        ax.legend(["Xoptfoil", "Reference"])
        plt.show()

        # Compare results
        np.testing.assert_almost_equal(foil.x, xtruth)
        np.testing.assert_almost_equal(foil.y, ytruth)

        # Now change domain and check that shape functions are 0 outside
        domaint = [0.1,0.6]
        domainb = [0.5,0.7]
        seed, xoffset, yoffset, foilscale = self.readAirfoil("data/ls417.dat", domaint, 
                                                 domainb, nshapest, nshapesb, "Hicks-Henne")
        foil = seed.createPerturbedAirfoil(modest, modesb, "Hicks-Henne", 1., 1., xoffset,
                                           yoffset, foilscale)

        nptt = seed.xt.shape[0]
        perturbt = np.zeros((nptt))
        for shape in seed.shapest:
            perturbt += shape.shape
        for i in range(nptt):
            x = seed.xt[i]
            if (x < domaint[0]) or (x > domaint[1]):
                self.assertAlmostEqual(perturbt[i], 0.)

        nptb = seed.xb.shape[0]
        perturbb = np.zeros((nptb))
        for shape in seed.shapesb:
            perturbb += shape.shape
        for i in range(nptb):
            x = seed.xb[i]
            if (x < domainb[0]) or (x > domainb[1]):
                self.assertAlmostEqual(perturbb[i], 0.)

        # Plot
        plt.clf()
        plt.close()
        fig, ax = plt.subplots()
        ax.set_xlabel("x")
        ax.set_ylabel("Hicks-Henne perturbation")
        ax.set_title("Top domain: [{:.2f},{:.2f}], bottom: [{:.2f},{:.2f}]"\
                     .format(domaint[0], domaint[1], domainb[0], domainb[1]))
        ax.grid()
        ax.plot(seed.xt, perturbt)
        ax.plot(seed.xb, perturbb)
        ax.legend(["Top", "Bottom"])
        plt.show()

    def test_NACA(self):
        # Shape function parameters
        modest = [0.1, -0.1, 0.05, 0.6, -0.3, -1.2, 0.03, -0.08, 0.18, -0.25, 0.3,
                  -0.15, 0.33, -0.33, 0.1]
        modesb = [0.12, 0.06, -0.15, -0.2, 0.27, -0.1, -0.06, 0.14, 0.43, 0.15, -0.7,
                  -0.05, 0.13, 0.35, -0.04]
        nshapest = len(modest)
        nshapesb = len(modesb)

        # Read and process seed airfoil
        domaint = [0.,1.]
        domainb = [0.,1.]
        seed, xoffset, yoffset, foilscale = self.readAirfoil("data/ag18_plain.dat", domaint,
                                            domainb, nshapest, nshapesb, "NACA")

        # Create a new airfoil
        foil = seed.createPerturbedAirfoil(modest, modesb, "NACA", 1., 1., xoffset,
                                           yoffset, foilscale)

        # Run Fortran reference code and read expected results
        truth_file = "data/nacaparam_truthdata.dat"
        if os.path.isfile(truth_file):
            os.remove(truth_file)
        os.chdir("reference_code")
        subprocess.run(["./test_nacaparam", "../data/ag18_plain.dat"])
        os.chdir("..")
        xtruth, ytruth = self.readTruthAirfoil(truth_file)

        # Plot
        fig, ax = plt.subplots()
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("NACA parametrization")
        ax.grid()
        ax.plot(foil.x, foil.y)
        ax.plot(xtruth, ytruth)
        ax.legend(["Xoptfoil", "Reference"])
        plt.show()

        # Compare results
        np.testing.assert_almost_equal(foil.x, xtruth)
        np.testing.assert_almost_equal(foil.y, ytruth)

        # Now change domain and check that shape functions are 0 outside
        domaint = [0.8,1.0]
        domainb = [0.1,0.4]
        seed, xoffset, yoffset, foilscale = self.readAirfoil("data/ls417.dat", domaint, 
                                                 domainb, nshapest, nshapesb, "NACA")
        foil = seed.createPerturbedAirfoil(modest, modesb, "NACA", 1., 1., xoffset,
                                           yoffset, foilscale)

        nptt = seed.xt.shape[0]
        perturbt = np.zeros((nptt))
        for shape in seed.shapest:
            perturbt += shape.shape
        for i in range(nptt):
            x = seed.xt[i]
            if (x < domaint[0]) or (x > domaint[1]):
                self.assertAlmostEqual(perturbt[i], 0.)

        nptb = seed.xb.shape[0]
        perturbb = np.zeros((nptb))
        for shape in seed.shapesb:
            perturbb += shape.shape
        for i in range(nptb):
            x = seed.xb[i]
            if (x < domainb[0]) or (x > domainb[1]):
                self.assertAlmostEqual(perturbb[i], 0.)

        # Plot
        plt.clf()
        plt.close()
        fig, ax = plt.subplots()
        ax.set_xlabel("x")
        ax.set_ylabel("NACA perturbation")
        ax.set_title("Top domain: [{:.2f},{:.2f}], bottom: [{:.2f},{:.2f}]"\
                     .format(domaint[0], domaint[1], domainb[0], domainb[1]))
        ax.grid()
        ax.plot(seed.xt, perturbt)
        ax.plot(seed.xb, perturbb)
        ax.legend(["Top", "Bottom"])
        plt.show()


if __name__ == "__main__":
    unittest.main()
