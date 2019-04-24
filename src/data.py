from airfoil import Airfoil

class Data:
    def __init__(self):
        self.seed_airfoil = Airfoil()
        self.current_airfoil = None

    def readSeedAirfoil(self, fname):
        """Reads seed airfoil from file.

        Inputs
            fname: file name or path
        Returns
            0 on success
            1 on file IO error
            2 on format error
        """
        return self.seed_airfoil.readFromFile(fname)

    def generate4DigitAirfoil(self, camber, xcamber, thickness):
        """Generates 4-digit NACA airfoil.

        Inputs
            camber: MCL max camber, as percent of chord
            xcamber: location of max camber, as x/c
            thickness: max thickness, as percent of chord
        Returns
            True
        """
        npointside = 100
        self.seed_airfoil.generate4Digit(camber, xcamber, thickness, npointside)
        return True

    def generate5DigitAirfoil(self, designation):
        """Generates 5-digit NACA airfoil.

        Inputs
            designation: 5-digit designation
        Returns
            True on success
            False on failure
        """
        npointside = 100
        return self.seed_airfoil.generate5Digit(designation, npointside)

data = Data()
