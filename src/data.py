from airfoil import Airfoil

class Data:
    def __init__(self):
        self.seed_airfoil = Airfoil()

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
