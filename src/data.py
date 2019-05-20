from airfoil import Airfoil
from settings import xfoilpanelingsettings

class Data:
    def __init__(self):
        self.seed_airfoil = Airfoil()
        self.current_airfoil = None
        self.xoffset = 0.
        self.yoffset = 0.
        self.foilscale = 1.

    def readSeedAirfoil(self, fname):
        '''Reads seed airfoil from file.

        Inputs
            fname: file name or path
        Returns
            retval: True on success, False on failure
            errmsg: message describing the error
        '''
        retval, errmsg = self.seed_airfoil.readFromFile(fname)
        return retval, errmsg

    def generate4DigitAirfoil(self, camber, xcamber, thickness):
        '''Generates 4-digit NACA airfoil.

        Inputs
            camber: MCL max camber, as percent of chord
            xcamber: location of max camber, as x/c
            thickness: max thickness, as percent of chord
        '''
        self.seed_airfoil.generate4Digit(camber, xcamber, thickness)

    def generate5DigitAirfoil(self, designation):
        '''Generates 5-digit NACA airfoil.

        Inputs
            designation: 5-digit designation
        Returns
            retval: True on success, False on failure
            errmsg: message describing the error
        '''
        errmsg = "Success"
        if not self.seed_airfoil.generate5Digit(designation):
            errmsg = "Invalid 5-digit designation."
            sys.stderr.write(errmsg+"\n")
            return False, errmsg

        return True, errmsg

    def processSeedAirfoil(self):
        '''Smoots, finds leading edge, and transforms seed airfoil
        
        Returns:
            retval: True on success, False if error
            errmsg: message describing the error
        '''
        retval, errmsg = self.seed_airfoil.smoothPaneling(xfoilpanelingsettings)
        if not retval:
            return retval, errmsg

        retval, errmsg = self.seed_airfoil.findLE()
        if not retval:
            return retval, errmsg

        self.xoffset = -self.seed_airfoil.xle
        self.yoffset = -self.seed_airfoil.yle
        chord, retval, errmsg = self.seed_airfoil.chord()
        if not retval:
            return retval, errmsg
        if chord > 0.:
            self.foilscale = 1./chord
        else:
            errmsg = "Chord must be greater than 0."
            return False, errmsg
        self.seed_airfoil.translate(self.xoffset, self.yoffset)
        self.seed_airfoil.scale(self.foilscale, self.foilscale)
        print(np.amax(self.seed_airfoil.x))

        return True, errmsg

data = Data()
