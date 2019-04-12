import sys
import numpy as np

class Airfoil:
    def __init__(self):
        self.x = []
        self.y = []

    def readFromFile(self, fname):
        """Reads airfoil from file.

        Inputs
            fname: file name or path
        Returns
            0 on success
            1 on file IO error
            2 on format error
        """
        
        print("Reading airfoil file {:s}".format(fname))
        try:
            f = open(fname)
        except IOError:
            sys.stderr.write("Unable to read {:s}\n.".format(fname))
            return 1

        # Determine if file has a label
        labeled = False
        splitline = f.readline().split()
        if len(splitline) < 2:
            labeled = True
        else:
            try:
                x = float(splitline[0])
                y = float(splitline[1])
            except ValueError:
                labeled = True
        if labeled:
            print("Labeled coordinates file.")
        else:
            print("Plain coordinates file.")
        f.seek(0)

        # Read coordinates
        linenum = 0
        x = []
        y = []
        for line in f:
            if not labeled or linenum > 0:
                splitline = line.split()
                if len(splitline) < 2:
                    f.close()
                    sys.stderr.write("Not enough data in line.\n")
                    return 2
                try:
                    x.append(float(splitline[0]))
                    y.append(float(splitline[1]))
                except ValueError:
                    f.close()
                    sys.stderr.write("Unable to convert to float.\n")
                    return 2
            linenum += 1
        if len(x) > 1:
            self.x = np.array(x)
            self.y = np.array(y)
        else:
            f.close()
            sys.stderr.write("Less than 2 vertices in coordinates file.\n")
            return 2
        f.close()
        return 0
