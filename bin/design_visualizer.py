import argparse
from matplotlib import pyplot as plt
import numpy as np

################################################################################
# Reads airfoil data (either coordinates or polars) from file
def read_airfoil_data(filename, zonetitle):

  ioerror = 0
  x = []
  y = []

  # Try to open the file

  try:
    f = open(filename) 
  except IOError:
    ioerror = 1
    return x, y, ioerror

  # Read lines until we get to the correct zone

  zonefound = False
  zonelen = len(zonetitle)

  for textline in f:

    if (not zonefound):

      # Check for the zone we are looking for

      if (textline[0:zonelen] == zonetitle):
    
        zonefound = True

    else:

      # Check to see if we've read all the coordinates

      if (textline[0:4] == "zone"): break

      # Otherwise read coordinates

      else:

        line = textline.split()
        x.append(float(line[0]))
        y.append(float(line[1]))

  # Error if zone has not been found after reading the file

  if (not zonefound):
    f.close()
    ioerror = 2

  # Return coordinate data

  return x, y, ioerror

################################################################################
# Main design_visualizer program
if __name__ == "__main__":

  # Read command line arguments

  parser = argparse.ArgumentParser()
  parser.add_argument("prefix", nargs='?', help="Output file prefix in "
        + "XoptFoil (preceding _design_coordinates.dat and _design_polars.dat)")
  args = parser.parse_args()

  # Determine names of data files

  if (args.prefix):
    coordfilename = args.prefix + '_design_coordinates.dat'
    polarfilename = args.prefix + '_design_polars.dat'
  else:
    coordfilename = 'design_coordinates.dat'
    polarfilename = 'design_polars.dat'

  # Read airfoil coordinate data

  print("Reading airfoil coordinates from file " + coordfilename + "...")

  zonetitle = 'zone t="Seed airfoil"'
  x, y, ioerror = read_airfoil_data(coordfilename, zonetitle)
  if (ioerror == 1):
    print("Error: file " + coordfilename + " not found.")
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + coordfilename
          + ".")

  # Read airfoil polar data

  print("Reading airfoil polars from file " + polarfilename + "...")

  zonetitle = 'zone t="Seed airfoil polar"'
  x, y, ioerror = read_airfoil_data(polarfilename, zonetitle)
  if (ioerror == 1):
    print("Error: file " + polarfilename + " not found.")
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + polarfilename
          + ".")
