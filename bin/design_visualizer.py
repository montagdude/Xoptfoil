import argparse
from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

# Default plottiong options

plotoptions = dict(show_seed_airfoil = True,
                   show_seed_airfoil_only = False,
                   plot_airfoils = True,
                   plot_polars = True,
                   save_animation = False,
                   axis_xmin = "auto",
                   axis_xmax = "auto",
                   axis_ymin = "auto",
                   axis_ymax = "auto",
                   axis_cdmin = "auto",
                   axis_cdmax = "auto",
                   axis_clmin = "auto",
                   axis_clmax = "auto")

################################################################################
#
# Airfoil class
#
################################################################################
class Airfoil:

  def __init__(self):
    self.x = np.zeros((0))
    self.y = np.zeros((0))
    self.cl = np.zeros((0))
    self.cd = np.zeros((0))
    self.npt = 0
    self.noper = 0
    
  def setCoordinates(self, x, y):
    self.x = x
    self.y = y
    self.npt = x.shape[0]

  def setPolars(self, cl, cd):
    self.cl = cl
    self.cd = cd
    self.noper = cl.shape[0]

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
# Loads airfoil coordinates and polars from files
def load_airfoils_from_file(coordfilename, polarfilename):

  # Initialize output data

  seedfoil = Airfoil()
  designfoils = []
  ioerror = 0

  # Read seed airfoil coordinates

  zonetitle = 'zone t="Seed airfoil"'
  x, y, ioerror = read_airfoil_data(coordfilename, zonetitle)
  if (ioerror == 1):
    print("Error: file " + coordfilename + " not found.")
    return seedfoil, designfoils, ioerror
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + coordfilename
          + ".")
    return seedfoil, designfoils, ioerror

  seedfoil.setCoordinates(np.array(x), np.array(y))

  # Read seed airfoil polars

  zonetitle = 'zone t="Seed airfoil polar"'
  cl, cd, ioerror = read_airfoil_data(polarfilename, zonetitle)
  if (ioerror == 1):
    print("Error: file " + polarfilename + " not found.")
    return seedfoil, designfoils, ioerror
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + polarfilename
          + ".")
    return seedfoil, designfoils, ioerror

  seedfoil.setPolars(np.array(cl), np.array(cd))

  # Read coordinate data for designs produced by optimizer

  print("Reading airfoil coordinates from file " + coordfilename + "...")

  read_finished = False
  counter = 1
  while (not read_finished):

    zonetitle = 'zone t="Airfoil", SOLUTIONTIME=' + str(counter)
    x, y, ioerror = read_airfoil_data(coordfilename, zonetitle)
    if (ioerror == 2):
      read_finished = True
      numfoils = counter - 1
      ioerror = 0
    else:
      currfoil = Airfoil()
      currfoil.setCoordinates(np.array(x), np.array(y))
      designfoils.append(currfoil)
      counter += 1

  print("Found " + str(numfoils) + " airfoil coordinates plus seed airfoil.")

  # Read polar data for designs produced by optimizer

  print("Reading airfoil polars from file " + polarfilename + "...")

  read_finished = False
  counter = 1
  while (not read_finished):

    zonetitle = 'zone t="Polars", SOLUTIONTIME=' + str(counter)
    cl, cd, ioerror = read_airfoil_data(polarfilename, zonetitle)
    if (ioerror == 2):
      read_finished = True
      ioerror = 0
    else:
      designfoils[counter-1].setPolars(np.array(cl), np.array(cd))
      counter += 1

  print("Found " + str(counter-1) + " airfoil polars plus seed airfoil.")
  if (counter != numfoils+1):
    print("Error: number of airfoil coordinates and polars does not match.")
    ioerror = 3 
    return seedfoil, designfoils, ioerror

  return seedfoil, designfoils, ioerror

################################################################################
# Plots an airfoil + polars
def plot_airfoil(seedfoil, designfoils, plotnum):
  global plotoptions

  # Select requested airfoil

  foil = designfoils[plotnum-1]

  # Set up plot

  plt.cla()
  plt.clf()
  plt.close()

  if ( (plotoptions["plot_airfoils"]) and (plotoptions["plot_polars"]) ):
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,2])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
  elif ( (plotoptions["plot_airfoils"]) and (not plotoptions["plot_polars"]) ):
    ax0 = plt.subplot(111)
  elif ( (not plotoptions["plot_airfoils"]) and (plotoptions["plot_polars"]) ):
    ax1 = plt.subplot(111)
  else:
    print("Error: must enable at least one of airfoils plot or polars plot.")
    return

  if (plotoptions["show_seed_airfoil_only"]):
    xmax = np.max(seedfoil.x)
    xmin = np.min(seedfoil.x)
  elif (plotoptions["show_seed_airfoil"]):
    xmax = max([np.max(seedfoil.x), np.max(foil.x)])
    xmin = min([np.min(seedfoil.x), np.min(foil.x)])
  else:
    xmax = np.max(foil.x)
    xmin = np.min(foil.x)
  xrng = xmax - xmin
  
  if (plotoptions["show_seed_airfoil_only"]):
    cdmax = np.max(seedfoil.cd)
    cdmin = np.min(seedfoil.cd)
  elif (plotoptions["show_seed_airfoil"]):
    cdmax = max([np.max(seedfoil.cd), np.max(foil.cd)])
    cdmin = min([np.min(seedfoil.cd), np.min(foil.cd)])
  else:
    cdmax = np.max(foil.cd)
    cdmin = np.min(foil.cd)
  cdrng = cdmax - cdmin

  # Plot airfoil coordinates

  if (plotoptions["plot_airfoils"]):
    if (plotoptions["show_seed_airfoil_only"]):
      ax0.plot(seedfoil.x, seedfoil.y, '-b')
    elif (plotoptions["show_seed_airfoil"]):
      ax0.plot(seedfoil.x, seedfoil.y, '-b', foil.x, foil.y, '-r')
    else:
      ax0.plot(foil.x, foil.y, '-r')
    ax0.set_aspect('equal', 'datalim')
    ax0.set_xlabel('x')
    ax0.set_ylabel('y')
    ax0.set_xlim([xmin-0.05*xrng,xmax+0.05*xrng])

  # Plot polars

  if (plotoptions["plot_polars"]):
    if (plotoptions["show_seed_airfoil_only"]):
      ax1.plot(seedfoil.cd, seedfoil.cl, '-ob')
    elif (plotoptions["show_seed_airfoil"]):
      ax1.plot(seedfoil.cd, seedfoil.cl, '-ob', foil.cd, foil.cl, '-sr') 
    else: 
      ax1.plot(foil.cd, foil.cl, '-sr') 
    ax1.set_xlabel('Drag coefficient')
    ax1.set_ylabel('Lift coefficient')
    ax1.set_xlim([cdmin-0.05*cdrng,cdmax+0.05*cdrng])
    ax1.grid()

  # Legend

  if (plotoptions["plot_airfoils"]): legendax = ax0
  else: legendax = ax1

  if (plotoptions["show_seed_airfoil_only"]):
    legendax.legend(['Seed airfoil'], loc="upper center", 
                    bbox_to_anchor=(0.5,1.4))
  elif (plotoptions["show_seed_airfoil"]):
    legendax.legend(['Seed airfoil', 'Design number ' + str(plotnum)], 
                    loc="upper center", bbox_to_anchor=(0.5,1.4))
  else:
    legendax.legend(['Design number ' + str(plotnum)], loc="upper center",
                    bbox_to_anchor=(0.5,1.4))

  plt.show()

################################################################################
# Plotting menu
def plotting_menu(seedfoil, designfoils):

  numfoils = len(designfoils)

  plotting_complete = False
  validchoice = False
  while (not validchoice):
    print()
    print("There are " + str(numfoils) + " designs.")
    plotnum = int(input("Enter design to plot (or 0 to return): "))

    # Return to main menu

    if (plotnum == 0): 
      validchoice = True
      plotting_complete = True

    # Check for bad index

    elif ( (plotnum < 1) or (plotnum > numfoils) ):
      validchoice = False
      print("Error: index out of bounds.")

    # Plot design

    else:
      validchoice = True
      plot_airfoil(seedfoil, designfoils, plotnum)
      plotting_complete = False

  return plotting_complete

################################################################################
# Gets boolean input from user
def get_boolean_input(key, keyval):

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval))
    print("Available choices: True, False\n")
    sel = input("Enter new value: ")
    if ( (sel == "True") or (sel == "true")):
      retval = True
      validchoice = True
    elif ( (sel == "False") or (sel == "false")):
      retval = False
      validchoice = True
    else:
      print("Please enter True or False.")
      validchoice = False

  return retval

################################################################################
# Options menu: allows user to change plot options
def options_menu():
  global plotoptions

  # Status variable

  options_complete = False
 
  # Print list of plotting options

  print()
  print("Available plotting options:")
  print()
  for key in sorted(plotoptions):
    print(key + " [" + str(plotoptions[key]) + "]")
  print()

  # Get user input

  key = input("Enter option to change (or 0 to return): ")

  # True/False settings

  if ( (key == "show_seed_airfoil") or (key == "show_seed_airfoil_only") or
       (key == "save_animation") or (key == "plot_airfoils") or 
       (key == "plot_polars") ):
    options_complete = False
    plotoptions[key] = get_boolean_input(key, plotoptions[key])

  # Change plot bounds

  elif ( (key == "axis_xmax") or (key == "axis_xmin") or
         (key == "axis_ymax") or (key == "axis_ymin") or
         (key == "axis_cdmax") or (key == "axis_cdmin") or
         (key == "axis_clmax") or (key == "axis_clmin") ):
    options_complete = False
    sel = input("Enter new value for " + key + " or enter 'auto': ")
    plotoptions[key] = sel 

  # Exit options menu

  elif (key == "0"): 
    options_complete = True

  # Error for invalid input

  else: 
    options_complete = False
    print("Unrecognized plot option.")

  return options_complete

################################################################################
# Main menu
def main_menu(seedfoil, designfoils):

  exitchoice = False
  while (not exitchoice):
    
    print()
    print("Options:")
    print("[0] Exit")
    print("[1] Plot a specific design")
    print("[2] Animate all designs")
    print("[3] Change plotting options")
    print()

    choice = input("Enter a choice [0-3]: ")

    if (choice == "0"):
      exitchoice = True

    elif (choice == "1"):
      exitchoice = False
      plotting_complete = False
      while (not plotting_complete): plotting_complete = plotting_menu(
                                                          seedfoil, designfoils)

    elif (choice == "3"):
      exitchoice = False
      options_complete = False
      while (not options_complete): options_complete = options_menu()

    else:
      print("Error: please enter a choice 0-3.")

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

  # Read airfoil coordinates and polars

  seedfoil, designfoils, ioerror = load_airfoils_from_file(
                                                   coordfilename, polarfilename)

  # Call main menu

  if (ioerror == 0): main_menu(seedfoil, designfoils)
