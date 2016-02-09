import argparse
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import rcParams
import numpy as np
from math import log10, floor
from sys import version_info
from os import remove

# Default plottiong options

plotoptions = dict(show_seed_airfoil = True,
                   show_seed_polar = True,
                   show_seed_airfoil_only = False,
                   show_seed_polar_only = False,
                   plot_airfoils = True,
                   plot_polars = True,
                   save_animation_frames = False,
                   axis_xmin = "auto",
                   axis_xmax = "auto",
                   axis_ymin = "auto",
                   axis_ymax = "auto",
                   axis_cdmin = "auto",
                   axis_cdmax = "auto",
                   axis_clmin = "auto",
                   axis_clmax = "auto",
                   color_for_seed = "blue",
                   color_for_new_designs = "red",
                   monitor_update_frequency = 10)

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
    ioerror = 2

  # Close the file

  f.close()

  # Return coordinate data

  return x, y, ioerror

################################################################################
# Loads airfoil coordinates and polars from files
def load_airfoils_from_file(coordfilename, polarfilename):

  # Initialize output data

  seedfoil = Airfoil()
  designfoils = []
  ioerror = 0

  # Check for seed airfoil coordinates

  print("Checking for airfoil coordinates file " + coordfilename + "...")

  zonetitle = 'zone t="Seed airfoil"'
  x, y, ioerror = read_airfoil_data(coordfilename, zonetitle)
  if (ioerror == 1):
    print("Warning: file " + coordfilename + " not found.")
    return seedfoil, designfoils, ioerror
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + coordfilename
          + ".")
    return seedfoil, designfoils, ioerror

  seedfoil.setCoordinates(np.array(x), np.array(y))

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

  # Read seed airfoil polars (note: negative error code means coordinates were
  # read but not polars)

  print("Checking for airfoil polars file " + coordfilename + "...")

  zonetitle = 'zone t="Seed airfoil polar"'
  cl, cd, ioerror = read_airfoil_data(polarfilename, zonetitle)
  if (ioerror == 1):
    print("Warning: file " + polarfilename + " not found.")
    return seedfoil, designfoils, 0 - ioerror
  elif (ioerror == 2):
    print("Error: zone labeled " + zonetitle + " not found in " + polarfilename
          + ".")
    return seedfoil, designfoils, 0 - ioerror

  seedfoil.setPolars(np.array(cl), np.array(cd))

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
def plot_airfoil(seedfoil, designfoils, plotnum, firsttime=True, 
                 animation=False, prefix=None):
  global plotoptions

  # Select requested airfoil

  if (plotnum > 0): foil = designfoils[plotnum-1]

  # Set up plot

  plt.cla()
  plt.clf()
  if (firsttime): plt.close()

  if (firsttime):
    fig = plt.figure()
  else:
    fig = plt.gcf()

  if ( (plotoptions["plot_airfoils"]) and (plotoptions["plot_polars"]) ):
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,2])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    numplots = 2
  elif ( (plotoptions["plot_airfoils"]) and 
         (not plotoptions["plot_polars"]) ):
    ax0 = plt.subplot(111)
    numplots = 1
  elif ( (not plotoptions["plot_airfoils"]) and 
         (plotoptions["plot_polars"]) ):
    ax1 = plt.subplot(111)
    numplots = 1
  else:
    print("Error: must enable at least one of airfoils plot or polars plot.")
    return

  # Set auto plotting bounds

  if (plotoptions["plot_airfoils"]):
    if ( (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
      xmax = np.max(seedfoil.x)
      xmin = np.min(seedfoil.x)
    elif (plotoptions["show_seed_airfoil"]):
      xmax = max([np.max(seedfoil.x), np.max(foil.x)])
      xmin = min([np.min(seedfoil.x), np.min(foil.x)])
    else:
      xmax = np.max(foil.x)
      xmin = np.min(foil.x)
    xrng = xmax - xmin
    xmaxauto = xmax + 0.1*xrng
    xminauto = xmin - 0.1*xrng
  
    if ( (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
      ymax = np.max(seedfoil.y)
      ymin = np.min(seedfoil.y)
    elif (plotoptions["show_seed_airfoil"]):
      ymax = max([np.max(seedfoil.y), np.max(foil.y)])
      ymin = min([np.min(seedfoil.y), np.min(foil.y)])
    else:
      ymax = np.max(foil.y)
      ymin = np.min(foil.y)
    yrng = ymax - ymin
    ymaxauto = ymax + 0.1*yrng
    yminauto = ymin - 0.1*yrng
  
  if (plotoptions["plot_polars"]):
    if ( (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
      cdmax = np.max(seedfoil.cd)
      cdmin = np.min(seedfoil.cd)
    elif (plotoptions["show_seed_polar"]):
      cdmax = max([np.max(seedfoil.cd), np.max(foil.cd)])
      cdmin = min([np.min(seedfoil.cd), np.min(foil.cd)])
    else:
      cdmax = np.max(foil.cd)
      cdmin = np.min(foil.cd)
    cdrng = cdmax - cdmin
    cdmaxauto = cdmax + 0.1*cdrng
    cdminauto = cdmin - 0.1*cdrng
  
    if ( (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
      clmax = np.max(seedfoil.cl)
      clmin = np.min(seedfoil.cl)
    elif (plotoptions["show_seed_polar"]):
      clmax = max([np.max(seedfoil.cl), np.max(foil.cl)])
      clmin = min([np.min(seedfoil.cl), np.min(foil.cl)])
    else:
      clmax = np.max(foil.cl)
      clmin = np.min(foil.cl)
    clrng = clmax - clmin
    clmaxauto = clmax + 0.1*clrng
    clminauto = clmin - 0.1*clrng

  # Set user-specified plotting bounds and check for errors

  if (plotoptions["plot_airfoils"]):
    if (plotoptions["axis_xmax"] != "auto"):
      xmax = float(plotoptions["axis_xmax"])
    else:
      xmax = xmaxauto
    if (plotoptions["axis_xmin"] != "auto"):
      xmin = float(plotoptions["axis_xmin"])
    else:
      xmin = xminauto
    if (xmin >= xmax):
      print("Warning: xmin must be less than xmax. " +
            "Reverting to auto x bounds.")
      xmin = xminauto
      xmax = xmaxauto
  
    if (plotoptions["axis_ymax"] != "auto"):
      ymax = float(plotoptions["axis_ymax"])
    else:
      ymax = ymaxauto
    if (plotoptions["axis_ymin"] != "auto"):
      ymin = float(plotoptions["axis_ymin"])
    else:
      ymin = yminauto
    if (ymin >= ymax):
      print("Warning: ymin must be less than ymax. " + 
            "Reverting to auto y bounds.")
      ymin = yminauto
      ymax = ymaxauto

  if (plotoptions["plot_polars"]):
    if (plotoptions["axis_cdmax"] != "auto"):
      cdmax = float(plotoptions["axis_cdmax"])
    else:
      cdmax = cdmaxauto
    if (plotoptions["axis_cdmin"] != "auto"):
      cdmin = float(plotoptions["axis_cdmin"])
    else:
      cdmin = cdminauto
    if (cdmin >= cdmax):
      print("Warning: cdmin must be less than cdmax. " +
            "Reverting to auto cd bounds.")
      cdmin = cdminauto
      cdmax = cdmaxauto
  
    if (plotoptions["axis_clmax"] != "auto"):
      clmax = float(plotoptions["axis_clmax"])
    else:
      clmax = clmaxauto
    if (plotoptions["axis_clmin"] != "auto"):
      clmin = float(plotoptions["axis_clmin"])
    else:
      clmin = clminauto
    if (clmin >= clmax):
      print("Warning: clmin must be less than clmax. " +
            "Reverting to auto cl bounds.")
      clmin = clminauto
      clmax = clmaxauto

  # Aliases for colors

  sc = plotoptions["color_for_seed"]
  nc = plotoptions["color_for_new_designs"]

  # Plot airfoil coordinates

  if (plotoptions["plot_airfoils"]):
    if ( (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
      ax0.plot(seedfoil.x, seedfoil.y, color=sc)
    elif (plotoptions["show_seed_airfoil"]):
      ax0.plot(seedfoil.x, seedfoil.y, color=sc)
      ax0.plot(foil.x, foil.y, color=nc)
    else:
      ax0.plot(foil.x, foil.y, color=nc)
    ax0.set_aspect('equal', 'datalim')
    ax0.set_xlabel('x')
    ax0.set_ylabel('z')
    ax0.set_xlim([xmin,xmax])
    ax0.set_ylim([ymin,ymax])

  # Plot polars

  if (plotoptions["plot_polars"]):
    if ( (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
      ax1.plot(seedfoil.cd, seedfoil.cl, linestyle='-', color=sc, marker='o')
    elif (plotoptions["show_seed_polar"]):
      ax1.plot(seedfoil.cd, seedfoil.cl, linestyle='-', color=sc, marker='o')
      ax1.plot(foil.cd, foil.cl, linestyle='-', color=nc, marker='s') 
    else: 
      ax1.plot(foil.cd, foil.cl, linestyle='-', color=nc, marker='s') 
    ax1.set_xlabel('Drag coefficient')
    ax1.set_ylabel('Lift coefficient')
    ax1.set_xlim([cdmin,cdmax])
    ax1.set_ylim([clmin,clmax])
    ax1.grid()

  # Legend 

  if (plotoptions["plot_airfoils"]): legendax = ax0
  else: legendax = ax1

  # Legend placement

  if (numplots == 2):
    bbox_loc = (0.5, 1.4)
  else:
    bbox_loc = (0.5, 1.1)

  # Fake lines for legend

  lines = []
  if ( (plotoptions["plot_airfoils"]) and (plotoptions["plot_polars"]) ):

    # Plotting airfoils and polars

    if ( (plotoptions["show_seed_airfoil"]) or 
         (plotoptions["show_seed_airfoil_only"]) or
         (plotoptions["show_seed_polar"]) or
         (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
      fakeline = plt.Line2D((0,1),(0,0), color=sc, label="Seed airfoil")
      lines.append(fakeline)
    if (plotnum != 0 and ( (not plotoptions["show_seed_airfoil_only"]) or
         (not plotoptions["show_seed_polar_only"]) )):
      fakeline = plt.Line2D((0,1),(0,0), color=nc, 
                            label="Design number " + str(plotnum))
      lines.append(fakeline)

  elif ( (plotoptions["plot_airfoils"]) and (not plotoptions["plot_polars"]) ):

    # Only plotting airfoils

    if ( (plotoptions["show_seed_airfoil"]) or 
         (plotoptions["show_seed_airfoil_only"]) or (plotnum == 0) ):
      fakeline = plt.Line2D((0,1),(0,0), color=sc, label="Seed airfoil")
      lines.append(fakeline)
    if ( (not plotoptions["show_seed_airfoil_only"]) and (plotnum != 0) ):
      fakeline = plt.Line2D((0,1),(0,0), color=nc, 
                            label="Design number " + str(plotnum))
      lines.append(fakeline)

  else:

    # Only plotting polars

    if ( (plotoptions["show_seed_polar"]) or 
         (plotoptions["show_seed_polar_only"]) or (plotnum == 0) ):
      fakeline = plt.Line2D((0,1),(0,0), linestyle='-', color=sc, marker='o',
                            label="Seed airfoil")
      lines.append(fakeline)
    if ( (not plotoptions["show_seed_polar_only"]) and (plotnum != 0) ):
      fakeline = plt.Line2D((0,1),(0,0), linestyle='-', color=nc, marker='s', 
                            label="Design number " + str(plotnum))
      lines.append(fakeline)

  # Create legend

  labels = [l.get_label() for l in lines]
  legendax.legend(lines, labels, loc="upper center", bbox_to_anchor=bbox_loc, 
                  numpoints=1)

  # Plot for the first time, or redraw

  if (not animation): plt.show()
  else:
    if (firsttime): fig.show()
    else: plt.pause(0.0001)
    fig.canvas.draw()

  # Save animation frames if requested

  if ( (animation) and (plotoptions["save_animation_frames"]) ):
    if (prefix == None):
      print("Error: no file prefix specified - cannot save animation frames.")
    else:
      imagefname = prefix + '.png'
      print("Saving image frame to file " + imagefname + ' ...')
      plt.savefig(imagefname)

################################################################################
# Input function that checks python version
def my_input(message):

  # Check python version

  python_version = version_info[0]

  # Issue correct input command

  if (python_version == 2):
    return raw_input(message)
  else:
    return input(message)

################################################################################
# Plotting menu
def plotting_menu(seedfoil, designfoils):

  numfoils = len(designfoils)

  plotting_complete = False
  validchoice = False
  while (not validchoice):
    print("")
    print("There are " + str(numfoils) + " designs.")
    plotnum = int(my_input("Enter design to plot (or 0 to return): "))

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
      plot_airfoil(seedfoil, designfoils, plotnum, firsttime=True)
      plotting_complete = False

  return plotting_complete

################################################################################
# Monitors airfoil coordinates and polar files for updates during optimization
def monitor_airfoil_data(seedfoil, designfoils, prefix):

  # Temporary airfoil struct

  foil = Airfoil()

  # Set up file names to monitor

  coordfilename = prefix + '_design_coordinates.dat'
  polarfilename = prefix + '_design_polars.dat'

  # Determine which design to read for coordinate file

  if (seedfoil.npt == 0):
    zonetitle = 'zone t="Seed airfoil"'
    foilstr = 'seed'
  else:
    nextdesign = len(designfoils) + 1
    zonetitle = 'zone t="Airfoil", SOLUTIONTIME=' + str(nextdesign)
    foilstr = 'design number ' + str(nextdesign)

  # Read data from coordinate file

  x, y, ioerror = read_airfoil_data(coordfilename, zonetitle)
  if (ioerror == 1):
    print("Airfoil coordinates file " + coordfilename + " not available yet.")
    return seedfoil, designfoils, ioerror
  elif (ioerror == 2):
    print("Waiting for next airfoil design.")
    return seedfoil, designfoils, ioerror
  else:
    print("Read coordinates for " + foilstr + ".")
    foil.setCoordinates(np.array(x), np.array(y))

  # Set zone title for polars

  if (foilstr == 'seed'):
    zonetitle = 'zone t="Seed airfoil polar"'
  else:
    zonetitle = 'zone t="Polars", SOLUTIONTIME=' + str(nextdesign)

  # Read data from polar file (not: negative error code means coordinates were
  # read but not polars)

  cl, cd, ioerror = read_airfoil_data(polarfilename, zonetitle)
  if ( (ioerror == 1) or (ioerror == 2) ):
    print("Warning: polars will not be available for this design.")
    ioerror *= -1
  else:
    print("Read polars for " + foilstr + ".")
    foil.setPolars(np.array(cl), np.array(cd))

  # Copy data to output objects

  if (foilstr == 'seed'): seedfoil = foil
  else: designfoils.append(foil)

  return seedfoil, designfoils, ioerror

################################################################################
# Gets boolean input from user
def get_boolean_input(key, keyval):

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval))
    print("Available choices: True, False\n")
    sel = my_input("Enter new value: ")
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
# Gets color input from user
def get_color_input(key, keyval):

  colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black"]

  validchoice = False
  while (not validchoice):
    print("Current value for " + key + ": " + str(keyval))
    print("Available choices: blue, green, red, cyan, magenta, yellow, black\n")
    sel = my_input("Enter new value: ")

    # Check for valid color

    for c in colors:
      if (sel == c):
        validchoice = True
        retval = sel
        break
    if (not validchoice):
      print("Invalid color specified.  Please enter a valid color.")
      validchoice = False

  return retval

################################################################################
# Options menu: allows user to change plot options
def options_menu():
  global plotoptions

  # Status variable

  options_complete = False
 
  # Print list of plotting options

  print("")
  print("Available plotting options:")
  print("")
  for key in sorted(plotoptions):
    print(key + " [" + str(plotoptions[key]) + "]")
  print("")

  # Get user input

  key = my_input("Enter option to change (or 0 to return): ")

  # True/False settings

  if ( (key == "show_seed_airfoil") or (key == "show_seed_airfoil_only") or
       (key == "show_seed_polar") or (key == "show_seed_polar_only") or
       (key == "save_animation_frames") or (key == "plot_airfoils") or 
       (key == "plot_polars") ):
    options_complete = False
    plotoptions[key] = get_boolean_input(key, plotoptions[key])

  # Change colors

  elif ( (key == "color_for_seed") or (key == "color_for_new_designs") ):
    options_complete = False
    plotoptions[key] = get_color_input(key, plotoptions[key])

  # Change plot bounds

  elif ( (key == "axis_xmax") or (key == "axis_xmin") or
         (key == "axis_ymax") or (key == "axis_ymin") or
         (key == "axis_cdmax") or (key == "axis_cdmin") or
         (key == "axis_clmax") or (key == "axis_clmin") ):
    options_complete = False
    print("Current value for " + key + ": " + str(plotoptions[key]) + '\n')
    sel = my_input("Enter new value for " + key + " or enter 'auto': ")
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
def main_menu(seedfoil, designfoils, prefix, menumode):
  global plotoptions

  exitchoice = False
  while (not exitchoice):
    
    print("")
    print("Options:")
    print("[0] Exit")
    print("[1] Plot a specific design")
    print("[2] Animate all designs")
    print("[3] Monitor an ongoing optimization")
    print("[4] Change plotting options")
    print("")

    choice = my_input("Enter a choice [0-4]: ")

    # Exit design_visualizer

    if (choice == "0"):
      exitchoice = True

    # Plot a single design

    elif (choice == "1"):
      exitchoice = False

      # Turn on matplotlib toolbar

      rcParams['toolbar'] = 'toolbar2'

      # Go to plotting menu

      plotting_complete = False
      while (not plotting_complete): plotting_complete = plotting_menu(
                                                          seedfoil, designfoils)

    # Animate all designs

    elif (choice == "2"):
      exitchoice = False

      # Number of digits in design counter string

      numfoils = len(designfoils)
      width = int(floor(log10(float(numfoils)))) - 1

      # Turn off matplotlib toolbar

      rcParams['toolbar'] = 'None'

      # Loop through designs, updating plot

      for i in range(0, numfoils):
        if (i == 0): init = True
        else: init = False

        if (plotoptions["save_animation_frames"]):

          # Determine number of zeroes to pad with and image file prefix
  
          currwidth = int(floor(log10(float(i+1)))) - 1
          numzeroes = width - currwidth
          imagepref = prefix + numzeroes*'0' + str(i+1)

        else: imagepref = None

        # Update plot

        plot_airfoil(seedfoil, designfoils, i+1, firsttime=init,
                     animation=True, prefix=imagepref)

    # Monitor optimization progress

    elif (choice == "3"):
      exitchoice = False
      print("Monitoring optimization progress. To stop, create a file called " +
            "stop_monitoring.")

      # Turn off matplotlib toolbar and temporarily disable saving images

      rcParams['toolbar'] = 'None'
      temp_save_frames = plotoptions["save_animation_frames"]
      plotoptions["save_animation_frames"] = False

      # Initial flag for plotting (will plot when ioerror <= 0)

      if (seedfoil.npt != 0): ioerror = 0
      else: ioerror = 1

      # Periodically read data and update plot

      init = True
      monitoring = True
      while (monitoring):

        # Update plot

        if (ioerror <= 0): 
          numfoils = len(designfoils)
          plot_airfoil(seedfoil, designfoils, numfoils, firsttime=init,
                       animation=True)
          init = False
        
        # Update airfoil data

        seedfoil, designfoils, ioerror = monitor_airfoil_data(seedfoil,
                                                            designfoils, prefix)

        # Pause for requested update interval

        plt.pause(plotoptions["monitor_update_frequency"])

        # Check for stop_monitoring file

        try:
          f = open('stop_monitoring')
        except IOError:
          continue

        # If stop_monitoring file was found, delete it and break loop
 
        print("stop_monitoring file found. Returning to main menu.")
        f.close()
        remove('stop_monitoring')
        monitoring = False
 
      # Change save_animation_frames back to original setting when done

      plotoptions["save_animation_frames"] = temp_save_frames

    # Change plotting options

    elif (choice == "4"):
      exitchoice = False
      options_complete = False
      while (not options_complete): options_complete = options_menu()

    # Invalid choice

    else:
      print("Error: please enter a choice 0-4.")

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
    prefix = args.prefix
  else:
    coordfilename = 'optfoil_design_coordinates.dat'
    polarfilename = 'optfoil_design_polars.dat'
    prefix = 'optfoil'

  # Read airfoil coordinates and polars

  seedfoil, designfoils, ioerror = load_airfoils_from_file(
                                                   coordfilename, polarfilename)

  # Warning if file is not found

  if (ioerror == 1):
    print("You will not be able to create plots until coordinate data is read.")
    menumode = "options and monitor"
  elif (ioerror < 0):
    print("Only airfoils are available for plotting (no polars).")
    menumode = "airfoils"
  else: menumode = "all"

  # Call main menu

  if (abs(ioerror) <= 1): main_menu(seedfoil, designfoils, prefix, menumode)
