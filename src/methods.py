import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import sys

from settings import (optimizationsettings, initializationsettings, particleswarmsettings,
                      geneticalgorithmsettings, simplexsettings, xfoilsettings,
                      xfoilpanelingsettings, plotsettings)
from operatingpoints import operatingpoints
from constraints import constraints
from data import data

def write_pretty_xml(elem, f, indentlevel=0, indent='  ', header=True):
    '''Writes an XML element to file with nice indentation. xml.dom.minidom is used to do the
       writing, because etree.ElementTree dumps it all out on a single line without indentation.

    Inputs:
        elem: an ElementTree XML element
        f: an open file object for writing
        indentlevel: how many tab levels the entire element should be indented
        indent: indent character(s)

    Reference:
    https://stackoverflow.com/questions/17402323/use-xml-etree-elementtree-to-print-nicely-formatted-xml-files
    '''
    ugly = ET.tostring(elem)
    pretty = minidom.parseString(ugly).toprettyxml(indent=indent)
    prettylist = pretty.split('\n')
    for line in prettylist:
        # Skip header if requested and blank lines
        if not header and line.startswith("<?xml version"):
            continue
        if len(line) > 0:
            f.write(indentlevel*indent + line + "\n")

def save_settings(fname):
    '''Saves settings, operating points, constraints, and seed airfoil source to XML file

    Inputs:
        fname: name or path of xml file to save
    '''

    # Ensure filename ends in .xml
    if not fname.endswith(".xml"):
        fname += ".xml"

    root = ET.Element("XoptfoilCaseSettings")
    settingselem = ET.SubElement(root, "Settings")
    settingselem.append(optimizationsettings.asXML("OptimizationSettings"))
    settingselem.append(initializationsettings.asXML("InitializationSettings"))
    settingselem.append(particleswarmsettings.asXML("ParticleswarmSettings"))
    settingselem.append(geneticalgorithmsettings.asXML("GeneticAlgorithmSettings"))
    settingselem.append(simplexsettings.asXML("SimplexSettings"))
    settingselem.append(xfoilsettings.asXML("XfoilSettings"))
    settingselem.append(xfoilpanelingsettings.asXML("XfoilPanelingSettings"))
    settingselem.append(plotsettings.asXML("PlotSettings"))
    root.append(operatingpoints.asXML("OperatingPoints"))
    root.append(constraints.asXML("Constraints"))
    foilelem = data.seed_airfoil.sourceAsXML("SeedAirfoil")
    if foilelem is not None:
        root.append(foilelem)

    try:
        f = open(fname, "w")
    except IOError:
        return False
    write_pretty_xml(root, f)
    f.close()

    return True

def read_settings(fname):
    '''Reads settings, operating points, constraints, and seed airfoil source from XML file

    Inputs:
        fname: name or path of xml file

    Returns:
        0 on success, 1 if file not found, 2 if I/O error, 3 if parse error
        msg: success/error message
    '''

    # Open file and catch errors
    msg = "Success"
    try:
        tree = ET.parse(fname)
    except FileNotFoundError:
        msg = "Error parsing {:s}: file not found.".format(fname)
        sys.stderr.write(msg+"\n")
        return 1, msg
    except IOError:
        msg = "Error parsing {:s}: I/O error.".format(fname)
        sys.stderr.write(msg+"\n")
        return 2, msg
    except ET.ParseError as e:
        msg = "Error parsing {:s}: {:s}.".format(fname, e.msg)
        sys.stderr.write(msg+"\n")
        return 3, msg

    # Read settings
    root = tree.getroot()
    for child in root:
        if child.tag == "Settings":
            for grandchild in child:
                if grandchild.tag == "OptimizationSettings":
                    optimizationsettings.fromXML(grandchild)
                elif grandchild.tag == "InitializationSettings":
                    initializationsettings.fromXML(grandchild)
                elif grandchild.tag == "ParticleswarmSettings":
                    particleswarmsettings.fromXML(grandchild)
                elif grandchild.tag == "GeneticAlgorithmSettings":
                    geneticalgorithmsettings.fromXML(grandchild)
                elif grandchild.tag == "SimplexSettings":
                    simplexsettings.fromXML(grandchild)
                elif grandchild.tag == "XfoilSettings":
                    xfoilsettings.fromXML(grandchild)
                elif grandchild.tag == "XfoilPanelingSettings":
                    xfoilpanelingsettings.fromXML(grandchild)
                elif grandchild.tag == "PlotSettings":
                    plotsettings.fromXML(grandchild)
                else:
                    sys.stderr.write("Warning: unrecognized element {:s}.\n"\
                                     .format(grandchild.tag))
        elif child.tag == "OperatingPoints":
            operatingpoints.fromXML(child)

    return 0, msg
