import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

from settings import (optimizationsettings, initializationsettings, particleswarmsettings,
                      geneticalgorithmsettings, simplexsettings, xfoilsettings,
                      xfoilpanelingsettings, plotsettings)

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
    #from settings import (optimizationsettings, initializationsettings, particleswarmsettings,
    #                      geneticalgorithmsettings, simplexsettings, xfoilsettings,
    #                      xfoilpanelingsettings, plotsettings)
    root = ET.Element("Settings")
    root.append(optimizationsettings.asXML("OptimizationSettings"))
    root.append(initializationsettings.asXML("InitializationSettings"))
    root.append(particleswarmsettings.asXML("ParticleswarmSettings"))
    root.append(geneticalgorithmsettings.asXML("GeneticAlgorithmSettings"))
    root.append(simplexsettings.asXML("SimplexSettings"))
    root.append(xfoilsettings.asXML("XfoilSettings"))
    root.append(xfoilpanelingsettings.asXML("XfoilPanelingSettings"))
    root.append(plotsettings.asXML("PlotSettings"))

    try:
        f = open(fname, "w")
    except IOError:
        return False
    write_pretty_xml(root, f)
    f.close()

    return True
