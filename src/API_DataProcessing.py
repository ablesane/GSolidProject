"""
author: Andria Lesane
This script provides functions to interface with MPDS API. Contains preliminary data processing and visualization functions.
"""
#MPDS API
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError
from mp_api.client import MPRester as New_MPRester
from urllib.parse import urlencode

#Materials Project API
from mp_api.client import MPRester as New_MPRester
from urllib.parse import urlencode
import httplib2
from pymatgen.ext.matproj import MPRester as Legacy_MPRester
from emmet.core.thermo import ThermoType
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element

#Data Visualization
%matplotlib inline
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

#Data Conversion
from IPython.display import HTML, SVG
import numpy as np
from svgpath2mpl import parse_path
from operator import itemgetter
import os
import time
import json

#Solving Equations
from sympy import Eq, Symbol, solve

MPDS_api_key = "SECRETKey"

# Below is the function used to query all available phase diagram data that contains the list of specified chemical components. 
# Further development of function will include methods to extract the most appropriate phase diagram of the specified binary system containing solid solutions.

def get_MPDS_data(components):
    """Retrive MPDS data from cache or API.
    Required: List of components. Must provide only two components if querying binary phase diagrams."""
    client = MPDSDataRetrieval(api_key=MPDS_api_key)
    client.dtype = MPDSDataTypes.PEER_REVIEWED
    sys = '-'.join(sorted(components))
    component_data = {}

    print("\nsearching for phase diagram JSON from MPDS...")
    # phase diagram properties
    # arity - num elements - 2
    # naxes - num axis - 2
    # diatype - phase diagram type - "binary"
    # comp_range - 0-100
    # reference - link to entry
    # shapes - phase boundary info
    # chemical_elements - alphabetized chemcial elements in system
    # temp - temp range of diagram

    sys_fields = {'C': ['chemical_elements', 'entry', 'comp_range', 'temp', 'labels', 'shapes', 'reference']}
    
    diagrams = [d for d in client.get_data(
                search={'elements': sys, 'classes': 'binary'}, fields=sys_fields)if d]
    
    valid_JSON=[]
    
    for d in diagrams:
                    dia_json = {}
                    for i in range(len(sys_fields['C'])):
                        dia_json[sys_fields['C'][i]] = d[i]
                        valid_JSON.append(dia_json)

    return valid_JSON





#The function below uses the specified material system to extract solidus curves into matplotlib compatible path. 
# Accepts phase diagram data in JSON format. Function uses the same error handling developed for the extract_MPDS_liquidus function.

def extract_MPDS_Solidus(MPDS_json, verbose=True):
    if MPDS_json['reference'] is None:
        if verbose:
            print("system JSON does not contain any data!\n")
        return None

    components = MPDS_json['chemical_elements']
    if verbose:
        print("reading MPDS solidus from entry at " + MPDS_json['reference']['entry'] + "...\n")
    
    # extract solidus curve svgpath from system JSON
    curves = ""
    for boundary in MPDS_json['shapes']:
        if 'label' in boundary and boundary['label'] != 'L' and boundary['label'] != 'G':
            curves = curves + boundary['svgpath']
            continue
    if not MPDS_json:
        if verbose:
            print("no solidus data found in JSON!")
        return None
    
    path = parse_path(curves)

    
    # return matplotlib path
    return path, path.vertices, path.codes


#Function to plot phase diagram, accepts json file and svgpath

def plot_phase_diagram(data, svg_path):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    patch = patches.PathPatch(svg_path)
    ax.add_patch(patch)
    ax.set_ylim(data['temp'][0],data['temp'][1])
    ax.set_xlim(data['comp_range'][0],data['comp_range'][1])
    
    return plt.show()


#The function below pulls the temperature and fractional composition at the rightmost point of each phase. 

def extract_solubility_limit(path):
    # convert to coordinates in fractional composition, kelvin
    vertices = [[float(i[0])/100 , float(i[1])+273.15] for i in path._vertices]
    codes = path._codes

    sections=[]
    ind = 0
    solubility_limits = []

    for i in range(0,len(vertices)-1):
        if (codes[i+1] == 1 or None):
            b = i + 1
            sections.append(vertices[ind:b])
            ind = b
        else: 
            continue

    sections.append(vertices[ind:])

    
    for section in sections:
        # Sort by X (descending), determinding the rightmost point of the given svgpath
        sorted_coordinates= sorted(section, key=lambda c: ( -c[0]))
        solubility_limits.append(sorted_coordinates[0])

    return solubility_limits


#Equations for solving for omega

#Formation energy is stored in units of eV/atom in Materials Project, creating conversion function to put in correct units.

def convert_formation_energy(eV_atom):
    J_mol= eV_atom*(1.6*10**-19)*(6.02*10**23)

    return J_mol

def calculate_gibbs_free_energy(Tm,molEntr,T,transEnth):
    #Tm is the phase transition temperature (units:  Kelvin)
    # molEntr is Molar Entropy (units: J/mol/K)
    # T is Temperature (units: Kelvin)
    #transEnth is transition enthalpy (units: J/mol)

    
    deltaS = transEnth/Tm
    G_CS1 = (-T)*(molEntr)
    deltaG = transEnth -(T*deltaS)
    G_CS2 = G_CS1 + deltaG

    return G_CS2

def solve_mixing_enthalpy(G, x, y, T, formEng, R=8.314):
    # G is Gibbs Free Energy (units: J/mol)
    # x is solubility limit (units:fractional composition)
    # y is the composition of competing oxide (units:fractional composition)
    # T is observed temperature (units: Kelvin)
    # formEng is Formation energy (units: J/mol)
    # R is ideal gas constant (units: J/mol_K)

    w = Symbol('w')

    eqn = (G + w*x*(1-x)+R*T*(x*np.log(x)+(1-x)*np.log(1-x))) + (y-x) *(w-2*w*x+R*T*(np.log(x)-np.log(1-x))) - formEng

    return solve(eqn)

