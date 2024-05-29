"""
author: Joshua Willwerth
this script provides functions to interface with APIs and locally cache data
"""

from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError
from mp_api.client import MPRester as New_MPRester
from urllib.parse import urlencode
import httplib2
from pymatgen.ext.matproj import MPRester as Legacy_MPRester
from emmet.core.thermo import ThermoType
from pymatgen.entries.mixing_scheme import MaterialsProjectDFTMixingScheme
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element

import os
import json
import time
import numpy as np

# input API keys here
MPDS_api_key = "29r3PRLsIW8YNzxxTCVgDmxP0ylZqfWSUvBVE3IVyrZvpqKG"
#New_MP_api_key = "Jcw46im7UV1xOfHzbZZ8nkq8BH00Pf6s"
New_MP_api_key ="S0l8OMrjeiqxN9y4nXNgzvVcAtCmbKWF"
Legacy_MP_api_key = "2d5wyVmhDCpPMAkq"

data_dir = "70x70_matrix_data"
enthalpies_file = data_dir + "\\fusion_enthalpies.json"
melt_temps_file = data_dir + "\\fusion_temperatures.json"

if not os.path.exists(enthalpies_file):
    enthalpies = {}
else:
    with open(enthalpies_file, "r") as file:
        enthalpies = json.load(file)

if not os.path.exists(melt_temps_file):
    melt_temps = {}
else:
    with open(melt_temps_file, "r") as file:
        melt_temps = json.load(file)


def t_at_boundary(t, boundary):
    return t <= boundary[0] + 2 or t >= boundary[1] - 2


def section_liquidus(points):
    sections = []
    current_section = []

    x1, y1 = points[0]
    x2, y2 = points[1]

    if x2 > x1:
        direction = "increasing"
        current_section.append(points[0])
    elif x2 < x1:
        direction = "decreasing"
        current_section.append(points[0])
    else:
        direction = None
        sections.append([points[0]])

    for i in range(1, len(points) - 1):
        x1, y1 = points[i]
        x2, y2 = points[i + 1]

        if x2 > x1:
            new_direction = "increasing"
        elif x2 < x1:
            new_direction = "decreasing"
        else:
            new_direction = None
        # add x1, y1
        current_section.append(points[i])

        if new_direction != direction or new_direction is None:
            if current_section:
                sections.append(current_section)
                current_section = []
        direction = new_direction

    if current_section:
        current_section.append(points[-1])
        sections.append(current_section)
    return sections


def within_tol_from_line(p1, p2, p3, tol):
    try:
        m = (p1[1] - p2[1]) / (p1[0] - p2[0])
    except ZeroDivisionError:
        return tol >= abs(p1[1] - p2[1])
    y_h = m * (p3[0] - p1[0]) + p1[1]
    return tol >= abs(p3[1] - y_h)


def fill_liquidus(p1, p2, max_interval):
    num_between = int(np.floor((p2[0] - p1[0]) / max_interval))
    filled_X = np.linspace(p1[0], p2[0], num_between + 2)
    filled_T = np.linspace(p1[1], p2[1], num_between + 2)
    filled_section = [[filled_X[i], filled_T[i]] for i in range(len(filled_X))]
    return filled_section[1:-1]


# load the MPDS data for system from the API or local cache
# returns dict, dict, 2D list
# ind is a testing feature used to select a phase diagram at a specific rank. highest ranked PD selected by default
def get_MPDS_data(components, pd=0, dtype='all'):
    """Retrive MPDS data from cache or API.
    Required: List of components. Must provide only two components if querying binary phase diagrams.
    pd -> index of PD sorted by selection criteria. Specifying 'None' will cache all PDs and return nothing.
    dtype -> type of data returned by function. 'comp' or 'all_thermo' will return only component data,
    with 'all_thermo' returning additional component data. Default return for 'all' is a tuple of three fields:
    MPDS JSON (dict), Component data (dict), MPDS Liquidus (2D list)"""
    client = MPDSDataRetrieval(api_key=MPDS_api_key)
    client.dtype = MPDSDataTypes.PEER_REVIEWED
    sys = '-'.join(sorted(components))
    component_data = {}

    for comp in components:
        if comp not in enthalpies or comp not in melt_temps:
            print("searching for " + comp + " data from MPDS...")

            # component properties (sample Al)
            #         "property": {
            #           "name": "enthalpy change at melting point",
            #           "units": "kJ g-at.-1",
            #           "domain": "thermal and thermodynamic properties",
            #           "scalar": 10.71,
            #           "category": "enthalpy change at phase transition"
            #         },
            #         "condition": [
            #           {
            #             "name": "Temperature",
            #             "units": "K",
            #             "scalar": 933.5

            comp_fields = {
                'P': ['sample.measurement[0].property.scalar', 'sample.measurement[0].condition[0].scalar']}
            data = client.get_data(search={'formulae': comp,
                                           'props': 'enthalpy change at melting point'}, fields=comp_fields)

            # either take the first value or median value
            if comp not in enthalpies:
                print("caching enthalpy of fusion data")
                enthalpies[comp] = float(data[0][0]) * 1000
                with open(enthalpies_file, "w") as f:
                    json.dump(enthalpies, f)
            if comp not in melt_temps:
                print("caching melting temperature data")
                melt_temps[comp] = float(data[0][1])
                with open(melt_temps_file, "w") as f:
                    json.dump(melt_temps, f)

            print(comp + ": enthalpy of fusion = " + "{:,.0f}".format(enthalpies[comp])
                  + " J/mol at " + "{:,.0f}".format(melt_temps[comp]) + " K")

        component_data[comp] = [enthalpies[comp], melt_temps[comp]]

        if dtype == 'all_thermo':
            comp_fields = {
                'P': ['sample.measurement[0].property.scalar', 'sample.measurement[0].condition[0].scalar']}
            data = client.get_data(search={'formulae': comp,
                                           'props': 'entropy'}, fields=comp_fields)

            best_entry = [float(data[0][0]), float(data[0][1])]
            for d in data:
                if float(d[1]) < 300 and float(d[0]) > best_entry[0]:
                    best_entry = [float(d[0]), float(d[1])]

            component_data[comp].extend(best_entry)

    # will skip system data if only component data is desired (comp) or additional data as specified above (all_thermo)
    if dtype == 'comp' or dtype == 'all_thermo':
        return component_data

    if isinstance(pd, str):
        ind = 0

        # look for cached json for the system
        sys_file = os.path.join(f"{data_dir}\\{sys}", f"{sys}_MPDS_PD_{ind}.json")
        while os.path.exists(sys_file):
            with open(sys_file, 'r') as f:
                MPDS_json = json.load(f)
                if MPDS_json['entry'] == pd:
                    print("\nloading JSON from cache...")
                    return MPDS_json, component_data, extract_MPDS_liquidus(MPDS_json, components)
                else:
                    ind += 1
                    sys_file = os.path.join(f"{data_dir}\\{sys}", f"{sys}_MPDS_PD_{ind}.json")

        print("\nspecified entry not found in cache, please download all entries first by specifying None as the index")
        return {}, {}, []

    elif isinstance(pd, int) or pd is None:
        ind = pd

        # look for cached json for the system
        sys_file = os.path.join(f"{data_dir}\\{sys}", f"{sys}_MPDS_PD_{ind}.json")
        if os.path.exists(sys_file):
            # get MPDS data from stored jsons in liquidus curves folder
            print("\nloading JSON from cache...")
            with open(sys_file, 'r') as f:
                MPDS_json = json.load(f)

            return MPDS_json, component_data, extract_MPDS_liquidus(MPDS_json, components)

        else:
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
            valid_JSONs = []

            # search for phase diagrams and filter out empty entries
            try:
                diagrams = [d for d in client.get_data(
                    search={'elements': sys, 'classes': 'binary'}, fields=sys_fields) if d]

                # step 1: determine which phase diagrams have usable liquidus data
                for d in diagrams:
                    dia_json = {}
                    for i in range(len(sys_fields['C'])):
                        dia_json[sys_fields['C'][i]] = d[i]

                    # remove any diagrams that don't span the entire composition range
                    if dia_json['comp_range'] != [0, 100]:
                        continue

                    # if liquidus curve is extracted from the diagram json without issue, assign score and add to list
                    dia_liquidus = extract_MPDS_liquidus(dia_json, verbose=False)
                    if dia_liquidus:
                        valid_JSONs.append(dia_json)

            except APIError:
                print(" Got 0 hits")

            if ind is not None and len(valid_JSONs) <= ind:
                print("specified index " + str(ind) + " exceeds number of valid JSONs found for the " + sys + " system")
                MPDS_json = {'reference': None}
                if ind == 0:
                    print("caching empty JSON...\n")
                    with open(sys_file, "w") as f:
                        json.dump(MPDS_json, f)
                return MPDS_json, component_data, extract_MPDS_liquidus(MPDS_json, verbose=False)

            # step 2: query additional information on valid PDs
            endpoint = "https://api.mpds.io/v0/search/facet"
            req = httplib2.Http()

            for dia_json in valid_JSONs:

                url = endpoint + '?' + urlencode(
                    {
                        'q': json.dumps({'entry': dia_json['entry']}),
                        'pagesize': 10,
                        'dtype': 1
                    })

                response, content = req.request(
                    uri=url,
                    method='GET',
                    headers={'Key': MPDS_api_key}
                )

                time.sleep(2)

                try:
                    info = json.loads(content)
                except json.decoder.JSONDecodeError as e:
                    print(e)
                    continue

                components = dia_json['chemical_elements']
                dia_json['year'] = int(info['out'][0][6])
                dia_json['unlocked'] = bool(info['out'][0][4]) or False
                dia_json['pdtype'] = info['out'][0][2].split('.')[0]
                dia_liquidus = extract_MPDS_liquidus(dia_json, verbose=False)
                dia_json['tm_agreement'] = [dia_liquidus[0][1] - melt_temps[components[0]],
                                            dia_liquidus[-1][1] - melt_temps[components[1]]]
                dia_json['num_solid_phases'] = 0
                for shape in dia_json['shapes']:
                    if shape['is_solid'] and shape['nphases'] == 1:
                        dia_json['num_solid_phases'] += 1


            # step 3: select the 'best' phase diagram from valid diagrams)

            # sort by 'unlocked' then by year, pdtype (preferred exp over comp), tm_agreement (preferred smallest diff)
            valid_JSONs.sort(
                key=lambda x: (not x['unlocked'], -x['year'], x['pdtype'], sum(map(abs, x['tm_agreement']))))

            # pd ind not specified or specified as a valid index -> caches PD at specified index and returns pd info
            if ind is not None:
                MPDS_json = valid_JSONs[ind]
                print(f"caching MPDS liquidus from entry at {MPDS_json['reference']['entry']} as {sys_file}...\n")
                with open(sys_file, "w") as f:
                    json.dump(MPDS_json, f)
                return MPDS_json, component_data, extract_MPDS_liquidus(MPDS_json, verbose=False)
            # pd ind specified as None -> caches all PDs and doesn't return anything
            else:
                ind = 0
                for MPDS_json in valid_JSONs:
                    sys_file = os.path.join(f"{data_dir}\\{sys}", f"{sys}_MPDS_PD_{ind}.json")
                    print(f"caching MPDS liquidus from entry at {MPDS_json['reference']['entry']} as {sys_file}...\n")
                    with open(sys_file, "w") as f:
                        json.dump(MPDS_json, f)
                    ind += 1
                return {}, {}, []


# pull liquidus curve data from MPDS json and convert to list of [X, T] coordinates in ascending composition order
# returns 2D list
def extract_MPDS_liquidus(MPDS_json, verbose=True):
    if MPDS_json['reference'] is None:
        if verbose:
            print("system JSON does not contain any data!\n")
        return None

    components = MPDS_json['chemical_elements']
    if verbose:
        print("reading MPDS liquidus from entry at " + MPDS_json['reference']['entry'] + "...\n")

    # extract liquidus curve svgpath from system JSON
    data = ""
    for boundary in MPDS_json['shapes']:
        if 'label' in boundary and boundary['label'] == 'L':
            data = boundary['svgpath']
            break
    if not data:
        if verbose:
            print("no liquidus data found in JSON!")
        return None

    # split svgpath into tags, ordered pairs
    data = data.split(' ')

    # remove 'L' and 'M' tags, so only ordered pairs remain
    data = [s for s in data if not (s == 'L' or s == 'M')]

    # remove points at the edge of the graph boundaries
    data = [s for s in data if not t_at_boundary(float(s.split(',')[1]), MPDS_json['temp'])]

    # convert string pairs into [X, T] float pairs and store as liquidus
    X = [float(i.split(',')[0]) / 100.0 for i in data]
    T = [float(i.split(',')[1]) + 273.15 for i in data]
    MPDS_liquidus = [[X[i], T[i]] for i in range(len(X))]

    # filter out duplicate values in the liquidus curve; greatly improves runtime efficiency
    for i in reversed(range(len(MPDS_liquidus) - 1)):
        if MPDS_liquidus[i][0] == 0 or MPDS_liquidus[i][1] == 0:
            continue
        if abs(1 - MPDS_liquidus[i + 1][0] / MPDS_liquidus[i][0]) < 0.0005 and \
                abs(1 - MPDS_liquidus[i + 1][1] / MPDS_liquidus[i][1]) < 0.0005:
            del (MPDS_liquidus[i + 1])

    if len(MPDS_liquidus) < 3:
        if verbose:
            print("MPDS liquidus does not span the entire composition range!")
        return None

    # split liquidus into segments of continuous points
    sections = section_liquidus(MPDS_liquidus)

    # sort sections by descending size
    sections.sort(key=len, reverse=True)

    # sort each section by ascending composition
    for section in sections:
        section.sort()

    # record endpoints of main section
    MPDS_liquidus = sections.pop(0)

    lhs = [0, melt_temps[components[0]]]
    rhs = [1, melt_temps[components[1]]]
    # print(components)

    # append sections to the liquidus if not overlapping in range
    for section in sections:

        # if section upper bound is less than main section lower bound
        if section[-1][0] <= MPDS_liquidus[0][0] and within_tol_from_line(MPDS_liquidus[0], lhs, section[-1], 100):
            MPDS_liquidus = section + MPDS_liquidus

        # if section lower bound is greater than main section upper bound
        elif section[0][0] >= MPDS_liquidus[-1][0] and within_tol_from_line(MPDS_liquidus[-1], rhs, section[0], 100):
            MPDS_liquidus.extend(section)

        # i'll admit it at this point I am feeling pretty dumb about not using a svgpath parser because I have
        # do all of these strange exceptions to make this work. don't worry about this, it's just some edge case
        elif len(section) == 2:
            if section[0][0] < MPDS_liquidus[0][0] and within_tol_from_line(MPDS_liquidus[0], lhs, section[0], 170):
                MPDS_liquidus = [section[0]] + MPDS_liquidus

            elif section[-1][0] > MPDS_liquidus[-1][0] and within_tol_from_line(MPDS_liquidus[-1], rhs, section[-1],
                                                                                170):
                MPDS_liquidus.extend(section)

    MPDS_liquidus.sort()

    # if the liquidus does not have endpoints near the ends of the composition range, melting temps won't be good
    if 100 * MPDS_liquidus[0][0] > 3 or 100 * MPDS_liquidus[-1][0] < 97:
        if verbose:
            print(f"MPDS liquidus does not span the entire composition range! "
                  f"({100 * MPDS_liquidus[0][0]}-{100 * MPDS_liquidus[-1][0]})")
        return None

    # fill in ranges with missing points
    for i in reversed(range(len(MPDS_liquidus) - 1)):
        if MPDS_liquidus[i + 1][0] - MPDS_liquidus[i][0] > 0.06:
            filler = fill_liquidus(MPDS_liquidus[i], MPDS_liquidus[i + 1], 0.03)
            for point in reversed(filler):
                MPDS_liquidus.insert(i + 1, point)

    return MPDS_liquidus


# dft_types = ["GGA", "GGA/GGA+U", "R2SCAN", "GGA/GGA+U/R2SCAN"]


# returns the DFT convex hull of a given system with specified functionals
def get_dft_convexhull(components, verbose=False):
    dft_type = "GGA/GGA+U"
    # if dft_type not in dft_types:
    #     print("invalid DFT type")
    #     exit(1)
    if 'Yb' in components:
        dft_type = "GGA"

    if verbose:
        print("using DFT entries solved with", dft_type, "functionals")
    sys = '-'.join(sorted(components))

    # dft_type_path = '_'.join(dft_type.split('/'))
    # dft_entries_file = os.path.join(f"{data_dir}\\{sys}", f"{sys}_MP_ENTRIES_{dft_type_path}.json")
    dft_entries_file = os.path.join(f"{data_dir}\\{sys}", f"{sys}_ENTRIES_MP_GGA.json")

    if os.path.exists(dft_entries_file):
        with open(dft_entries_file, "r") as f:
            dft_entries = json.load(f)

        try:
            pd = PhaseDiagram(elements=[Element(c) for c in components],
                              entries=[ComputedEntry.from_dict(e) for e in dft_entries])
            if verbose:
                print(len(pd.stable_entries) - 2, "stable line compound(s) on the DFT convex hull\n")
            return pd
        except ValueError as e:
            print(f"error loading DFT entries from cache: {e}")

    # no cache or invalid cached data
    entries = []

    # using legacy MP energies (GGA)
    if dft_type == "GGA":
        with Legacy_MPRester(Legacy_MP_api_key) as MPR:
            entries = MPR.get_entries_in_chemsys(components, inc_structure=True)

    # using new MP energies (GGA/GGA+U, R2SCAN, GGA/GGA+U/R2SCAN)
    else:
        with New_MPRester(New_MP_api_key) as MPR:
            # if dft_type == "R2SCAN" or dft_type == "GGA/GGA+U/R2SCAN":
            #     scan_entries = MPR.get_entries_in_chemsys(components,
            #                                               additional_criteria={
            #                                                   'thermo_types': [ThermoType.R2SCAN]})
            if dft_type == "GGA/GGA+U" or dft_type == "GGA/GGA+U/R2SCAN":
                gga_entries = MPR.get_entries_in_chemsys(components,
                                                         additional_criteria={
                                                             'thermo_types': [ThermoType.GGA_GGA_U]})

        # if dft_type == "GGA/GGA+U/R2SCAN":
        #     entries = MaterialsProjectDFTMixingScheme().process_entries(scan_entries + gga_entries,
        #                                                                 verbose=verbose)
        if dft_type == "GGA/GGA+U":
            entries = MaterialsProjectDFTMixingScheme().process_entries(gga_entries, verbose=verbose)
        # elif dft_type == "R2SCAN":
        #     entries = MaterialsProjectDFTMixingScheme().process_entries(scan_entries, verbose=verbose)

    if verbose:
        print(f"caching DFT entry data as {dft_entries_file}...")
    dft_entries = [e.as_dict() for e in entries]
    for e in dft_entries:
        e.pop('structure')
        e.pop('data')
    with open(dft_entries_file, "w") as f:
        json.dump(dft_entries, f)

    try:
        pd = PhaseDiagram(elements=[Element(c) for c in components],
                          entries=[ComputedEntry.from_dict(e) for e in dft_entries])
        if verbose:
            print(len(pd.stable_entries) - 2, "stable line compound(s) on the DFT convex hull\n")
        return pd
    except ValueError as e:
        print(f"error with DFT entries downloaded from API: {e}")
        return None


def identify_MPDS_phases(MPDS_json, verbose=False):
    if MPDS_json['reference'] is None:
        if verbose:
            print("system JSON does not contain any data!\n")
        return None

    misc_gap_labels = []
    for label in MPDS_json['labels']:
        delim_label = label[0].split(' ')
        if len(delim_label) == 3 and delim_label[0][0] == 'L' and delim_label[2][0] == 'L':
            misc_gap_labels.append([label[1][0] / 100.0, label[1][1] + 273.15])

    phases = []
    data = ""
    for shape in MPDS_json['shapes']:

        if 'nphases' in shape and 'is_solid' in shape:
            # indentify line compounds and single-phase solid solutions
            if shape['nphases'] == 1 and shape['is_solid'] and 'label' in shape:
                if '(' in shape['label'].split(' ')[0]:
                    continue
                # split svgpath into tags, ordered pairs
                data = shape['svgpath'].split(' ')
                # remove 'L' and 'M' tags, so only ordered pairs remain
                data = [s for s in data if not (s == 'L' or s == 'M')]
                if not data:
                    if verbose:
                        print(f"no point data found for phase {shape['label']} in JSON!")
                    continue
                # convert to coordinates in fractional composition, kelvin
                data = [[float(i.split(',')[0]) / 100.0, float(i.split(',')[1]) + 273.15] for i in data]
                # sort by ascending X value and filter out components with solid solubility
                data.sort(key=lambda x: x[0])
                if data[0][0] < 0.01 or data[-1][0] > 0.99:
                    continue
                # sort by ascending T value
                data.sort(key=lambda x: x[1])
                tbounds = [data[0], data[-1]]

                if shape['kind'] == 'phase':
                    data.sort(key=lambda x: x[0])
                    cbounds = [data[0], data[-1]]

                    phases.append({'name': shape['label'].split(' ')[0], 'comp': tbounds[1][0],
                                   'cbounds': cbounds, 'tbounds': tbounds, 'type': 'ss'})
                else:  # kind == compound
                    phases.append({'name': shape['label'].split(' ')[0], 'comp': tbounds[1][0],
                                   'tbounds': tbounds, 'type': 'lc'})

            # identify miscibility gap regions, can be erraneous when the phase diagram is underconstrained
            if shape['nphases'] == 2 and misc_gap_labels:
                # for each label, check each boundary to see if composition endpoints are on opposite sides,
                # and the highest temp is above the label
                for mgl in misc_gap_labels:
                    # split svgpath into tags, ordered pairs
                    data = shape['svgpath'].split(' ')
                    # remove 'L' and 'M' tags, so only ordered pairs remain
                    data = [s for s in data if not (s == 'L' or s == 'M')]
                    if not data:
                        continue
                    # convert to coordinates in fractional composition, kelvin
                    data = [[float(i.split(',')[0]) / 100.0, float(i.split(',')[1]) + 273.15] for i in data]
                    # sort by ascending X value
                    data.sort(key=lambda x: x[0])
                    if not data[0][0] < mgl[0] < data[-1][0]:
                        continue
                    cbounds = [data[0], data[-1]]
                    # sort by ascending T value
                    data.sort(key=lambda x: x[1])
                    if not data[-1][1] > mgl[1]:
                        continue
                    tbounds = [data[0], data[-1]]

                    phases.append({'comp': tbounds[1][0], 'cbounds': cbounds,
                                   'tbounds': tbounds, 'type': 'mig'})
    if not data:
        if verbose:
            print("no phase data found in JSON!")
        return phases

    phases.sort(key=lambda x: x['comp'])
    return phases
