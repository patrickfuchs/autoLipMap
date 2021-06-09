"""This file contains some useful functions for autolipmap

TODO: complete docstrings
"""

def gen_buildH_json_lipid(df, G, dic_PN2MN, DEBUG=False):
    """Takes a pdb dataframe, a graph with mapping names and a dict PN2MN --> generates a buildH json file.
    """
    jsonf = ""
    jsonfMN = ""
    jsonf += "{\n"
    jsonfMN += "{\n"
    # Get residue name.
    resname = df["resname"][0]
    jsonf += f'    "resname": ["{resname}"],\n'
    jsonfMN += f'    "resname": ["{resname}"],\n'
    # Create a list with all carbons names.
    list_Cnodes = [node for node in G.nodes() if node.startswith("C")]
    # Loop over Carbons.
    for Cnode in list_Cnodes:
        nbH = 0
        nbheavy = 0
        for atom in G[Cnode]:
            if atom.startswith("H"):
                nbH += 1
        for atom in G[Cnode]:
            if atom.startswith("C") or atom.startswith("O"):
                nbheavy += 1
        # CH2.
        if nbH == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH2 --> {G[Cnode]}")
            helper1, helper2 = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O")]
            jsonf += f'    "{Cnode}": ["CH2", "{helper1}", "{helper2}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CH2", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}"],\n'
            if DEBUG:
                print()
        # CH3.
        if nbH == 3:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH3 --> {G[Cnode]}")
            (helper1,) = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O")]
            helper2 = [atom for atom in G[helper1] if atom.startswith("C") or atom.startswith("O")][0]
            jsonf += f'    "{Cnode}": ["CH3", "{helper1}", "{helper2}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CH3", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}"],\n'
            if DEBUG:
                print()
        # CH.
        if nbH == 1 and nbheavy == 3:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH --> {G[Cnode]}")
            helper1, helper2, helper3 = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O")]
            jsonf += f'    "{Cnode}": ["CH", "{helper1}", "{helper2}", "{helper3}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CH", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}", "{dic_PN2MN[helper3]}"],\n'
            if DEBUG:
                print()
        # CH double bond.
        if nbH == 1 and nbheavy == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CHdoublebond --> {G[Cnode]}")
            helper1, helper2 = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O")]
            jsonf += f'    "{Cnode}": ["CHdoublebond", "{helper1}", "{helper2}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CHdoublebond", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}"],\n'
            if DEBUG:
                print()       
        # No H reconstruction.
        if nbH == 0:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> no H reconstruction! --> {G[Cnode]}")
                print()
    jsonf += "}\n"
    jsonfMN += "}\n"
    return jsonf, jsonfMN


def gen_buildH_def_file(df, G, basename="cholesterol", DEBUG=False):
    """Takes a pdb dataframe, a graph with mapping names --> generates a buildH def file.
    """
    deff = ""
    # Get residue name.
    resname = df["resname"][0]
    # Create a list with all carbons names.
    list_Cnodes = [node for node in G.nodes() if node.startswith("C")]
    # Loop over Carbons.
    for Cnode in list_Cnodes:
        nbH = 0
        nbheavy = 0
        for atom in G[Cnode]:
            if atom.startswith("H"):
                nbH += 1
        for atom in G[Cnode]:
            if atom.startswith("C") or atom.startswith("O"):
                nbheavy += 1
        # CH2.
        if nbH == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH2 --> {G[Cnode]}")
            carbon_name_woC = Cnode[1:]
            H1_name = f"H{carbon_name_woC}1"
            H2_name = f"H{carbon_name_woC}2"
            deff += f"{basename}{carbon_name_woC}a {resname} {Cnode} {H1_name}\n"
            deff += f"{basename}{carbon_name_woC}b {resname} {Cnode} {H2_name}\n"
            if DEBUG:
                print()
        # CH3.
        if nbH == 3:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH3 --> {G[Cnode]}")
            carbon_name_woC = Cnode[1:]
            H1_name = f"H{carbon_name_woC}1"
            H2_name = f"H{carbon_name_woC}2"
            H3_name = f"H{carbon_name_woC}3"
            deff += f"{basename}{carbon_name_woC}a {resname} {Cnode} {H1_name}\n"
            deff += f"{basename}{carbon_name_woC}b {resname} {Cnode} {H2_name}\n"
            deff += f"{basename}{carbon_name_woC}c {resname} {Cnode} {H3_name}\n"
            if DEBUG:
                print()
        # CH.
        if nbH == 1 and nbheavy == 3:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH --> {G[Cnode]}")
            carbon_name_woC = Cnode[1:]
            H1_name = f"H{carbon_name_woC}1"
            deff += f"{basename}{carbon_name_woC}a {resname} {Cnode} {H1_name}\n"
            if DEBUG:
                print()
        # CH double bond.
        if nbH == 1 and nbheavy == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CHdoublebond --> {G[Cnode]}")
            carbon_name_woC = Cnode[1:]
            H1_name = f"H{carbon_name_woC}1"
            deff += f"{basename}{carbon_name_woC}a {resname} {Cnode} {H1_name}\n"
            if DEBUG:
                print()       
        # No H reconstruction.
        if nbH == 0:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> no H reconstruction! --> {G[Cnode]}")
                print()
    return deff
