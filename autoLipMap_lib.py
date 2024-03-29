"""This file contains some useful functions for autolipmap

TODO: complete docstrings
"""

def gen_buildH_json_lipid(resname, G, dic_PN2MN, DEBUG=False):
    """Takes a resname, a graph with mapping names and a dict PN2MN --> generates a buildH json file.
    """
    jsonf = ""
    jsonfMN = ""
    jsonf += "{\n"
    jsonfMN += "{\n"
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
            if atom.startswith("C") or atom.startswith("O") or atom.startswith("N"):
                nbheavy += 1
        # CH2.
        if nbH == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH2 --> {G[Cnode]}")
            helper1, helper2 = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O") or atom.startswith("N")]
            jsonf += f'    "{Cnode}": ["CH2", "{helper1}", "{helper2}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CH2", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}"],\n'
            if DEBUG:
                print()
        # CH3.
        if nbH == 3:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH3 --> {G[Cnode]}")
            (helper1,) = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O") or atom.startswith("N")]
            helper2 = [atom for atom in G[helper1] if atom.startswith("C") or atom.startswith("O") or atom.startswith("N")][0]
            jsonf += f'    "{Cnode}": ["CH3", "{helper1}", "{helper2}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CH3", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}"],\n'
            if DEBUG:
                print()
        # CH.
        if nbH == 1 and nbheavy == 3:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CH --> {G[Cnode]}")
            helper1, helper2, helper3 = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O") or atom.startswith("N")]
            jsonf += f'    "{Cnode}": ["CH", "{helper1}", "{helper2}", "{helper3}"],\n'
            jsonfMN += f'    "{dic_PN2MN[Cnode]}": ["CH", "{dic_PN2MN[helper1]}", "{dic_PN2MN[helper2]}", "{dic_PN2MN[helper3]}"],\n'
            if DEBUG:
                print()
        # CH double bond.
        if nbH == 1 and nbheavy == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CHdoublebond --> {G[Cnode]}")
            helper1, helper2 = [atom for atom in G[Cnode] if atom.startswith("C") or atom.startswith("O") or atom.startswith("N")]
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


def gen_buildH_def_file(resname, G, basename="cholesterol", DEBUG=False, CHARMMway=False):
    """Takes a resname, a graph with mapping names --> generates a buildH def file.
    """
    deff = ""
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
            if CHARMMway:
                carbon_name_woC = Cnode[1:]
                carbon_firstnumber = carbon_name_woC[0]
                carbon_lastnumbers = carbon_name_woC[1:]
                if DEBUG:
                    print("Cnode", Cnode)
                    print("carbon_name_woC", carbon_name_woC)
                    print("carbon_firstnumber", carbon_firstnumber)
                    print("carbon_lastnumbers", carbon_lastnumbers)
                if carbon_firstnumber == "2":
                    H1_name = f"H{carbon_lastnumbers}R"
                    H2_name = f"H{carbon_lastnumbers}S"
                elif carbon_firstnumber == "3":
                    H1_name = f"H{carbon_lastnumbers}X"
                    H2_name = f"H{carbon_lastnumbers}Y"
                elif carbon_firstnumber == "1":
                    H1_name = f"H{carbon_lastnumbers}A"
                    H2_name = f"H{carbon_lastnumbers}B"
                else:
                    exit("Error in CHARMM naming")
                if DEBUG:
                    print(f"--> Will name H1: {H1_name} and H2: {H2_name}")
            else:
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
            if CHARMMway:
                carbon_name_woC = Cnode[1:]
                carbon_firstnumber = carbon_name_woC[0]
                carbon_lastnumbers = carbon_name_woC[1:]
                if DEBUG:
                    print("Cnode", Cnode)
                    print("carbon_name_woC", carbon_name_woC)
                    print("carbon_firstnumber", carbon_firstnumber)
                    print("carbon_lastnumbers", carbon_lastnumbers)
                if carbon_firstnumber == "2":
                    H1_name = f"H{carbon_lastnumbers}R"
                    H2_name = f"H{carbon_lastnumbers}S"
                    H3_name = f"H{carbon_lastnumbers}T"
                elif carbon_firstnumber == "3":
                    H1_name = f"H{carbon_lastnumbers}X"
                    H2_name = f"H{carbon_lastnumbers}Y"
                    H3_name = f"H{carbon_lastnumbers}Z"
                elif carbon_firstnumber == "1":
                    H1_name = f"H{carbon_lastnumbers}A"
                    H2_name = f"H{carbon_lastnumbers}B"
                    H3_name = f"H{carbon_lastnumbers}C"
                else:
                    exit("Error in CHARMM naming")
                if DEBUG:
                    print(f"--> Will name H1: {H1_name} and H2: {H2_name} and H3: {H3_name}")
            else:
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
            if CHARMMway:
                carbon_name_woC = Cnode[1:]
                carbon_firstnumber = carbon_name_woC[0]
                carbon_lastnumbers = carbon_name_woC[1:]
                if DEBUG:
                    print("Cnode", Cnode)
                    print("carbon_name_woC", carbon_name_woC)
                    print("carbon_firstnumber", carbon_firstnumber)
                    print("carbon_lastnumbers", carbon_lastnumbers)
                if carbon_firstnumber == "2":
                    H1_name = f"H{carbon_lastnumbers}S"
                elif carbon_firstnumber == "3":
                    H1_name = f"H{carbon_lastnumbers}X"
                elif carbon_firstnumber == "1":
                    H1_name = f"H{carbon_lastnumbers}A"
                else:
                    exit("Error in CHARMM naming")
                if DEBUG:
                    print(f"--> Will name H1: {H1_name}")
            else:
                carbon_name_woC = Cnode[1:]
                H1_name = f"H{carbon_name_woC}1"
            deff += f"{basename}{carbon_name_woC}a {resname} {Cnode} {H1_name}\n"
            if DEBUG:
                print()
        # CH double bond.
        if nbH == 1 and nbheavy == 2:
            if DEBUG:
                print(f"{Cnode} has {nbH} Hs, and {nbheavy} heavy atoms ==> CHdoublebond --> {G[Cnode]}")
            if CHARMMway:
                carbon_name_woC = Cnode[1:]
                carbon_firstnumber = carbon_name_woC[0]
                carbon_lastnumbers = carbon_name_woC[1:]
                if DEBUG:
                    print("Cnode", Cnode)
                    print("carbon_name_woC", carbon_name_woC)
                    print("carbon_firstnumber", carbon_firstnumber)
                    print("carbon_lastnumbers", carbon_lastnumbers)
                if carbon_firstnumber == "2":
                    H1_name = f"H{carbon_lastnumbers}R"
                elif carbon_firstnumber == "3":
                    H1_name = f"H{carbon_lastnumbers}X"
                elif carbon_firstnumber == "1":
                    H1_name = f"H{carbon_lastnumbers}A"
                else:
                    exit("Error in CHARMM naming")
                if DEBUG:
                    print(f"--> Will name H1: {H1_name}")
            else:
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
