import sys

import argparse

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import lipids_info

# This is the distance below which a bond is defined between any 2 atoms.
MAXBOND = 1.7 # in A


def pdb2pandasdf(pdb_filename):
    """Reads a PDB file and returns a pandas data frame.

    Arguments
    ---------
    pdb_filename : string

    Returns
    -------
    pandas dataframe
        The col index are: atnum, atname, resname, resnum, x, y, z
    """
    rows = []
    with open(pdb_filename, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atnum = int(line[6:11])
                atname = line[12:16].strip()
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                rows.append((atnum, atname, resname, resnum, x, y, z))
    df_atoms = pd.DataFrame(rows, columns=["atnum", "atname", "resname",
                                           "resnum", "x", "y", "z"])
    # Warning if duplicate atnames.
    list_atnames = list(df_atoms["atname"])
    for atname in list_atnames:
        n = list_atnames.count(atname)
        if  n > 1:
            print(f"!!! Warning, {atname} appears {n} times in pdb file !!!")
    return df_atoms


def calc_distance(x1, y1, z1, x2, y2, z2):
    return ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**.5


def build_graph_from_pdb(df, theoretical_nb_nodes, theoretical_nb_edges):
    """Reads a pdb and builds a graph of the lipid molecule.

    Arguments
    ---------
    df : pandas dataframe
        Contains the pdb lines organized in columns "atnum", "atname",
        "resname", "resnum", "x", "y", "z"
    theoretical_nb_nodes : int
        Number of nodes of the graph built from mapping names.
    theoretical_nb_edges : int
        Number of edges of the graph built from mapping names.

    Returns
    -------
    networkx graph instance
        The graph describing the molecule found in the pdb where each node is 
        an atom (with label atom name) and each edge is a bond.
    """
    G = nx.Graph()
    # Add nodes.
    for atname in df["atname"]:
        # Check each atname is unique.
        if atname in G.nodes():
            exit(f"{atname} already exists in lipid read from PDB. "
                 "Can't build a graph with 2 identical node names.")
        else:
            G.add_node(atname)
    # Check nb of nodes OK.
    if G.number_of_nodes() != theoretical_nb_nodes:
        #raise(Exception, f"Actual POPC has {G.number_of_nodes()} atoms, while "
        #      "it should have {theoretical_nb_nodes}")
        exit(f"Actual POPC has {G.number_of_nodes()} atoms, while it should "
             f"have {theoretical_nb_nodes}")
    # Now create edges for each chemical bond.
    for i in range(G.number_of_nodes()-1):
        for j in range(i+1, G.number_of_nodes()):
            namei = df.iloc[i]["atname"]
            namej = df.iloc[j]["atname"]
            numi = df.iloc[i]["atnum"]
            numj = df.iloc[j]["atnum"]
            # If both atoms are H, skip edge creation.
            if "H" in namei and "H" in namej:
                continue
            x1, y1, z1 = df.iloc[i][["x", "y", "z"]]
            x2, y2, z2 = df.iloc[j][["x", "y", "z"]]
            #print(i, j, calc_distance(x1, y1, z1, x2, y2, z2))
            if calc_distance(x1, y1, z1, x2, y2, z2) < MAXBOND:
                G.add_edge(namei, namej)
                print(f"Add bond #{G.number_of_edges()} between atom {numi} "
                      f"({namei}) and atom {numj} ({namej})")
    # Check nb of edges OK.
    if G.number_of_edges() != theoretical_nb_edges:
        exit(f"Actual POPC has {G.number_of_edges()} bonds, while it should "
             f"have {theoretical_nb_edges}")
    print(f"PDB graph has {G.number_of_nodes()} nodes (atoms) and "
          f"{G.number_of_edges()} edges (bonds)")
    return G


def write_mapping_file(filename, list_atoms_MN, dict_MN2PN, last_lines):
    """Write mapping file.

    A mapping file maps the mapping names to pdb names, for example:
    -----------------
    #Individual atoms
    M_G1_M        C32
    M_G1H1_M      H322
    M_G1H2_M      H321
    [...]
    -----------------

    Arguments
    ---------
    filename : str
        The output file name.
    list_atoms_MN : list
        This list contains the mapping names (str) of each atom.
    dict_MN2PN : dict
        This dict contains mapping names to pdb names (for all atoms).
    last_lines : str
        The last lines in the mapping file which contain some more info on the
        lipid ("#Water" and "#Whole molecule") which cannot be guessed
        automatically.

    Returns
    -------
    None
    """
    warning = False
    with open(filename, "w") as f:
        beginning = "#Individual atoms"
        f.write(f"{beginning}\n")
        for key in list_atoms_MN:
            if key in dict_MN2PN:
                f.write(f"{key:<13s} {dict_MN2PN[key]:<s}\n")
            else:
                warning = True
                print(f"!!! Atom {key} is not in dict which maps mapping "
                      f"names 2 pdb names !!!")
        f.write(f"{last_lines}\n")
    if warning:
        print("!!! Beware some Hs were missing !!!")


def write_def_file(filename, dict_MN2GN, dict_MN2PN, pdb_G):
    """Write def file.

    Arguments
    ---------
    filename : str
        The output file name.
    dict_MN2GN : dict
        This dict contains mapping names to generic names (for H only).
    dict_MN2PN : dict
        This dict contains mapping names to pdb names (for all atoms).
    pdb_G : networkx graph instance
        This the graph with the pdb names.

    Returns
    -------
    None
    """
    with open(filename, "w") as f:
        # Loop over all H.
        for key in dict_MN2GN:
            H_generic_name = dict_MN2GN[key]
            H_pdb_name = dict_MN2PN[key]
            C_neighbour_pdb_name = list(pdb_G.neighbors(H_pdb_name))[0]
            resname = "POPC"
            f.write(f"{H_generic_name} {resname} {C_neighbour_pdb_name} "
                    f"{H_pdb_name}\n")


if __name__ == "__main__":
    #############################
    # 0) Parse command line arguments
    #############################
    message = """This program generates a mapping file and a .def file from a pdb file containing one lipid only."""
    parser = argparse.ArgumentParser(description=message)
    parser.add_argument("-p", "--pdb", required=True, type=str,
                        help="pdb file containing a single lipid.")
    parser.add_argument("-l", "--lipid", required=True, type=str,
                        help="Name of lipid (e.g. POPC).")
    parser.add_argument("-om", "--omap", required=True, type=str,
                        help="Output mapping file name.")
    parser.add_argument("-od", "--odef", required=True, type=str,
                        help="Output .def file name.")
    parser.add_argument("--graph", action="store_true", help="Draw the graphs.")
    options = parser.parse_args()
    # Extract lists and dists from lipids_info.py.
    lipid_name = options.lipid # e.g. POPC
    list_atoms_MN = getattr(lipids_info, options.lipid+"_atoms_MN")
    list_bonds_MN = getattr(lipids_info, options.lipid+"_bonds_MN")
    dict_MN2GN = getattr(lipids_info, options.lipid+"_MN2GN")
    last_lines_in_mapping_file = getattr(lipids_info, options.lipid+
                                         "_last_lines_in_mapping_file")
    
    # In the following: MN = mapping name, PN = pdb name, GN = generic name.
    #############################
    # 1) Build graph G_MN with names from mapping file.
    # This graph is built from a dictionnary (lipids_info.py)
    #############################
    G_MN = nx.Graph()
    G_MN.add_nodes_from(list_atoms_MN)
    # Build connectivity.
    G_MN.add_edges_from(list_bonds_MN)
    # Plot graph properties.
    print(f"{lipid_name} graph from mapping names has {G_MN.number_of_nodes()} "
          f" nodes (atoms) and {G_MN.number_of_edges()} edges (bonds)")
    # These 2 vars will be used to check our 2nd graph is similar (same nb of nodes, connectivity, etc.).
    theoretical_nb_nodes = G_MN.number_of_nodes()
    theoretical_nb_edges = G_MN.number_of_edges()
    # Draw the graph G_MN and check manually all nodes and edges are correct.
    if options.graph:
        #nx.draw(G)
        #nx.draw_random(G)
        #nx.draw_circular(G)
        #nx.draw_spectral(G)
        nx.draw_kamada_kawai(G_MN, with_labels=True)
        #plt.savefig("graph_mapping_names.png")
        plt.show()

    #############################
    # 2) Build graph from a pdb file.
    # Now we read a pdb file, build a graph from the atom connectivity.
    #############################
    # Now read pdb in a pandas dataframe.
    df = pdb2pandasdf(options.pdb)
    # Build graph.
    print("Now infering the graph from the pdb structure.")
    G_PN = build_graph_from_pdb(df, theoretical_nb_nodes, theoretical_nb_edges)
    # Plot graph.
    if options.graph:
        nx.draw_kamada_kawai(G_PN, with_labels=True)
        plt.show()

    ##########################
    # 3) Match the 2 graphs  #
    ##########################
    graph_matching_obj = nx.isomorphism.GraphMatcher(G_MN, G_PN)
    if graph_matching_obj.is_isomorphic():
        print("Both graphs as isomorphic :-) (they match !).")
    else:
        exit("Graphs are not isomorphic (i.e. they do not match).")
    dict_MN2PN = graph_matching_obj.mapping
    #print("The mapping dictionnary is the following")
    #print(dict_MN2PN)

    ##########################
    # 4) Output mapping file #
    ##########################
    print("Writing mapping file")
    write_mapping_file(options.omap, list_atoms_MN, dict_MN2PN, \
                       last_lines_in_mapping_file)
    
    #######################
    # 5) Output .def file #
    #######################
    print("Writing .def file")
    write_def_file(options.odef, dict_MN2GN, dict_MN2PN, G_PN)
    
    print("Finished :-) !")
