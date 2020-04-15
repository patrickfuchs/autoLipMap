"""This module file contains a number of definitions.

We are using 3 different atoms names:

- mapping names (e.g. "M_G1_M", "M_G1H1_M", etc.) -> unique for a given lipid
  --> will be called MN
- pdb names (e.g. C1, H11, H12, etc.) -> force field dependent
  --> will be called PN
- generic names (beta_1, beta_2, etc.) -> for hydrogens only (unique for a given lipid)
  --> will be called GN

For each lipid, 3 objects have to be defined, for example on POPC:

- POPC_atoms_MN: The list of atoms using mapping names (the order matters, it is used for ordering  the list of atoms in the .def and mapping file).
- POPC_bonds_MN: The list of bonds using mapping names.
- POPC_MN2GN: The dict which maps mapping names to generic names (for hydrogens only).

For using this script on a new lipid, e.g. DPPC, one has to create the 3 objects DPPC_atoms_MN, DPPC_bonds_MN and DPPC_MN2GN using the same rules.
"""

# This is the list of atoms using mapping names (MN).
# !!! Warning !!! The order of the atoms in this list is used to write the different files.
POPC_atoms_MN = ["M_G1_M", "M_G1H1_M", "M_G1H2_M", "M_G1O1_M",
                 "M_G1C2_M", "M_G1C2O1_M", "M_G1C3_M", "M_G1C3H1_M",
                 "M_G1C3H2_M", "M_G1C4_M", "M_G1C4H1_M", "M_G1C4H2_M",
                 "M_G1C5_M", "M_G1C5H1_M", "M_G1C5H2_M", "M_G1C6_M",
                 "M_G1C6H1_M", "M_G1C6H2_M", "M_G1C7_M", "M_G1C7H1_M",
                 "M_G1C7H2_M", "M_G1C8_M", "M_G1C8H1_M", "M_G1C8H2_M",
                 "M_G1C9_M", "M_G1C9H1_M", "M_G1C9H2_M", "M_G1C10_M",
                 "M_G1C10H1_M", "M_G1C10H2_M", "M_G1C11_M", "M_G1C11H1_M",
                 "M_G1C11H2_M", "M_G1C12_M", "M_G1C12H1_M", "M_G1C12H2_M",
                 "M_G1C13_M", "M_G1C13H1_M", "M_G1C13H2_M", "M_G1C14_M",
                 "M_G1C14H1_M", "M_G1C14H2_M", "M_G1C15_M", "M_G1C15H1_M",
                 "M_G1C15H2_M", "M_G1C16_M", "M_G1C16H1_M", "M_G1C16H2_M",
                 "M_G1C17_M", "M_G1C17H1_M", "M_G1C17H2_M", "M_G1C17H3_M",
                 "M_G2_M", "M_G2H1_M", "M_G2O1_M", "M_G2C2_M", "M_G2C2O1_M",
                 "M_G2C3_M", "M_G2C3H1_M", "M_G2C3H2_M", "M_G2C4_M",
                 "M_G2C4H1_M", "M_G2C4H2_M", "M_G2C5_M", "M_G2C5H1_M",
                 "M_G2C5H2_M", "M_G2C6_M", "M_G2C6H1_M", "M_G2C6H2_M",
                 "M_G2C7_M", "M_G2C7H1_M", "M_G2C7H2_M", "M_G2C8_M",
                 "M_G2C8H1_M", "M_G2C8H2_M", "M_G2C9_M", "M_G2C9H1_M",
                 "M_G2C9H2_M", "M_G2C10_M", "M_G2C10H1_M", "M_G2C11_M",
                 "M_G2C11H1_M", "M_G2C12_M", "M_G2C12H1_M", "M_G2C12H2_M",
                 "M_G2C13_M", "M_G2C13H1_M", "M_G2C13H2_M", "M_G2C14_M",
                 "M_G2C14H1_M", "M_G2C14H2_M", "M_G2C15_M", "M_G2C15H1_M",
                 "M_G2C15H2_M", "M_G2C16_M", "M_G2C16H1_M", "M_G2C16H2_M",
                 "M_G2C17_M", "M_G2C17H1_M", "M_G2C17H2_M", "M_G2C18_M",
                 "M_G2C18H1_M", "M_G2C18H2_M", "M_G2C19_M", "M_G2C19H1_M",
                 "M_G2C19H2_M", "M_G2C19H3_M", "M_G3_M", "M_G3H1_M",
                 "M_G3H2_M", "M_G3O1_M", "M_G3P2_M", "M_G3P2O1_M",
                 "M_G3P2O2_M", "M_G3O3_M", "M_G3C4_M", "M_G3C4H1_M",
                 "M_G3C4H2_M", "M_G3C5_M", "M_G3C5H1_M", "M_G3C5H2_M",
                 "M_G3N6_M", "M_G3N6C1_M", "M_G3N6C1H1_M", "M_G3N6C1H2_M",
                 "M_G3N6C1H3_M", "M_G3N6C2_M", "M_G3N6C2H1_M", "M_G3N6C2H2_M",
                 "M_G3N6C2H3_M", "M_G3N6C3_M", "M_G3N6C3H1_M", "M_G3N6C3H2_M",
                 "M_G3N6C3H3_M"]

# This is the list of bonds using mapping names (MN).
                 # glycerol
POPC_bonds_MN = [("M_G1_M", "M_G2_M"), ("M_G2_M", "M_G3_M"),
                 ("M_G1_M", "M_G1H1_M"), ("M_G1_M", "M_G1H2_M"), ("M_G1_M", "M_G1O1_M"),
                 ("M_G2_M", "M_G2H1_M"), ("M_G2_M", "M_G2O1_M"),
                 ("M_G3_M", "M_G3H1_M"), ("M_G3_M", "M_G3H2_M"), ("M_G3_M", "M_G3O1_M"),
                 # sn-1
                 ("M_G1O1_M", "M_G1C2_M"), ("M_G1C2_M", "M_G1C2O1_M"),
                 ("M_G1C2_M", "M_G1C3_M"), ("M_G1C3_M", "M_G1C3H1_M"), ("M_G1C3_M", "M_G1C3H2_M"),
                 ("M_G1C3_M", "M_G1C4_M"), ("M_G1C4_M", "M_G1C4H1_M"), ("M_G1C4_M", "M_G1C4H2_M"),
                 ("M_G1C4_M", "M_G1C5_M"), ("M_G1C5_M", "M_G1C5H1_M"), ("M_G1C5_M", "M_G1C5H2_M"),
                 ("M_G1C5_M", "M_G1C6_M"), ("M_G1C6_M", "M_G1C6H1_M"), ("M_G1C6_M", "M_G1C6H2_M"),
                 ("M_G1C6_M", "M_G1C7_M"), ("M_G1C7_M", "M_G1C7H1_M"), ("M_G1C7_M", "M_G1C7H2_M"),
                 ("M_G1C7_M", "M_G1C8_M"), ("M_G1C8_M", "M_G1C8H1_M"), ("M_G1C8_M", "M_G1C8H2_M"),
                 ("M_G1C8_M", "M_G1C9_M"), ("M_G1C9_M", "M_G1C9H1_M"), ("M_G1C9_M", "M_G1C9H2_M"),
                 ("M_G1C9_M", "M_G1C10_M"), ("M_G1C10_M", "M_G1C10H1_M"), ("M_G1C10_M", "M_G1C10H2_M"),
                 ("M_G1C10_M", "M_G1C11_M"), ("M_G1C11_M", "M_G1C11H1_M"), ("M_G1C11_M", "M_G1C11H2_M"),
                 ("M_G1C11_M", "M_G1C12_M"), ("M_G1C12_M", "M_G1C12H1_M"), ("M_G1C12_M", "M_G1C12H2_M"),
                 ("M_G1C12_M", "M_G1C13_M"), ("M_G1C13_M", "M_G1C13H1_M"), ("M_G1C13_M", "M_G1C13H2_M"),
                 ("M_G1C13_M", "M_G1C14_M"), ("M_G1C14_M", "M_G1C14H1_M"), ("M_G1C14_M", "M_G1C14H2_M"),
                 ("M_G1C14_M", "M_G1C15_M"), ("M_G1C15_M", "M_G1C15H1_M"), ("M_G1C15_M", "M_G1C15H2_M"),
                 ("M_G1C15_M", "M_G1C16_M"), ("M_G1C16_M", "M_G1C16H1_M"),("M_G1C16_M", "M_G1C16H2_M"),
                 ("M_G1C16_M", "M_G1C17_M"),("M_G1C17_M", "M_G1C17H1_M"), ("M_G1C17_M", "M_G1C17H2_M"), ("M_G1C17_M", "M_G1C17H3_M"),
                 # sn-2
                 ("M_G2O1_M", "M_G2C2_M"), ("M_G2C2_M", "M_G2C2O1_M"),
                 ("M_G2C2_M", "M_G2C3_M"), ("M_G2C3_M", "M_G2C3H1_M"), ("M_G2C3_M", "M_G2C3H2_M"),
                 ("M_G2C3_M", "M_G2C4_M"), ("M_G2C4_M", "M_G2C4H1_M"), ("M_G2C4_M", "M_G2C4H2_M"),
                 ("M_G2C4_M", "M_G2C5_M"), ("M_G2C5_M", "M_G2C5H1_M"), ("M_G2C5_M", "M_G2C5H2_M"),
                 ("M_G2C5_M", "M_G2C6_M"), ("M_G2C6_M", "M_G2C6H1_M"), ("M_G2C6_M", "M_G2C6H2_M"),
                 ("M_G2C6_M", "M_G2C7_M"), ("M_G2C7_M", "M_G2C7H1_M"), ("M_G2C7_M", "M_G2C7H2_M"),
                 ("M_G2C7_M", "M_G2C8_M"), ("M_G2C8_M", "M_G2C8H1_M"), ("M_G2C8_M", "M_G2C8H2_M"),
                 ("M_G2C8_M", "M_G2C9_M"), ("M_G2C9_M", "M_G2C9H1_M"), ("M_G2C9_M", "M_G2C9H2_M"),
                 ("M_G2C9_M", "M_G2C10_M"), ("M_G2C10_M", "M_G2C10H1_M"),
                 ("M_G2C10_M", "M_G2C11_M"), ("M_G2C11_M", "M_G2C11H1_M"),
                 ("M_G2C11_M", "M_G2C12_M"), ("M_G2C12_M", "M_G2C12H1_M"), ("M_G2C12_M", "M_G2C12H2_M"),
                 ("M_G2C12_M", "M_G2C13_M"), ("M_G2C13_M", "M_G2C13H1_M"), ("M_G2C13_M", "M_G2C13H2_M"),
                 ("M_G2C13_M", "M_G2C14_M"), ("M_G2C14_M", "M_G2C14H1_M"), ("M_G2C14_M", "M_G2C14H2_M"),
                 ("M_G2C14_M", "M_G2C15_M"), ("M_G2C15_M", "M_G2C15H1_M"), ("M_G2C15_M", "M_G2C15H2_M"),
                 ("M_G2C15_M", "M_G2C16_M"), ("M_G2C16_M", "M_G2C16H1_M"), ("M_G2C16_M", "M_G2C16H2_M"),
                 ("M_G2C16_M", "M_G2C17_M"), ("M_G2C17_M", "M_G2C17H1_M"), ("M_G2C17_M", "M_G2C17H2_M"),
                 ("M_G2C17_M", "M_G2C18_M"), ("M_G2C18_M", "M_G2C18H1_M"), ("M_G2C18_M", "M_G2C18H2_M"),
                 ("M_G2C18_M", "M_G2C19_M"), ("M_G2C19_M", "M_G2C19H1_M"), ("M_G2C19_M", "M_G2C19H2_M"), ("M_G2C19_M", "M_G2C19H3_M"),
                 # sn-3
                 ("M_G3O1_M", "M_G3P2_M"), ("M_G3P2_M", "M_G3P2O1_M"),
                 ("M_G3P2_M", "M_G3P2O2_M"), ("M_G3P2_M", "M_G3O3_M"),
                 ("M_G3O3_M", "M_G3C4_M"), ("M_G3C4_M", "M_G3C4H1_M"), ("M_G3C4_M", "M_G3C4H2_M"),                
                 ("M_G3C4_M", "M_G3C5_M"), ("M_G3C5_M", "M_G3C5H1_M"), ("M_G3C5_M", "M_G3C5H2_M"),
                 ("M_G3C5_M", "M_G3N6_M"),
                 ("M_G3N6_M", "M_G3N6C1_M"), ("M_G3N6C1_M", "M_G3N6C1H1_M"), ("M_G3N6C1_M", "M_G3N6C1H2_M"), ("M_G3N6C1_M", "M_G3N6C1H3_M"),
                 ("M_G3N6_M", "M_G3N6C2_M"), ("M_G3N6C2_M", "M_G3N6C2H1_M"), ("M_G3N6C2_M", "M_G3N6C2H2_M"), ("M_G3N6C2_M", "M_G3N6C2H3_M"),
                 ("M_G3N6_M", "M_G3N6C3_M"), ("M_G3N6C3_M", "M_G3N6C3H1_M"), ("M_G3N6C3_M", "M_G3N6C3H2_M"), ("M_G3N6C3_M", "M_G3N6C3H3_M")]

# This dictionnary maps from mapping names (MN) to generic names (GN). (hydrogens only)
POPC_MN2GN = {"M_G1H1_M": "g1_1",
              "M_G1H2_M": "g1_2",
              "M_G1C3H1_M": "palmytoyl_C2a",
              "M_G1C3H2_M": "palmytoyl_C2b",
              "M_G1C4H1_M": "palmytoyl_C3a",
              "M_G1C4H2_M": "palmytoyl_C3b",
              "M_G1C5H1_M": "palmytoyl_C4a",
              "M_G1C5H2_M": "palmytoyl_C4b",
              "M_G1C6H1_M": "palmytoyl_C5a",
              "M_G1C6H2_M": "palmytoyl_C5b",
              "M_G1C7H1_M": "palmytoyl_C6a",
              "M_G1C7H2_M": "palmytoyl_C6b",
              "M_G1C8H1_M": "palmytoyl_C7a",
              "M_G1C8H2_M": "palmytoyl_C7b",
              "M_G1C9H1_M": "palmytoyl_C8a",
              "M_G1C9H2_M": "palmytoyl_C8b",
              "M_G1C10H1_M": "palmytoyl_C9a",
              "M_G1C10H2_M": "palmytoyl_C9b",
              "M_G1C11H1_M": "palmytoyl_C10a",
              "M_G1C11H2_M": "palmytoyl_C10b",
              "M_G1C12H1_M": "palmytoyl_C11a",
              "M_G1C12H2_M": "palmytoyl_C11b",
              "M_G1C13H1_M": "palmytoyl_C12a",
              "M_G1C13H2_M": "palmytoyl_C12b",
              "M_G1C14H1_M": "palmytoyl_C13a",
              "M_G1C14H2_M": "palmytoyl_C13b",
              "M_G1C15H1_M": "palmytoyl_C14a",
              "M_G1C15H2_M": "palmytoyl_C14b",
              "M_G1C16H1_M": "palmytoyl_C15a",
              "M_G1C16H2_M": "palmytoyl_C15b",
              "M_G1C17H1_M": "palmytoyl_C16a",
              "M_G1C17H2_M": "palmytoyl_C16b",
              "M_G1C17H3_M": "palmytoyl_C16c",
              "M_G2H1_M": "g2_1",
              "M_G2C3H1_M": "oleoyl_C2a",
              "M_G2C3H2_M": "oleoyl_C2b",
              "M_G2C4H1_M": "oleoyl_C3a",
              "M_G2C4H2_M": "oleoyl_C3b",
              "M_G2C5H1_M": "oleoyl_C4a",
              "M_G2C5H2_M": "oleoyl_C4b",
              "M_G2C6H1_M": "oleoyl_C5a",
              "M_G2C6H2_M": "oleoyl_C5b",
              "M_G2C7H1_M": "oleoyl_C6a",
              "M_G2C7H2_M": "oleoyl_C6b",
              "M_G2C8H1_M": "oleoyl_C7a",
              "M_G2C8H2_M": "oleoyl_C7b",
              "M_G2C9H1_M": "oleoyl_C8a",
              "M_G2C9H2_M": "oleoyl_C8b",
              "M_G2C10H1_M": "oleoyl_C9a",
              "M_G2C11H1_M": "oleoyl_C10a",
              "M_G2C12H1_M": "oleoyl_C11a",
              "M_G2C12H2_M": "oleoyl_C11b",
              "M_G2C13H1_M": "oleoyl_C12a",
              "M_G2C13H2_M": "oleoyl_C12b",
              "M_G2C14H1_M": "oleoyl_C13a",
              "M_G2C14H2_M": "oleoyl_C13b",
              "M_G2C15H1_M": "oleoyl_C14a",
              "M_G2C15H2_M": "oleoyl_C14b",
              "M_G2C16H1_M": "oleoyl_C15a",
              "M_G2C16H2_M": "oleoyl_C15b",
              "M_G2C17H1_M": "oleoyl_C16a",
              "M_G2C17H2_M": "oleoyl_C16b",
              "M_G2C18H1_M": "oleoyl_C17a",
              "M_G2C18H2_M": "oleoyl_C17b",
              "M_G2C19H1_M": "oleoyl_C18a",
              "M_G2C19H2_M": "oleoyl_C18b",
              "M_G2C19H3_M": "oleoyl_C18c",
              "M_G3H1_M": "g3_1",
              "M_G3H2_M": "g3_1",
              "M_G3C4H1_M": "alpha1",
              "M_G3C4H2_M": "alpha2",
              "M_G3C5H1_M": "beta1",
              "M_G3C5H2_M": "beta2",
              "M_G3N6C1H1_M": "gamma1_1",
              "M_G3N6C1H2_M": "gamma1_2",
              "M_G3N6C1H3_M": "gamma1_3",
              "M_G3N6C2H1_M": "gamma2_1",
              "M_G3N6C2H2_M": "gamma2_2",
              "M_G3N6C2H3_M": "gamma2_3",
              "M_G3N6C3H1_M": "gamma3_1",
              "M_G3N6C3H2_M": "gamma3_2",
              "M_G3N6C3H3_M": "gamma3_3"}

POPC_last_lines_in_mapping_file = """#Water
M_OW_M           OW
M_HW1_M          HW1
M_HW2_M          HW2
#Whole molecules
M_POPC_M         POPC
M_NA_M           NA
M_CHOL_M         CHO
M_SOL_M          SOL"""
