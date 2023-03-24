#!/usr/bin/env python

from simtk.openmm import app
from simtk import unit
import simtk.openmm as omm

#def get_atom_index(res):
#    i = 0
#    j = 0
#    for atom in res.atoms():
#        if atom.name == 'B':
#            i = atom.index
#        else:
#            j = atom.index
#    return [i, j]
#
#def get_res_from_index(idx, topology):
#    for res in topology.residues():
#        if res.id == idx:
#            return res

#########################################
def build_by_seq(seq, forcefield):
    step = 0.22
    geometry = {'ALA': {'B': [0, 0, 0], 'A': [0.313, -0.313, 0.]},
                'ARG': {'B': [0, 0, 0], 'R': [0.366, -0.366, 0.]},
                'ASN': {'B': [0, 0, 0], 'N': [0.335, -0.335, 0.]},
                'ASP': {'B': [0, 0, 0], 'D': [0.332, -0.332, 0.]},
                'CYS': {'B': [0, 0, 0], 'C': [0.328, -0.328, 0.]},
                'GLU': {'B': [0, 0, 0], 'E': [0.344, -0.344, 0.]},
                'GLN': {'B': [0, 0, 0], 'Q': [0.347, -0.347, 0.]},
                'GLY': {'B': [0, 0, 0]},# 'G': [0.077, -0.077, 0.]},
                'HIS': {'B': [0, 0, 0], 'H': [0.349, -0.349, 0.]},
                'ILE': {'B': [0, 0, 0], 'I': [0.353, -0.353, 0.]},
                'LEU': {'B': [0, 0, 0], 'L': [0.353, -0.353, 0.]},
                'LYS': {'B': [0, 0, 0], 'K': [0.359, -0.359, 0.]},
                'MET': {'B': [0, 0, 0], 'M': [0.353, -0.353, 0.]},
                'PHE': {'B': [0, 0, 0], 'F': [0.359, -0.359, 0.]},
                'PRO': {'B': [0, 0, 0], 'P': [0.331, -0.331, 0.]},
                'SER': {'B': [0, 0, 0], 'S': [0.317, -0.317, 0.]},
                'THR': {'B': [0, 0, 0], 'T': [0.333, -0.333, 0.]},
                'TRP': {'B': [0, 0, 0], 'W': [0.374, -0.374, 0.]},
                'TYR': {'B': [0, 0, 0], 'Y': [0.363, -0.363, 0.]},
                'VAL': {'B': [0, 0, 0], 'V': [0.342, -0.342, 0.]}}
    name_map = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    topo = omm.app.Topology()
    chain = topo.addChain('X')
    atoms = []
    positions = []
    idx_offset = 0
    transformStack = [[0, 0, 0]]
    stackOffsets = []
    last_idx = None

    for i, resSymbol in enumerate(seq):
        symbol = name_map[resSymbol]
        geometry_for_res = geometry[symbol]
        if i == 0 or i == len(seq) - 1:
            symbol = symbol + "T"

        res = topo.addResidue(symbol, chain)
        for atom in forcefield._templates[symbol].atoms:
            atoms.append(topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res))
            if atom.name in geometry_for_res:
                if i % 2 == 0:
                    positions.append(geometry_for_res[atom.name])
                else:
                    new_geo = [-element for element in geometry_for_res[atom.name]]
                    positions.append(new_geo)
            else:
                print ("Residue % not found!!" % atom.name)
                exit()

        for bond in forcefield._templates[symbol].bonds:
#            print "Adding bond  %s - %s" % (atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset])
            topo.addBond(atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset])

        curr_idx = None
        for bond in forcefield._templates[symbol].externalBonds:
            a = atoms[bond + idx_offset]
            curr_idx = bond + idx_offset

        if curr_idx != None:
            loc_to = positions[curr_idx]
        else:
            loc_to = None

        if last_idx != None:
            loc_from = positions[last_idx]
        else:
            loc_from = [0, 0, 0]

        if loc_to != None and loc_from != None:
            curr_t = [loc_from[i] - loc_to[i] for i in range(0, 3)]
            curr_t[0] += step
            curr_t[1] += step
            curr_t[2] += step
            transformStack.append(curr_t)

        if last_idx != None and curr_idx != None:
#            print "Adding external bond  to %s  %s" % (atoms[last_idx], atoms[curr_idx])
            topo.addBond(atoms[last_idx], atoms[curr_idx])
        if curr_idx != None:
            last_idx = curr_idx

        idx_offset += len(forcefield._templates[symbol].atoms)
        stackOffsets.append(len(forcefield._templates[symbol].atoms))

    def transform_geometry(atoms, transform):
        return [atoms[i] + transform[i] for i in range(0, 3)]

    transformed_positions = []
    idx = 0
    currStackIdx = 0
    transformStack.append([0, 0, 0])
    curr_t = transformStack[0]
    for offset in stackOffsets:
        for i in range(0, offset):
            transformed_positions.append(transform_geometry(positions[idx], curr_t))
            idx += 1
        currStackIdx += 1
        curr_t = [curr_t[i] + transformStack[currStackIdx][i] for i in range(0, 3)]

    return topo, transformed_positions
