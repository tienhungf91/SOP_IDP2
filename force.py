#!/usr/bin/env python

import simtk.openmm as omm
from simtk import unit
import itertools as it

def get_sigma (x):
    if x == "B":
        return 1.90
    elif x == "A":
        return 2.52
    elif x == "R":
        return 3.28
    elif x == "N":
        return 2.84
    elif x == "D":
        return 2.79
    elif x == "C":
        return 2.74
    elif x == "E":
        return 2.96
    elif x == "Q":
        return 3.01
    elif x == "H":
        return 3.04
    elif x == "I" or x == "L" or x == "M":
        return 3.09
    elif x == "K" or x == "F":
        return 3.18
    elif x == "P":
        return 2.78
    elif x == "S":
        return 2.59
    elif x == "T":
        return 2.81
    elif x == "W":
        return 3.39
    elif x == "Y":
        return 3.23
    elif x == "V":
        return 2.93

############################################
########## bond force
############################################
def add_bond_force(topology, system, forcegroup):
    energy_function =  '- k_b * R0_2 * log(1 - (r-r_ref)^2 / R0_2)'
    bondforce = omm.CustomBondForce(energy_function)
    bondforce.addGlobalParameter('k_b', 10.0*unit.kilocalorie_per_mole/(unit.angstroms**2))
    bondforce.addGlobalParameter('R0_2', 4.0*unit.angstroms**2)
    bondforce.addPerBondParameter('r_ref')

    for bond in topology.bonds():
        bond_length = get_sigma (bond[0].name) + get_sigma (bond[1].name)
#        print("Add bond for %s  %s,  %d  %d" % (bond[0], bond[1], bond[0].index, bond[1].index))
#        print("Bond length %f" % bond_length)
        bondforce.addBond (bond[0].index, bond[1].index, [bond_length*unit.angstroms])

    bondforce.setForceGroup(forcegroup)
    system.addForce(bondforce)

##############################################
def take_3_at_a_time(iterable):
    i, nxt1, nxt2 = it.tee(iterable, 3)
    j = it.chain(it.islice(nxt1, 1, None), [None])
    k = it.chain(it.islice(nxt2, 2, None), [None, None])
    return zip(i, j, k)

###############################################
###### excluded volume force
###############################################
def add_exV_force(topology, system, forcegroup):
    energy_function = 'ep_loc * (sigma / r)^6'

    exVforce = omm.CustomBondForce(energy_function)
    exVforce.addGlobalParameter('ep_loc', 1.*unit.kilocalorie_per_mole)
    exVforce.addPerBondParameter('sigma')

    def add_f (atm1, atm2, exVforce):
        bond_length = get_sigma(atm1.name) + get_sigma(atm2.name)
#        print("Add exV for %s   %s,  %d  %d" % (atm1, atm2, atm1.index, atm2.index))
#        print("Bond length %f" % bond_length)
        exVforce.addBond (atm1.index, atm2.index, [bond_length*unit.angstroms])

    for res1, res2, res3 in take_3_at_a_time(topology.residues()):
        for atm1 in res1.atoms():
            if atm1.name == "B":
                if res2 != None:
                    for atm2 in res2.atoms():
                        if atm2.name != "B":
                            add_f (atm1, atm2, exVforce)
                if res3 != None:
                    for atm2 in res3.atoms():
                        if atm2.name == "B":
                            add_f (atm1, atm2, exVforce)
            else:
                if res2 != None:
                    for atm2 in res2.atoms():
                        if atm2.name == "B":
                            add_f (atm1, atm2, exVforce)

    exVforce.setForceGroup(forcegroup)
    system.addForce(exVforce)

#######################################
####### Debye-Huckel
#######################################
def add_DH_force(topology, system, simu, forcegroup):
    DHforce = omm.CustomNonbondedForce("U_ee * q1*q2 * exp(-kappa*r)/r")
    DHforce.addGlobalParameter("U_ee", simu.l_Bjerrum * unit.kilocalorie_per_mole / unit.elementary_charge**2)
    DHforce.addGlobalParameter("kappa", simu.kappa)
    DHforce.addPerParticleParameter('q')

    for atom in topology.atoms():
        if atom.name == "R" or atom.name == "K":
            DHforce.addParticle([1*unit.elementary_charge])
        elif atom.name == "D" or atom.name == "E":
            DHforce.addParticle([-1*unit.elementary_charge])
        else:
            DHforce.addParticle([0*unit.elementary_charge])

#    2    4    6    8
#    |    |    |    |
#    1 -- 3 -- 5 -- 7 ...
#    where 1, 3, 5, and 7 are backbone beads, and 2, 4, 6, and 8 are side-chain beads.
#    Then, the following interactions are excluded (meaning, no electrostatics, and only repulsive vdw):
#    1-2 ; 1-3; 1-4; 1-5
#    2-3
#    3-4; 3-5; 3-6 ; 3-7
#    4-5
#    ...
#    This leaves 2-4, 2-5, 2-6 as attractive interactions

    for res1, res2, res3 in take_3_at_a_time(topology.residues()):
        for atm in res1.atoms():
            if atm.name == "B":
                if "GLY" not in res1.name:
                    DHforce.addExclusion (atm.index, atm.index+1)
#                    print("Excluding   %d  %d" % (atm.index, atm.index+1))
                if res2 != None:
                    for atm2 in res2.atoms():
                        DHforce.addExclusion (atm.index, atm2.index)
#                        print("Excluding   %d  %d" % (atm.index, atm2.index))
                if res3 != None:
                    for atm2 in res3.atoms():
                        if atm2.name == "B":
                            DHforce.addExclusion (atm.index, atm2.index)
#                            print("Excluding   %d  %d" % (atm.index, atm2.index))
            elif res2 != None:
                DHforce.addExclusion (atm.index, atm.index+1)
#                print("Excluding   %d  %d" % (atm.index, atm.index+1))

    DHforce.setCutoffDistance(simu.cutoff)
    DHforce.setForceGroup(forcegroup)
    DHforce.setNonbondedMethod(omm.CustomNonbondedForce.CutoffNonPeriodic)
    system.addForce(DHforce)

#################################################################
####   Knowledge-based potential for inter-residue interactions
################################################################
def add_statistical_force(topology, system, forcegroup):
    #### Betancourt-Thirumalai potential
    ## Betancourt, Thirumalai. Prot Sci 1999, 8, 361
    epsilon = [
      # G     A      R      N      D      C      E      Q      H      I      L      K      M      F      P      S      T      W      Y      V     B   (internal unit - kJ/mol)
      #GLY   ALA    ARG    ASN    ASP    CYS    GLU    GLN    HIS    ILE    LEU    LYS    MET    PHE    PRO    SER    THR    TRP    TYR    VAL    BBone
     -0.50, -0.07,  0.35,  0.25,  0.42, -0.22,  1.19,  0.50,  0.57,  0.52,  0.35,  0.30,  0.20,  0.27, -0.03,  0.25,  0.00, -0.59, -0.10,  0.10,  0.74,  # G  GLY
     -0.07, -0.50,  0.67,  0.59,  0.74, -0.64,  1.07,  0.52,  0.52, -0.87, -0.92,  0.50, -0.57, -0.82,  0.17,  0.37,  0.00, -0.99, -0.37, -0.94,  0.74,  # A  ALA
      0.35,  0.67,  0.32,  0.05, -1.76,  0.79, -1.86, -0.30,  0.10,  0.45,  0.22,  1.24,  0.42,  0.20, -0.05,  0.30,  0.00, -1.02, -0.92,  0.42,  0.74,  # R  ARG
      0.25,  0.59,  0.05, -0.10, -0.30,  0.69, -0.02, -0.12,  0.25,  1.36,  0.89, -0.35,  0.79,  0.72,  0.32,  0.35,  0.00, -0.22,  0.02,  0.97,  0.74,  # N  ASN
      0.42,  0.74, -1.76, -0.30,  0.67,  0.94,  0.99,  0.30, -0.55,  1.34,  1.54, -1.71,  1.54,  1.19,  0.62,  0.02,  0.00,  0.15, -0.17,  1.64,  0.74,  # D  ASP
     -0.22, -0.64,  0.79,  0.69,  0.94, -3.32,  1.14,  0.10, -0.47, -1.19, -1.24,  0.87, -1.21, -1.31, -0.45,  0.22,  0.00, -1.83, -0.40, -1.26,  0.74,  # C  CYS
      1.19,  1.07, -1.86, -0.02,  0.99,  1.14,  1.12,  0.25, -0.27,  0.94,  0.92, -2.16,  0.59,  0.84,  0.64,  0.25,  0.00, -0.37, -0.40,  1.02,  0.74,  # E  GLU
      0.50,  0.52, -0.30, -0.12,  0.30,  0.10,  0.25,  0.35,  0.55,  0.35,  0.20, -0.50, -0.02, -0.10, -0.12,  0.62,  0.00, -0.27, -0.45,  0.42,  0.74,  # Q  GLN
      0.57,  0.52,  0.10,  0.25, -0.55, -0.47, -0.27,  0.55, -0.82,  0.47,  0.25,  0.64, -0.42, -0.47, -0.12,  0.37,  0.00, -1.14, -0.52,  0.45,  0.74,  # H  HIS
      0.52, -0.87,  0.45,  1.36,  1.34, -1.19,  0.94,  0.35,  0.47, -1.49, -1.96,  0.52, -1.49, -1.61,  0.12,  0.87,  0.00, -1.61, -0.82, -1.69,  0.74,  # I  ILE
      0.35, -0.92,  0.22,  0.89,  1.54, -1.24,  0.92,  0.20,  0.25, -1.96, -2.01,  0.40, -1.69, -1.93, -0.20,  0.64,  0.00, -1.74, -1.09, -1.98,  0.74,  # L  LEU
      0.30,  0.50,  1.24, -0.35, -1.71,  0.87, -2.16, -0.50,  0.64,  0.52,  0.40,  0.94,  0.55,  0.27,  0.30,  0.25,  0.00, -0.69, -0.99,  0.40,  0.74,  # K  LYS
      0.20, -0.57,  0.42,  0.79,  1.54, -1.21,  0.59, -0.02, -0.42, -1.49, -1.69,  0.55, -1.39, -2.21, -0.40,  0.79,  0.00, -2.33, -1.26, -1.17,  0.74,  # M  MET
      0.27, -0.82,  0.20,  0.72,  1.19, -1.31,  0.84, -0.10, -0.47, -1.61, -1.93,  0.27, -2.21, -2.03, -0.47,  0.25,  0.00, -1.93, -1.21, -1.66,  0.74,  # F  PHE
     -0.03,  0.17, -0.05,  0.32,  0.62, -0.45,  0.64, -0.12, -0.12,  0.12, -0.20,  0.30, -0.40, -0.47, -0.17,  0.42,  0.00, -1.81, -0.99, -0.20,  0.74,  # P  PRO
      0.25,  0.37,  0.30,  0.35,  0.02,  0.22,  0.25,  0.62,  0.37,  0.87,  0.64,  0.25,  0.79,  0.25,  0.42,  0.32,  0.00,  0.17,  0.17,  0.62,  0.74,  # S  SER
      0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.74,  # T  THR
     -0.59, -0.99, -1.02, -0.22,  0.15, -1.83, -0.37, -0.27, -1.14, -1.61, -1.74, -0.69, -2.33, -1.93, -1.81,  0.17,  0.00, -1.83, -1.36, -1.54,  0.74,  # W  TRP
     -0.10, -0.37, -0.92,  0.02, -0.17, -0.40, -0.40, -0.45, -0.52, -0.82, -1.09, -0.99, -1.26, -1.21, -0.99,  0.17,  0.00, -1.36, -0.67, -0.67,  0.74,  # Y  TYR
      0.10, -0.94,  0.42,  0.97,  1.64, -1.26,  1.02,  0.42,  0.45, -1.69, -1.98,  0.40, -1.17, -1.66, -0.20,  0.62,  0.00, -1.54, -0.67, -1.78,  0.74,  # V  VAL
      0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74,  0.74]  # B  BBone

              #  BB    SC
#    scaling = [0.5021, 1.004,  # BB   Convert BB-BB and BB-SC to kJ/mol by scaling a factor of 4.184
#               1.004,  0.3035] # SC   Convert SC-SC to kJ/mol by scaling a factor of kT=0.593 (4.184 was already in the matrix)

    # after round 1
#    scaling = [0.4554, 0.4692,
#               0.4692, 0.3035]

    # after round 2
    scaling = [0.4352, 0.5731,
               0.5731, 0.3035]

    energy_function = 'scale * abs(eps-1.74) * R6*(R6 - 2);'
    energy_function += 'R6=((sig1 + sig2) / r)^6; eps=epsilon(type1, type2); scale=scaling(sc_type1, sc_type2);'

    interresidueforce = omm.CustomNonbondedForce(energy_function)
    interresidueforce.addPerParticleParameter('sig')
    interresidueforce.addTabulatedFunction('epsilon', omm.Discrete2DFunction(21, 21, epsilon))
    interresidueforce.addTabulatedFunction('scaling', omm.Discrete2DFunction(2, 2, scaling))
    interresidueforce.addPerParticleParameter('type')
    interresidueforce.addPerParticleParameter('sc_type')

    for atom in topology.atoms():
        sigma = get_sigma(atom.name)
        if atom.name == "B":
            interresidueforce.addParticle([sigma*unit.angstroms, 20, 0])
        elif atom.name == "A":
            interresidueforce.addParticle([sigma*unit.angstroms, 1, 1])
        elif atom.name == "R":
            interresidueforce.addParticle([sigma*unit.angstroms, 2, 1])
        elif atom.name == "N":
            interresidueforce.addParticle([sigma*unit.angstroms, 3, 1])
        elif atom.name == "D":
            interresidueforce.addParticle([sigma*unit.angstroms, 4, 1])
        elif atom.name == "C":
            interresidueforce.addParticle([sigma*unit.angstroms, 5, 1])
        elif atom.name == "E":
            interresidueforce.addParticle([sigma*unit.angstroms, 6, 1])
        elif atom.name == "Q":
            interresidueforce.addParticle([sigma*unit.angstroms, 7, 1])
        elif atom.name == "H":
            interresidueforce.addParticle([sigma*unit.angstroms, 8, 1])
        elif atom.name == "I":
            interresidueforce.addParticle([sigma*unit.angstroms, 9, 1])
        elif atom.name == "L":
            interresidueforce.addParticle([sigma*unit.angstroms, 10, 1])
        elif atom.name == "K":
            interresidueforce.addParticle([sigma*unit.angstroms, 11, 1])
        elif atom.name == "M":
            interresidueforce.addParticle([sigma*unit.angstroms, 12, 1])
        elif atom.name == "F":
            interresidueforce.addParticle([sigma*unit.angstroms, 13, 1])
        elif atom.name == "P":
            interresidueforce.addParticle([sigma*unit.angstroms, 14, 1])
        elif atom.name == "S":
            interresidueforce.addParticle([sigma*unit.angstroms, 15, 1])
        elif atom.name == "T":
            interresidueforce.addParticle([sigma*unit.angstroms, 16, 1])
        elif atom.name == "W":
            interresidueforce.addParticle([sigma*unit.angstroms, 17, 1])
        elif atom.name == "Y":
            interresidueforce.addParticle([sigma*unit.angstroms, 18, 1])
        elif atom.name == "V":
            interresidueforce.addParticle([sigma*unit.angstroms, 19, 1])

    for res1, res2, res3 in take_3_at_a_time(topology.residues()):
        for atm in res1.atoms():
            if atm.name == "B":
                if "GLY" not in res1.name:
                    interresidueforce.addExclusion (atm.index, atm.index+1)
                if res2 != None:
                    for atm2 in res2.atoms():
                        interresidueforce.addExclusion (atm.index, atm2.index)
                if res3 != None:
                    for atm2 in res3.atoms():
                        if atm2.name == "B":
                            interresidueforce.addExclusion (atm.index, atm2.index)
            elif res2 != None:
                interresidueforce.addExclusion (atm.index, atm.index+1)

    interresidueforce.setCutoffDistance(30.*unit.angstroms)
    interresidueforce.setForceGroup(forcegroup)
    interresidueforce.setNonbondedMethod(omm.CustomNonbondedForce.CutoffNonPeriodic)
    system.addForce(interresidueforce)
