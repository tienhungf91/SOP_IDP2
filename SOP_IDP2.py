#!/usr/bin/env python

from simtk.openmm import app
import simtk.openmm as omm
from simtk import unit
import time
import sys
import argparse
import build
import force

KELVIN_TO_KT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB / unit.kilocalorie_per_mole

parser = argparse.ArgumentParser(description='Coarse-grained SOP_IDP simulation using OpenMM')
parser.add_argument('-f','--sequence', type=str, help='input sequence')
parser.add_argument('-K','--monovalent_concentration', type=float, default='150.',
                    help='Monovalent concentration (mM) [150.0]')
parser.add_argument('-c','--cutoff', type=float, default='40.',
                    help='Cutoff distance for electrostatics (A) [40.0]')
parser.add_argument('-T','--temperature', type=float, default='20.',
                    help='Temperature (oC) [20.0]')
parser.add_argument('-t','--traj', type=str, default='md.dcd',
                    help='trajectory output')
parser.add_argument('-o','--output', type=str, default='md.out',
                    help='status and energy output')
parser.add_argument('-x','--frequency', type=int, default='10000',
                    help='output frequency')
parser.add_argument('-s','--step', type=int, default='10000',
                    help='Number of step [10000]')
args = parser.parse_args()

class simu:    ### structure to group all simulation parameter
    temp = 0.
    Kconc = 0.
    Nstep = 0
    epsilon = 0.
    cutoff = 40.

simu.temp = (args.temperature + 273.15)*unit.kelvin
simu.Nstep = args.step
simu.Kconc = args.monovalent_concentration
simu.cutoff = args.cutoff

T_unitless = simu.temp * KELVIN_TO_KT
print("T_unitless  ", T_unitless)
simu.epsilon = 296.0736276 - 619.2813716 * T_unitless + 531.2826741 * T_unitless**2 - 180.0369914 * T_unitless**3;
print("epsilon  ", simu.epsilon)
simu.l_Bjerrum = 332.0637*unit.angstroms / simu.epsilon
print("Bjerrum length  ", simu.l_Bjerrum / T_unitless)
simu.kappa = unit.sqrt (4*3.14159 * simu.l_Bjerrum * 2*simu.Kconc*6.022e-7 / (T_unitless * unit.angstrom**3))
print("kappa   ", simu.kappa)

forcefield = app.ForceField('SOP_IDP2.xml')
topology = None
positions = None

if args.sequence != None:
    print("Building from sequence %s ..." % args.sequence)
    topology, positions = build.build_by_seq(args.sequence, forcefield)
else:
    print("Need sequence !!!")
    sys.exit()

system = forcefield.createSystem(topology)

########## add force
force.add_bond_force (topology, system, 0)
force.add_exV_force  (topology, system, 1)
force.add_DH_force   (topology, system, simu, 2)
force.add_statistical_force (topology, system, 3)
totalforcegroup = 3

########## Simulation ############
integrator = omm.LangevinIntegrator(simu.temp, 0.5/unit.picosecond, 10*unit.femtoseconds)
#integrator = omm.GeodesicBAOABIntegrator(K_r=3, temperature=simu.temp, collision_rate=0.5/unit.picosecond, timestep=20.*unit.femtoseconds)            
properties = {'CudaPrecision': 'mixed'}

simulation = app.Simulation(topology, system, integrator)
#simulation = app.Simulation(topology, system, integrator, platform, properties)

simulation.context.setPositions(positions)
print("Initial energy   %f   kcal/mol" % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / unit.kilocalorie_per_mole))

state = simulation.context.getState(getPositions=True)
app.PDBFile.writeFile(topology, state.getPositions(), open("input.pdb", "w"), keepIds=True)

# print('Minimizing ...')
# simulation.minimizeEnergy(1.*unit.kilocalorie_per_mole, 5000)
simulation.context.setVelocitiesToTemperature(simu.temp)

simulation.reporters.append(app.DCDReporter(args.traj, args.frequency))
simulation.reporters.append(app.StateDataReporter(args.output, args.frequency, step=True, potentialEnergy=True, temperature=True, remainingTime=True, totalSteps=simu.Nstep, separator='  '))
#simulation.reporters.append(app.CheckpointReporter(args.res_file, int(args.frequency)*100))

print('Running simulation ...')
t0 = time.time()
simulation.step(simu.Nstep)
prodtime = time.time() - t0
print("Simulation speed: % .2e steps/day" % (86400*simu.Nstep/(prodtime)))
