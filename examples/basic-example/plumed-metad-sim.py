from openmm.app import *
from openmm import *                    
from openmm.unit import *
import bz2
from openmmplumed import PlumedForce

sim_temp = 300.0 * kelvin
H_mass = 4.0 * amu
time_step = 0.002 * picosecond  
nb_cutoff = 10.0 * angstrom                                                                                                          
box_padding = 12.0 * angstrom
salt_conc = 0.15 * molar
current_file="apo-met-meta"

# Misc parameters                                                 
restraint_distance = 0.0 * angstroms                              
restart_freq = 10
log_freq = 10                                                                                                                         
prd_steps = 100  
iterations = 500000

# unpack the equilbirated system
bz2_file = f'system.xml.bz2'
with bz2.open(bz2_file, 'rb') as infile:
    system = XmlSerializer.deserialize(infile.read().decode())

# unpack the state
bz2_file = f'state.xml.bz2'
with bz2.open(bz2_file, 'rb') as infile:
    state = XmlSerializer.deserialize(infile.read().decode())

# get topology info from equilibrated protein
pdb = PDBFile('receptor_equilibrated.pdb')
positions=pdb.positions
topology=pdb.topology

# Add barostat
system.addForce(openmm.MonteCarloBarostat(1*bar, sim_temp))

# Generate collective variables
# here we have 2 distances together using the PLUMED syntax
# Remember: Add +1 to each atom index when adding to PLUMED
script = """
d1: DISTANCE ATOMS=2465,975
d2: DISTANCE ATOMS=2465,642
METAD ARG=d1,d2 SIGMA=0.1,0.1 HEIGHT=0.3 PACE=50 FILE=colvar"""
system.addForce(PlumedForce(script))

integrator = LangevinIntegrator(sim_temp, 1.0/picosecond, time_step)

# put everything together into a simulation object
simulation = Simulation(pdb.topology, system, integrator)

# Only use the below line if you're loading from a previous checkpoint
# simulation.loadCheckpoint("%s.chk" % current_file)

# Reporters
simulation.reporters.append(DCDReporter(''+current_file+''+ '.dcd', restart_freq))
simulation.reporters.append(CheckpointReporter(''+current_file+''+ '.chk', 
                                                min(prd_steps, 10*restart_freq)))
simulation.reporters.append(StateDataReporter(open('log.' + current_file+'', 'w'),
                            log_freq, 
                            step=True,
                            potentialEnergy=True,
                            kineticEnergy=True,
                            totalEnergy=True,
                            temperature=True,
                            volume=True,
                            density=True,
                            speed=True))

# Set positions
simulation.context.setPositions(positions)

# minimize
print('  initial : %s' % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))
simulation.minimizeEnergy()
print('  final : %s' % (simulation.context.getState(getEnergy=True).getPotentialEnergy()))

for i in range(iterations):
    simulation.step(10)


