import simtk.openmm as mm
import simtk.unit as unit
from openmmplumed import PlumedForce
import numpy as np
import unittest


class TestTorchForce(unittest.TestCase):

    def testForce(self):
        # Create a System that applies a force based on the distance between two atoms.

        numParticles = 4
        system = mm.System()
        positions = np.empty((numParticles, 3))
        for i in range(numParticles):
            system.addParticle(1.0)
            positions[i] = [i, 0.1*i, -0.3*i]
        script = '''
            d: DISTANCE ATOMS=1,3
            BIASVALUE ARG=d
        '''
        force = PlumedForce(script)
        system.addForce(force)
        integ = mm.LangevinIntegrator(300.0, 1.0, 1.0)
        context = mm.Context(system, integ, mm.Platform.getPlatformByName('Reference'))
        context.setPositions(positions)

        # Compute the forces and energy.

        state = context.getState(getEnergy=True, getForces=True)
        delta = positions[0] - positions[2] 
        dist = np.sqrt(np.sum(delta**2))
        zero = np.zeros(3)
        self.assertAlmostEqual(dist, state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))
        self.assertTrue(np.allclose(-delta/dist, state.getForces(asNumpy=True)[0]))
        self.assertTrue(np.allclose(zero, state.getForces(asNumpy=True)[1]))
        self.assertTrue(np.allclose(delta/dist, state.getForces(asNumpy=True)[2]))
        self.assertTrue(np.allclose(zero, state.getForces(asNumpy=True)[3]))


if __name__ == '__main__':
    unittest.main()