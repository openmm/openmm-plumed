/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This tests the Reference implementation of PlumedForce.
 */

#include "PlumedForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerPlumedOpenCLKernelFactories();

void testForce() {
    // Create a chain of particles connected by bonds.
    
    const int numBonds = 10;
    const int numParticles = numBonds+1;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i] = Vec3(i, 0.1*i, -0.3*i);
    }
    PlumedForce* force = new PlumedForce();
    system.addForce(force);
    for (int i = 0; i < numBonds; i++)
        force->addBond(i, i+1, 1.0+sin(0.8*i), cos(0.3*i));
    
    // Compute the forces and energy.

    VerletIntegrator integ(1.0);
    Platform& platform = Platform::getPlatformByName("OpenCL");
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    
    // See if the energy is correct.
    
    double expectedEnergy = 0;
    for (int i = 0; i < numBonds; i++) {
        double length = 1.0+sin(0.8*i);
        double k = cos(0.3*i);
        Vec3 delta = positions[i+1]-positions[i];
        double dr = sqrt(delta.dot(delta))-length;
        expectedEnergy += k*dr*dr*dr*dr;
    }
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);

    // Validate the forces by moving each particle along each axis, and see if the energy changes by the correct amount.
    
    double offset = 1e-3;
    for (int i = 0; i < numParticles; i++)
        for (int j = 0; j < 3; j++) {
            vector<Vec3> offsetPos = positions;
            offsetPos[i][j] = positions[i][j]-offset;
            context.setPositions(offsetPos);
            double e1 = context.getState(State::Energy).getPotentialEnergy();
            offsetPos[i][j] = positions[i][j]+offset;
            context.setPositions(offsetPos);
            double e2 = context.getState(State::Energy).getPotentialEnergy();
            ASSERT_EQUAL_TOL(state.getForces()[i][j], (e1-e2)/(2*offset), 1e-3);
        }
}

void testChangingParameters() {
    const double k = 1.5;
    const double length = 0.5;
    Platform& platform = Platform::getPlatformByName("OpenCL");
    
    // Create a system with one bond.
    
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    PlumedForce* force = new PlumedForce();
    force->addBond(0, 1, length, k);
    system.addForce(force);
    vector<Vec3> positions(2);
    positions[0] = Vec3(1, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    
    // Check the energy.
    
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(k*pow(1.0-length, 4), state.getPotentialEnergy(), 1e-5);
    
    // Modify the parameters.
    
    const double k2 = 2.2;
    const double length2 = 0.9;
    force->setBondParameters(0, 0, 1, length2, k2);
    force->updateParametersInContext(context);
    state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(k2*pow(1.0-length2, 4), state.getPotentialEnergy(), 1e-5);
}

int main(int argc, char* argv[]) {
    try {
        registerPlumedOpenCLKernelFactories();
        if (argc > 1)
            Platform::getPlatformByName("OpenCL").setPropertyDefaultValue("OpenCLPrecision", string(argv[1]));
        testForce();
        testChangingParameters();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
