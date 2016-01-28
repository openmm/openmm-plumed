/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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
#include "openmm/CustomExternalForce.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <iostream>
#include <string>
#include <vector>

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerPlumedReferenceKernelFactories();

void testForce() {
    // Create a System that applies a force based on the distance between two atoms.

    const int numParticles = 4;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions[i] = Vec3(i, 0.1*i, -0.3*i);
    }
    string script =
        "d: DISTANCE ATOMS=1,3\n"
        "BIASVALUE ARG=d";
    PlumedForce* plumed = new PlumedForce(script);
    system.addForce(plumed);
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Compute the forces and energy.

    State state = context.getState(State::Energy | State::Forces);
    Vec3 delta = positions[0]-positions[2];
    double dist = sqrt(delta.dot(delta));
    ASSERT_EQUAL_TOL(dist, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(-delta/dist, state.getForces()[0], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(), state.getForces()[1], 1e-5);
    ASSERT_EQUAL_VEC(delta/dist, state.getForces()[2], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(), state.getForces()[3], 1e-5);
}

void testMetadynamics() {
    // Create a System that does metadynamics within a one dimensional harmonic well.

    System system;
    system.addParticle(1.0);
    CustomExternalForce* external = new CustomExternalForce("x^2");
    external->addParticle(0);
    system.addForce(external);
    string script =
        "p: POSITION ATOM=1\n"
        "METAD ARG=p.x SIGMA=0.5 HEIGHT=0.1 PACE=1";
    PlumedForce* plumed = new PlumedForce(script);
    system.addForce(plumed);
    vector<Vec3> positions;
    positions.push_back(Vec3());
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Run a short simulation and check the energy at each step.

    vector<double> centers;
    for (int i = 0; i < 100; i++) {
        integ.step(1);
        State state = context.getState(State::Positions | State::Energy);
        double x = state.getPositions()[0][0];
        double expected = x*x;
        for (int j = 0; j < centers.size(); j++)
            expected += 0.1*exp(-(x-centers[j])*(x-centers[j])/(2*0.5*0.5));
        ASSERT_EQUAL_TOL(expected, state.getPotentialEnergy(), 1e-3);
        if (i > 0)
            centers.push_back(x);
    }
}

int main() {
    try {
        registerPlumedReferenceKernelFactories();
        testForce();
        testMetadynamics();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
