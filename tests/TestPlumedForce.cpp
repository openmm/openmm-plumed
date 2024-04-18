/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016-2023 Stanford University and the Authors.      *
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
#include "openmm/NonbondedForce.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

FILE* logstream;

void testForce(Platform& platform) {
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
    plumed->setLogStream(logstream);
    system.addForce(plumed);
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);

    // Compute the forces and energy.

    State state = context.getState(State::Energy | State::Forces);
    Vec3 delta = positions[0]-positions[2];
    double dist = sqrt(delta.dot(delta));
    Vec3 zero;
    ASSERT_EQUAL_TOL(dist, state.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(-delta/dist, state.getForces()[0], 1e-5);
    ASSERT_EQUAL_VEC(zero, state.getForces()[1], 1e-5);
    ASSERT_EQUAL_VEC(delta/dist, state.getForces()[2], 1e-5);
    ASSERT_EQUAL_VEC(zero, state.getForces()[3], 1e-5);
}

void testMetadynamics(Platform& platform) {
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
    plumed->setLogStream(logstream);
    system.addForce(plumed);
    vector<Vec3> positions;
    positions.push_back(Vec3());
    LangevinIntegrator integ(300.0, 1.0, 1.0);
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
        ASSERT(fabs(expected - state.getPotentialEnergy()) < 0.2);
        centers.push_back(x);
    }
}

void testWellTemperedMetadynamics(Platform& platform) {

    // Simulation parameters
    const double height0 = 0.1;
    const double sigma = 0.5;
    const double temperatue = 300.0;
    const double delta_temperature = 30.0;
    // Note: BIASFACTOR = temperature + delta_temperature / temperature,
    //       so PLUMED has to know a temperature to compupute
    //       delta_temperature from BIASFACTOR
    //       (https://www.plumed.org/doc-master/user-doc/html/belfast-6.html).
    const string script =
        "p: POSITION ATOM=1\n"
        "METAD ARG=p.x HEIGHT=0.1 SIGMA=0.5 BIASFACTOR=1.1 PACE=1 TEMP=300";

    // Create a system within a one dimensional harmonic potential
    System system;
    system.addParticle(1.0);
    CustomExternalForce* external = new CustomExternalForce("x^2");
    external->addParticle(0);
    system.addForce(external);

    // Create a well-tempered metadynamics simulation
    PlumedForce* plumed = new PlumedForce(script);
    plumed->setLogStream(logstream);
    plumed->setTemperature(temperatue); // This is tested here!
    system.addForce(plumed);
    LangevinIntegrator integ(temperatue, 1.0, 1.0);
    Context context(system, integ, platform);
    context.setPositions({Vec3()});

    // Run the simulation and compare potential energy
    vector<double> centers, heights;
    for (int i = 0; i < 100; i++) {
        integ.step(1);
        State state = context.getState(State::Positions | State::Energy);
        double x = state.getPositions()[0][0];

        // Compute bias
        double bias = 0;
        for (int j = 0; j < centers.size(); j++)
            bias += heights[j]*exp(-(x-centers[j])*(x-centers[j])/(2*sigma*sigma));
        if (i > 0) {
            centers.push_back(x);
            heights.push_back(height0*exp(-bias/(delta_temperature*BOLTZ)));
        }
        ASSERT(fabs(bias + x*x - state.getPotentialEnergy()) < 0.1);
    }
}

void testMassesCharges(Platform& platform) {

    // Create a system with one paticle
    System system;
    system.addParticle(3.8); // Set mass
    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->addParticle(-2.1, 0.0, 0.0); // Set charge
    system.addForce(nonbonded);

    // Setup PLUMED to write the mass and chage of the particles to a file
    const string script = "DUMPMASSCHARGE ATOMS=@mdatoms FILE=mass_charge.txt";
    PlumedForce* plumed = new PlumedForce(script);
    plumed->setLogStream(logstream);
    system.addForce(plumed);

    // Setup simulation
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Context context(system, integ, platform);

    ifstream stream;
    char header[100];
    double _, mass, charge;

    context.setPositions({Vec3()});
    integ.step(2); // Need at least 2 step for dumping to work

    // Parse the dumped file
    stream.open("mass_charge.txt");
    stream.getline(&header[0], 100);
    stream >> _ >> mass >> charge;
    stream.close();

    // Chekc if the mass and change from System is used
    ASSERT_EQUAL(mass, 3.8);
    ASSERT_EQUAL(charge, -2.1);

    // Set the PLUMED masses
    plumed->setMasses({7.5});
    context.reinitialize(true);
    integ.step(2); // Need at least 2 step for dumping

    // Parse the dumped file
    stream.open("mass_charge.txt");
    stream.getline(&header[0], 100);
    stream >> _ >> mass >> charge;
    stream.close();

    // Chekc if the mass from PLUMED is used
    ASSERT_EQUAL(mass, 7.5);
    ASSERT_EQUAL(charge, -2.1);

    // Reset the PLUMED masses
    plumed->setMasses({});
    context.reinitialize(true);
    integ.step(2); // Need at least 2 step for dumping

    // Parse the dumped file
    stream.open("mass_charge.txt");
    stream.getline(&header[0], 100);
    stream >> _ >> mass >> charge;
    stream.close();

    // Chekc if the mass and change from System is used again
    ASSERT_EQUAL(mass, 3.8);
    ASSERT_EQUAL(charge, -2.1);
}

void testScript(Platform& platform) {

    // Create a system
    System system;
    system.addParticle(0);

    // Setup PLUMED
    const string script = "#Comment\n"
                          "p: POSITION ATOM=1\n"
                          "\n"
                          "# More comments and empty lines\n"
                          "PRINT ...\n"
                          "\n"
                          "  ARG=p.x,p.y,p.z\n"
                          "  # A comment in the middle\n"
                          "  STRIDE=10\n"
                          "...";
    PlumedForce* plumed = new PlumedForce(script);
    plumed->setLogStream(logstream);
    system.addForce(plumed);

    // Setup simulation
    LangevinIntegrator integ(300.0, 1.0, 1.0);
    Context context(system, integ, platform);

    // If the parser fails, an exception is thrown during the context creation
}

void testPlatform(Platform& platform) {
    testForce(platform);
    testMetadynamics(platform);
    testWellTemperedMetadynamics(platform);
    testMassesCharges(platform);
    testScript(platform);
}

int main() {
    try {
        logstream = fopen("/dev/null", "w");
        Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
        for (int i = 0; i < Platform::getNumPlatforms(); i++) {
            Platform& platform = Platform::getPlatform(i);
            printf("Testing %s\n", platform.getName().c_str());
            testPlatform(platform);
        }
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
