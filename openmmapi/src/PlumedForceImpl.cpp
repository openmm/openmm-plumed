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

#include "internal/PlumedForceImpl.h"
#include "PlumedKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/NonbondedForce.h"
#include "openmm/reference/SimTKOpenMMRealType.h"

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

PlumedForceImpl::PlumedForceImpl(const PlumedForce& owner) : CustomCPPForceImpl(owner), owner(owner), hasInitialized(false) {
}

PlumedForceImpl::~PlumedForceImpl() {
    if (hasInitialized)
        plumed_finalize(plumedmain);
}

void PlumedForceImpl::initialize(ContextImpl& context) {
    CustomCPPForceImpl::initialize(context);
    const OpenMM::System& system = context.getSystem();

    // Construct and initialize the PLUMED interface object.

    plumedmain = plumed_create();
    hasInitialized = true;
    int apiVersion;
    plumed_cmd(plumedmain, "getApiVersion", &apiVersion);
    if (apiVersion < 4)
        throw OpenMMException("Unsupported API version.  Upgrade PLUMED to a newer version.");
    int precision = 8;
    plumed_cmd(plumedmain, "setRealPrecision", &precision);
    double conversion = 1.0;
    plumed_cmd(plumedmain, "setMDEnergyUnits", &conversion);
    plumed_cmd(plumedmain, "setMDLengthUnits", &conversion);
    plumed_cmd(plumedmain, "setMDTimeUnits", &conversion);
    plumed_cmd(plumedmain, "setMDEngine", "OpenMM");
    plumed_cmd(plumedmain, "setLog", owner.getLogStream());
    int numParticles = system.getNumParticles();
    plumed_cmd(plumedmain, "setNatoms", &numParticles);
    double dt = context.getIntegrator().getStepSize();
    plumed_cmd(plumedmain, "setTimestep", &dt);
    double kT = owner.getTemperature() * BOLTZ;
    if (kT >= 0.0)
        plumed_cmd(plumedmain, "setKbT", &kT);
    int restart = owner.getRestart();
    plumed_cmd(plumedmain, "setRestart", &restart);
    plumed_cmd(plumedmain, "init", NULL);
    if(apiVersion > 7) {
        plumed_cmd(plumedmain, "readInputLines", owner.getScript().c_str());
    } else {
        // NOTE: the comments and line continuation does not works
        //       (https://github.com/plumed/plumed2/issues/571)
        // TODO: remove this when PLUMED 2.6 support is dropped
        vector<char> scriptChars(owner.getScript().size()+1);
        strcpy(&scriptChars[0], owner.getScript().c_str());
        char* line = strtok(&scriptChars[0], "\r\n");
        while (line != NULL) {
            plumed_cmd(plumedmain, "readInputLine", line);
            line = strtok(NULL, "\r\n");
        }
    }
    usesPeriodic = system.usesPeriodicBoundaryConditions();

    // Record the particle masses.

    masses.resize(numParticles);
    const auto& plumedMasses = owner.getMasses();
    if (plumedMasses.size() == 0) // User System masses
        for (int i = 0; i < numParticles; i++)
            masses[i] = system.getParticleMass(i);
    else if (plumedMasses.size() == numParticles) // User PLUMED masses
        masses = plumedMasses;
    else
        throw OpenMMException("The number of PLUMED masses is different from the number of particles!");

    // If there's a NonbondedForce, get charges from it.

    for (int j = 0; j < system.getNumForces(); j++) {
        const NonbondedForce* nonbonded = dynamic_cast<const NonbondedForce*>(&system.getForce(j));
        if (nonbonded != NULL) {
            charges.resize(numParticles);
            double sigma, epsilon;
            for (int i = 0; i < numParticles; i++)
                nonbonded->getParticleParameters(i, charges[i], sigma, epsilon);
        }
    }
}

double PlumedForceImpl::computeForce(ContextImpl& context, const vector<Vec3>& positions, vector<Vec3>& forces) {
    // Pass the current state to PLUMED.

    int step = context.getStepCount();
    plumed_cmd(plumedmain, "setStep", &step);
    plumed_cmd(plumedmain, "setMasses", masses.data());
    if (charges.size() > 0)
        plumed_cmd(plumedmain, "setCharges", charges.data());
    vector<Vec3>& pos = const_cast<vector<Vec3>&>(positions);
    plumed_cmd(plumedmain, "setPositions", &pos[0][0]);
    plumed_cmd(plumedmain, "setForces", &forces[0][0]);
    if (usesPeriodic) {
        Vec3 boxVectors[3];
        context.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        plumed_cmd(plumedmain, "setBox", &boxVectors[0][0]);
    }
    vector<double> virial(9, 0.0);
    plumed_cmd(plumedmain, "setVirial", virial.data());

    // Calculate the forces and energy.

    memset(&forces[0], 0, forces.size()*sizeof(Vec3));
    plumed_cmd(plumedmain, "prepareCalc", NULL);
    plumed_cmd(plumedmain, "performCalcNoUpdate", NULL);
    if (step != lastStepIndex) {
        plumed_cmd(plumedmain, "update", NULL);
        lastStepIndex = step;
    }
    double energy = 0;
    plumed_cmd(plumedmain, "getBias", &energy);
    return energy;
}
