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

#include "CudaPlumedKernels.h"
#include "CudaPlumedKernelSources.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include <cstring>
#include <map>

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

CudaCalcPlumedForceKernel::~CudaCalcPlumedForceKernel() {
    cu.setAsCurrent();
    if (plumedForces != NULL)
        delete plumedForces;
}

void CudaCalcPlumedForceKernel::initialize(const System& system, const PlumedForce& force) {
    cu.setAsCurrent();
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    plumedForces = new CudaArray(cu, 3*system.getNumParticles(), elementSize, "plumedForces");
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaPlumedKernelSources::plumedForce, defines);
    addForcesKernel = cu.getKernel(module, "addForces");

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
    int numParticles = system.getNumParticles();
    plumed_cmd(plumedmain, "setNatoms", &numParticles);
    double dt = contextImpl.getIntegrator().getStepSize();
    plumed_cmd(plumedmain, "setTimestep", &dt);
    plumed_cmd(plumedmain, "init", NULL);
    vector<char> scriptChars(force.getScript().size()+1);
    strcpy(&scriptChars[0], force.getScript().c_str());
    char* line = strtok(&scriptChars[0], "\r\n");
    while (line != NULL) {
        plumed_cmd(plumedmain, "readInputLine", line);
        line = strtok(NULL, "\r\n");
    }
    usesPeriodic = system.usesPeriodicBoundaryConditions();

    // Record the particle masses.

    masses.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
        masses[i] = system.getParticleMass(i);

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

double CudaCalcPlumedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    int numParticles = context.getSystem().getNumParticles();
    int step = cu.getStepCount();
    plumed_cmd(plumedmain, "setStep", &step);
    plumed_cmd(plumedmain, "setMasses", &masses[0]);
    if (charges.size() > 0)
        plumed_cmd(plumedmain, "setCharges", &charges[0]);
    context.getPositions(positions);
    plumed_cmd(plumedmain, "setPositions", &positions[0][0]);
    forces.resize(numParticles);
    Vec3 zero;
    for (int i = 0; i < numParticles; i++)
        forces[i] = zero;
    plumed_cmd(plumedmain, "setForces", &forces[0][0]);
    if (usesPeriodic) {
        Vec3 boxVectors[3];
        context.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        plumed_cmd(plumedmain, "setBox", &boxVectors[0][0]);
    }
    double virial[9];
    plumed_cmd(plumedmain, "setVirial", &virial);

    // Calculate the forces and energy.

    plumed_cmd(plumedmain, "prepareCalc", NULL);
    plumed_cmd(plumedmain, "performCalcNoUpdate", NULL);
    if (step != lastStepIndex) {
        plumed_cmd(plumedmain, "update", NULL);
        lastStepIndex = step;
    }
    double energy = 0;
    plumed_cmd(plumedmain, "getBias", &energy);
    
    // Upload the forces to the device and add them in.
    
    if (cu.getUseDoublePrecision()) {
        double* buffer = (double*) cu.getPinnedBuffer();
        for (int i = 0; i < numParticles; ++i) {
            const Vec3& p = forces[i];
            buffer[3*i] = p[0];
            buffer[3*i+1] = p[1];
            buffer[3*i+2] = p[2];
        }
        plumedForces->upload(buffer);
    }
    else {
        float* buffer = (float*) cu.getPinnedBuffer();
        for (int i = 0; i < numParticles; ++i) {
            const Vec3& p = forces[i];
            buffer[3*i] = (float) p[0];
            buffer[3*i+1] = (float) p[1];
            buffer[3*i+2] = (float) p[2];
        }
        plumedForces->upload(buffer);
    }
    void* args[] = {&plumedForces->getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer()};
    cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    return energy;
}

