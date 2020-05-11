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

#include "CudaPlumedKernels.h"
#include "CudaPlumedKernelSources.h"
#include "openmm/NonbondedForce.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include <cstring>
#include <map>

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

class CudaCalcPlumedForceKernel::StartCalculationPreComputation : public CudaContext::ForcePreComputation {
public:
    StartCalculationPreComputation(CudaCalcPlumedForceKernel& owner) : owner(owner) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.beginComputation(includeForces, includeEnergy, groups);
    }
    CudaCalcPlumedForceKernel& owner;
};

class CudaCalcPlumedForceKernel::ExecuteTask : public CudaContext::WorkTask {
public:
    ExecuteTask(CudaCalcPlumedForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.executeOnWorkerThread();
    }
    CudaCalcPlumedForceKernel& owner;
};

class CudaCalcPlumedForceKernel::CopyForcesTask : public ThreadPool::Task {
public:
    CopyForcesTask(CudaContext& cu, vector<Vec3>& forces) : cu(cu), forces(forces) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        // Copy the forces applied by PLUMED to a buffer for uploading.  This is done in parallel for speed.
        
        int numParticles = cu.getNumAtoms();
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (cu.getUseDoublePrecision()) {
            double* buffer = (double*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                const Vec3& p = forces[i];
                buffer[3*i] = p[0];
                buffer[3*i+1] = p[1];
                buffer[3*i+2] = p[2];
            }
        }
        else {
            float* buffer = (float*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                const Vec3& p = forces[i];
                buffer[3*i] = (float) p[0];
                buffer[3*i+1] = (float) p[1];
                buffer[3*i+2] = (float) p[2];
            }
        }
    }
    CudaContext& cu;
    vector<Vec3>& forces;
};

class CudaCalcPlumedForceKernel::AddForcesPostComputation : public CudaContext::ForcePostComputation {
public:
    AddForcesPostComputation(CudaCalcPlumedForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return owner.addForces(includeForces, includeEnergy, groups);
    }
    CudaCalcPlumedForceKernel& owner;
};

CudaCalcPlumedForceKernel::~CudaCalcPlumedForceKernel() {
    cu.setAsCurrent();
    if (plumedForces != NULL)
        delete plumedForces;
    cuStreamDestroy(stream);
    cuEventDestroy(syncEvent);
}

void CudaCalcPlumedForceKernel::initialize(const System& system, const PlumedForce& force) {
    cu.setAsCurrent();
    cuStreamCreate(&stream, CU_STREAM_NON_BLOCKING);
    cuEventCreate(&syncEvent, CU_EVENT_DISABLE_TIMING);
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    plumedForces = new CudaArray(cu, 3*system.getNumParticles(), elementSize, "plumedForces");
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaPlumedKernelSources::plumedForce, defines);
    addForcesKernel = cu.getKernel(module, "addForces");
    forceGroupFlag = (1<<force.getForceGroup());
    cu.addPreComputation(new StartCalculationPreComputation(*this));
    cu.addPostComputation(new AddForcesPostComputation(*this));

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
    plumed_cmd(plumedmain, "setLog", force.getLogStream());
    int numParticles = system.getNumParticles();
    plumed_cmd(plumedmain, "setNatoms", &numParticles);
    double dt = contextImpl.getIntegrator().getStepSize();
    plumed_cmd(plumedmain, "setTimestep", &dt);
    double kT = force.getKbT();
    if (kT >= 0.0)
        plumed_cmd(plumedmain,"setKbT",&kT);
    int restart = force.getRestart();
    plumed_cmd(plumedmain, "setRestart", &restart);
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
    // This method does nothing.  The actual calculation is started by the pre-computation, continued on
    // the worker thread, and finished by the post-computation.
    
    return 0;
}

void CudaCalcPlumedForceKernel::beginComputation(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return;
    contextImpl.getPositions(positions);
    
    // The actual force computation will be done on a different thread.
    
    cu.getWorkThread().addTask(new ExecuteTask(*this));
}

void CudaCalcPlumedForceKernel::executeOnWorkerThread() {
    // Configure the PLUMED interface object.
    
    int numParticles = contextImpl.getSystem().getNumParticles();
    int step = cu.getStepCount();
    plumed_cmd(plumedmain, "setStep", &step);
    plumed_cmd(plumedmain, "setMasses", &masses[0]);
    if (charges.size() > 0)
        plumed_cmd(plumedmain, "setCharges", &charges[0]);
    plumed_cmd(plumedmain, "setPositions", &positions[0][0]);
    forces.resize(numParticles);
    memset(&forces[0], 0, numParticles*sizeof(Vec3));
    plumed_cmd(plumedmain, "setForces", &forces[0][0]);
    if (usesPeriodic) {
        Vec3 boxVectors[3];
        contextImpl.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
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
    
    // Upload the forces to the device.
    
    CopyForcesTask task(cu, forces);
    cu.getPlatformData().threads.execute(task);
    cu.getPlatformData().threads.waitForThreads();
    cu.setAsCurrent();
    cuMemcpyHtoDAsync(plumedForces->getDevicePointer(), cu.getPinnedBuffer(), plumedForces->getSize()*plumedForces->getElementSize(), stream);
    cuEventRecord(syncEvent, stream);
    
}

double CudaCalcPlumedForceKernel::addForces(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return 0;

    // Wait until executeOnWorkerThread() is finished.
    
    cu.getWorkThread().flush();
    cuStreamWaitEvent(cu.getCurrentStream(), syncEvent, 0);

    // Add in the forces.
    
    if (includeForces) {
        void* args[] = {&plumedForces->getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer()};
        cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    }
    
    // Return the energy.
    
    double energy = 0;
    plumed_cmd(plumedmain, "getBias", &energy);
    return energy;
}
