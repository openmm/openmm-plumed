/* -------------------------------------------------------------------------- *
 *                                OpenMMPlumed                                 *
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

#include "PlumedForceProxy.h"
#include "PlumedForce.h"
#include "openmm/serialization/SerializationNode.h"

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

PlumedForceProxy::PlumedForceProxy() : SerializationProxy("PlumedForce") {
}

void PlumedForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 4);
    const PlumedForce& force = *reinterpret_cast<const PlumedForce*>(object);
    node.setStringProperty("script", force.getScript());
    node.setDoubleProperty("temperature", force.getTemperature());
    auto& particles = node.createChildNode("particles");
    for (const auto& mass: force.getMasses())
        particles.createChildNode("particle").setDoubleProperty("mass", mass);
    node.setBoolProperty("restart", force.getRestart());
}

void* PlumedForceProxy::deserialize(const SerializationNode& node) const {
    const int version = node.getIntProperty("version");
    if (version < 1 || version > 4)
        throw OpenMMException("Unsupported version number");

    PlumedForce* force = new PlumedForce(node.getStringProperty("script"));
    if (version > 1)
        force->setRestart(node.getBoolProperty("restart"));
    if (version > 2)
        force->setTemperature(node.getDoubleProperty("temperature"));
    if (version > 3) {
        std::vector<double> masses;
        for (const auto& particle: node.getChildNode("particles").getChildren())
            masses.push_back(particle.getDoubleProperty("mass"));
        force->setMasses(masses);
    }

    return force;
}
