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

#include "PlumedForce.h"
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

extern "C" void registerPlumedSerializationProxies();

void testSerialization() {
    // Create a Force.

    string script = "d: DISTANCE ATOMS=1,3\n"
                    "BIASVALUE ARG=d";
    bool restart = true;
    double temperature = 42.0;
    const std::vector<double> masses = {3.1, 4.1, 5.9};
    PlumedForce force(script);
    force.setRestart(restart);
    force.setTemperature(temperature);
    force.setMasses(masses);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<PlumedForce>(&force, "Force", buffer);
    PlumedForce* copy = XmlSerializer::deserialize<PlumedForce>(buffer);

    // Compare the two forces to see if they are identical.

    PlumedForce& force2 = *copy;
    ASSERT_EQUAL(script, force2.getScript());
    ASSERT_EQUAL(restart, force2.getRestart());
    ASSERT_EQUAL(temperature, force2.getTemperature());
    ASSERT_EQUAL_CONTAINERS(masses, force2.getMasses());
}

int main() {
    try {
        registerPlumedSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
