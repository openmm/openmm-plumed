%module openmmplumed

%include "factory.i"
%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"

%{
#include "PlumedForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
%}

namespace PlumedPlugin {

class PlumedForce : public OpenMM::Force {
public:
    PlumedForce(const std::string& script);
    const std::string& getScript() const;
};

}
