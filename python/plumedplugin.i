%module openmmplumed

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

%{
#include "PlumedForce.h"
#include "OpenMM.h"
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
