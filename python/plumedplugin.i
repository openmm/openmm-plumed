%module openmmplumed

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"
%include "std_vector.i"

namespace std {
  %template(vectord) vector<double>;
}

%{
#include "PlumedForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include <numpy/arrayobject.h>

int isNumpyAvailable() {
    return true;
}
%}

%pythoncode %{
import openmm as mm
%}

namespace PlumedPlugin {

class PlumedForce : public OpenMM::Force {
public:
    PlumedForce(const std::string& script);
    const std::string& getScript() const;
    bool usesPeriodicBoundaryConditions() const;
    void setTemperature(double temperature);
    double getTemperature() const;
    void setMasses(const std::vector<double>& masses);
    const std::vector<double>& getMasses() const;
    void setRestart(bool restart);
    bool getRestart() const;
};

}
