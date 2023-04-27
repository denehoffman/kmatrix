#if !defined(PHASESPACE)
#define PHASESPACE

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>


using std::complex;
using namespace std;

class Kinematics;

class PhaseSpace: public UserAmplitude<PhaseSpace> {  
    public:

        PhaseSpace(): UserAmplitude<PhaseSpace>() {}
        PhaseSpace(const vector<string>& args);

        enum UserVars {kM = 0, kRhoReal, kRhoImag, kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~PhaseSpace(){}

        string name() const {return "PhaseSpace";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:

        pair<string, string> m_daughters;
        AmpParameter m_slope;
        
        complex<GDouble> rho(double s, double m1, double m2) const;
};
#endif
