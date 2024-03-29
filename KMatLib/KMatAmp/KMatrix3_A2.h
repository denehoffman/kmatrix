#if !defined(KMATRIX3_A2)
#define KMATRIX3_A2

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/AmpParameter.h"
#include "IUAmpTools/UserAmplitude.h"
#include "GPUManager/GPUCustomTypes.h"

#include "Math/SVector.h"
#include "Math/SMatrix.h"

#include <utility>
#include <string>
#include <complex>
#include <vector>

using std::complex;
using namespace std;
using namespace ROOT::Math;
typedef SMatrix<complex<GDouble>, 3> SMatrix3;
typedef SMatrix<complex<GDouble>, 3, 3, ROOT::Math::MatRepSym<complex<GDouble>, 3>> SMatrix3Sym;
typedef SVector<complex<GDouble>, 3> SVector3;


class Kinematics;

class KMatrix3_A2: public UserAmplitude<KMatrix3_A2> {
    public:
        KMatrix3_A2(): UserAmplitude<KMatrix3_A2>() {}
        KMatrix3_A2(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re,
            k0im, k1im, k2im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KMatrix3_A2(){}
    
        string name() const {return "KMatrix3_A2";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter ba21320_re;
        AmpParameter ba21320_im;
        AmpParameter ba21700_re;
        AmpParameter ba21700_im;
        // Channel order is PiEta, KK, PiEta'
        SVector3 ga21320;
        SVector3 ga21700;
        vector<SVector3> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix3 mat_bkg;

        SMatrix3 inverse3(SMatrix3 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;

};
#endif
