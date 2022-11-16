#if !defined(KMATRIX2_A0)
#define KMATRIX2_A0

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
typedef SMatrix<complex<GDouble>, 2> SMatrix2;
typedef SMatrix<complex<GDouble>, 2, 2, ROOT::Math::MatRepSym<complex<GDouble>, 2>> SMatrix2Sym;
typedef SVector<complex<GDouble>, 2> SVector2;


class Kinematics;

class KMatrix2_A0: public UserAmplitude<KMatrix2_A0> {
    public:
        KMatrix2_A0(): UserAmplitude<KMatrix2_A0>() {}
        KMatrix2_A0(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k00re, k01re,
            k10re, k11re,
            k00im, k01im,
            k10im, k11im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KMatrix2_A0(){}
    
        string name() const {return "KMatrix2_A0";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        AmpParameter ba0980_re;
        AmpParameter ba0980_im;
        AmpParameter ba01450_re;
        AmpParameter ba01450_im;
        // Channel order is PiEta, KK
        SVector2 ga0980;
        SVector2 ga01450;
        vector<SVector2> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix2 mat_bkg;

        SMatrix2 inverse2(SMatrix2 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;

};
#endif
