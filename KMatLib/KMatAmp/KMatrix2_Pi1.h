#if !defined(KMATRIX2_PI1)
#define KMATRIX2_PI1

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

class KMatrix2_Pi1: public UserAmplitude<KMatrix2_Pi1> {
    public:
        KMatrix2_Pi1(): UserAmplitude<KMatrix2_Pi1>() {}
        KMatrix2_Pi1(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re,
            k0im, k1im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KMatrix2_Pi1(){}
    
        string name() const {return "KMatrix2_Pi1";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter bpi11600_re;
        AmpParameter bpi11600_im;
        // Channel order is PiEta, PiEta'
        SVector2 gpi11600;
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
