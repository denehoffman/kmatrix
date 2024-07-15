#if !defined(KMATRIX5_F0)
#define KMATRIX5_F0

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
typedef SMatrix<complex<GDouble>, 5> SMatrix5;
typedef SMatrix<complex<GDouble>, 5, 5, ROOT::Math::MatRepSym<complex<GDouble>, 5>> SMatrix5Sym;
typedef SVector<complex<GDouble>, 5> SVector5;


class Kinematics;

class KMatrix5_F0: public UserAmplitude<KMatrix5_F0> {
    public:
        KMatrix5_F0(): UserAmplitude<KMatrix5_F0>() {}
        KMatrix5_F0(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re, k3re, k4re,
            k0im, k1im, k2im, k3im, k4im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}
        
        ~KMatrix5_F0(){}
    
        string name() const {return "KMatrix5_F0";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter bf0500_re;
        AmpParameter bf0500_im;
        AmpParameter bf0980_re;
        AmpParameter bf0980_im;
        AmpParameter bf01370_re;
        AmpParameter bf01370_im;
        AmpParameter bf01500_re;
        AmpParameter bf01500_im;
        AmpParameter bf01710_re;
        AmpParameter bf01710_im;

        // Channel order is PiPi, 2Pi2Pi, KK, EtaEta, EtaEta'
        SVector5 gf0500;
        SVector5 gf0980;
        SVector5 gf01370;
        SVector5 gf01500;
        SVector5 gf01710;
        vector<SVector5> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix5 mat_bkg;

        SMatrix5 inverse5(SMatrix5 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;
        complex<GDouble> poleProduct(double s) const;
        complex<GDouble> poleProductRemainder(double s, size_t index) const;
};
#endif
