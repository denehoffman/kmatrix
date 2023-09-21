#if !defined(KMATRIX5_F0B)
#define KMATRIX5_F0B

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

class KMatrix5_F0b: public UserAmplitude<KMatrix5_F0b> {
    public:
        KMatrix5_F0b(): UserAmplitude<KMatrix5_F0b>() {}
        KMatrix5_F0b(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k00re, k01re, k02re, k03re, k04re,
            k10re, k11re, k12re, k13re, k14re,
            k20re, k21re, k22re, k23re, k24re,
            k30re, k31re, k32re, k33re, k34re,
            k40re, k41re, k42re, k43re, k44re,
            k00im, k01im, k02im, k03im, k04im,
            k10im, k11im, k12im, k13im, k14im,
            k20im, k21im, k22im, k23im, k24im,
            k30im, k31im, k32im, k33im, k34im,
            k40im, k41im, k42im, k43im, k44im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}
        
        ~KMatrix5_F0b(){}
    
        string name() const {return "KMatrix5_F0b";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
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
        AmpParameter c0_re;
        AmpParameter c0_im;

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
};
#endif
