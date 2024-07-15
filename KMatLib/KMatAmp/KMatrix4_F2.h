#if !defined(KMATRIX4_F2)
#define KMATRIX4_F2

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
typedef SMatrix<complex<GDouble>, 4> SMatrix4;
typedef SMatrix<complex<GDouble>, 4, 4, ROOT::Math::MatRepSym<complex<GDouble>, 4>> SMatrix4Sym;
typedef SVector<complex<GDouble>, 4> SVector4;


class Kinematics;

class KMatrix4_F2: public UserAmplitude<KMatrix4_F2> {
    public:
        KMatrix4_F2(): UserAmplitude<KMatrix4_F2>() {}
        KMatrix4_F2(const vector<string> &args);
	    enum UserVars {kM = 0, kS,
            k0re, k1re, k2re, k3re,
            k0im, k1im, k2im, k3im,
            kNumUserVars};
        unsigned int numUserVars() const {return kNumUserVars;}
        void calcUserVars(GDouble** pKin, GDouble* userVars) const;
        bool needsUserVarsOnly() const {return true;}
        bool areUserVarsStatic() const {return true;}

        ~KMatrix4_F2(){}
    
        string name() const {return "KMatrix4_F2";}

        complex<GDouble> calcAmplitude(GDouble** pKin, GDouble* userVars) const;
        void updatePar(const AmpParameter &par);

    private:
        pair<string, string> m_daughters;
        int channel;
        AmpParameter bf21270_re;
        AmpParameter bf21270_im;
        AmpParameter bf21525_re;
        AmpParameter bf21525_im;
        AmpParameter bf21810_re;
        AmpParameter bf21810_im;
        AmpParameter bf21950_re;
        AmpParameter bf21950_im;
        // Channel order is PiPi, 2Pi2Pi, KK, EtaEta
        SVector4 gf21270;
        SVector4 gf21525;
        SVector4 gf21810;
        SVector4 gf21950;
        vector<SVector4> couplings;
        vector<GDouble> masses;
        vector<GDouble> m1s;
        vector<GDouble> m2s;
        vector<GDouble> a_bkg;
        SMatrix4 mat_bkg;

        SMatrix4 inverse4(SMatrix4 mat) const;
        complex<GDouble> rho(double s, double m1, double m2) const;
        complex<GDouble> xi(double s, double m1, double m2) const;
        complex<GDouble> chew(double s, double m1, double m2) const;
        complex<GDouble> poleProduct(double s) const;
        complex<GDouble> poleProductRemainder(double s, size_t index) const;
};
#endif
