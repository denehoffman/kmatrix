#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <math.h>

#include "TLorentzVector.h"
#include "barrierFactor.h"
#include "breakupMomentum.h"

#include "IUAmpTools/Kinematics.h"
#include "KMatrix5_F0.h"


KMatrix5_F0::KMatrix5_F0(const vector<string> &args): UserAmplitude<KMatrix5_F0>(args) {
    /*
     * Usage: KMatrix5_F0 <daughter 1> <daughter 2> <Re[f_0(500)]> <Im[f_0(500)]> ...
     */
	m_daughters = pair<string, string>(args[0], args[1]);
	bf0500_re = AmpParameter(args[2]);
	bf0500_im = AmpParameter(args[3]);
    bf0980_re = AmpParameter(args[4]);
    bf0980_im = AmpParameter(args[5]);
    bf01370_re = AmpParameter(args[6]);
    bf01370_im = AmpParameter(args[7]);
    bf01500_re = AmpParameter(args[8]);
    bf01500_im = AmpParameter(args[9]);
    bf01710_re = AmpParameter(args[10]);
    bf01710_im = AmpParameter(args[11]);
    registerParameter(bf0500_re);
    registerParameter(bf0500_im);
    registerParameter(bf0980_re);
    registerParameter(bf0980_im);
    registerParameter(bf01370_re);
    registerParameter(bf01370_im);
    registerParameter(bf01500_re);
    registerParameter(bf01500_im);
    registerParameter(bf01710_re);
    registerParameter(bf01710_im);
    gf0500 = SVector5(0.74987, -0.01257, 0.27536, -0.15102, 0.36103);
    gf0980 = SVector5(0.06401, 0.00204, 0.77413, 0.50999, 0.13112);
    gf01370 = SVector5(-0.23417, -0.01032, 0.72283, 0.11934, 0.36792);
    gf01500 = SVector5(0.01270, 0.26700, 0.09214, 0.02742, -0.04025);
    gf01710 = SVector5(-0.14242, 0.22780, 0.15981, 0.16272, -0.17397);
    masses = {0.51461, 0.90630, 1.23089, 1.46104, 1.69611};
    couplings = {gf0500, gf0980, gf01370, gf01500, gf01710};
    m1s = {0.1349768, 2*0.1349768, 0.493677, 0.547862, 0.547862};
    m2s = {0.1349768, 2*0.1349768, 0.497611, 0.547862, 0.95778};
    a_bkg = {
        0.03728,
        0.00000, 0.00000,
        -0.01398, 0.00000, 0.02349,
        -0.02203, 0.00000, 0.03101, -0.13769,
        0.01397, 0.00000, -0.04003, -0.06722, -0.28401};
    mat_bkg = SMatrix5Sym(a_bkg.begin(), a_bkg.end());
}


void KMatrix5_F0::calcUserVars(GDouble** pKin, GDouble* userVars) const {
    TLorentzVector pTemp, pTot;
    //TLorentzVector p1, p2;
    /* This allows us to input something like
     * "12 34" for the daughter particle parameters
     * to indicate that the first daughter is the
     * sum of the final state particles indexed
     * 1 and 2 and the second daughter is created
     * from the sum of 3 and 4
     */
    for(unsigned int i=0; i < m_daughters.first.size(); ++i) {
        string num; num += m_daughters.first[i];
        int index = atoi(num.c_str());
        pTemp.SetPxPyPzE(pKin[index][1],
                         pKin[index][2],
                         pKin[index][3],
                         pKin[index][0]);
        //p1 += pTemp;
        pTot += pTemp;
    }
    for(unsigned int i=0; i < m_daughters.second.size(); ++i) {
        string num; num += m_daughters.second[i];
        int index = atoi(num.c_str());
        pTemp.SetPxPyPzE(pKin[index][1],
                         pKin[index][2],
                         pKin[index][3],
                         pKin[index][0]);
        //p2 += pTemp;
        pTot += pTemp;
    }
    // Grab our hypothesis mass
    GDouble m = pTot.M();
    userVars[kM] = m;
    GDouble s = pTot.M2();
    userVars[kS] = s;
    // Calculate K-Matrix
    // Initialize K-Matrix
    SMatrix5 mat_K; // Initialized as a 5x5 0-matrix
    SMatrix5 mat_C;
    // Loop over resonances
    for (int i = 0; i < 5; i++) {
        SMatrix5 temp_K; // Initialized as a 5x5 0-matrix
        SMatrix5 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_K /= denominator;
        temp_K += mat_bkg;
        // Loop over channels:
        SVector5 B_factor;
        for (int j = 0; j < 5; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 0);
            GDouble B_den = barrierFactor(q_alpha, 0);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        for (int k = 0; k < 5; k++) {
            for (int l = 0; l < 5; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // Needed element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 5; j++) {
        mat_C(j, j) = KMatrix5_F0::chew(s, m1s[j], m2s[j]);
    }
    // Now mat_K should be done (ignore (s-s_0)/s_norm for now)
    SMatrix5 temp = SMatrixIdentity();
    mat_K *= mat_C;
    mat_K *= ((s - 0.0091125) / 1); // Adler zero term
    temp += mat_K;
    // Now temp is the stuff that we want to invert before multiplying by the P-vector
    SMatrix5 temp_inv = KMatrix5_F0::inverse5(temp);
    // Now we cache the results
    // Note that it doesn't matter if I store it (i, j) or (j, i) as
    // long as I'm consistent later when unpacking it...
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            userVars[i * 5 + j + 2] = temp_inv(i, j).real(); // +2 because kM and kS are first in the enum
            userVars[i * 5 + j + 2 + 25] = temp_inv(i, j).imag(); // +25 to skip the real parts 
        }
    }
}

complex<GDouble> KMatrix5_F0::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector5 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(bf0500_re, bf0500_im),
        complex<GDouble>(bf0980_re, bf0980_im),
        complex<GDouble>(bf01370_re, bf01370_im),
        complex<GDouble>(bf01500_re, bf01500_im),
        complex<GDouble>(bf01710_re, bf01710_im)
    };
    // Loop over resonances
    for (int i = 0; i < 5; i++) {
        SVector5 temp_P;
        SMatrix5 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_P /= denominator;
        // Loop over channels:
        SVector5 B_factor;
        for (int j = 0; j < 5; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 0);
            GDouble B_den = barrierFactor(q_alpha, 0);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        temp_P = temp_P * B_factor;
        vec_P += temp_P;
    }
    SMatrix5 temp_inv;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            temp_inv(i, j) = complex<GDouble>(userVars[i * 5 + j + 2],
                                              userVars[i * 5 + j + 2 + 25]);
        }
    }
    SVector5 res = temp_inv * vec_P;
    return res[2]; // return the KK channel contribution
}

SMatrix5 KMatrix5_F0::inverse5(SMatrix5 A) const {
    SMatrix5 M;
    SMatrix5 I = SMatrixIdentity();
    SMatrix5 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 5; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KMatrix5_F0::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KMatrix5_F0::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KMatrix5_F0::chew(double s, double m1, double m2) const {
    complex<GDouble> tot = 0;
    tot += (KMatrix5_F0::rho(s, m1, m2) / Pi()) * log((KMatrix5_F0::xi(s, m1, m2) + KMatrix5_F0::rho(s, m1, m2)) / (KMatrix5_F0::xi(s, m1, m2) - KMatrix5_F0::rho(s, m1, m2)));
    tot -= (KMatrix5_F0::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

void KMatrix5_F0::updatePar(const AmpParameter &par) {}
