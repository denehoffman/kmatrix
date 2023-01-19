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
#include "KMatrix4_F2.h"


KMatrix4_F2::KMatrix4_F2(const vector<string> &args): UserAmplitude<KMatrix4_F2>(args) {
    /*
     * Usage: KMatrix4_F2 <daughter 1> <daughter 2> <Re[f_0(400)]> <Im[f_0(400)]> ...
     */
	m_daughters = pair<string, string>(args[0], args[1]);
	bf21270_re = AmpParameter(args[2]);
	bf21270_im = AmpParameter(args[3]);
    bf21525_re = AmpParameter(args[4]);
    bf21525_im = AmpParameter(args[5]);
    bf21810_re = AmpParameter(args[6]);
    bf21810_im = AmpParameter(args[7]);
    bf21950_re = AmpParameter(args[8]);
    bf21950_im = AmpParameter(args[9]);
    registerParameter(bf21270_re);
    registerParameter(bf21270_im);
    registerParameter(bf21525_re);
    registerParameter(bf21525_im);
    registerParameter(bf21810_re);
    registerParameter(bf21810_im);
    registerParameter(bf21950_re);
    registerParameter(bf21950_im);
    gf21270 = SVector4(0.40033, 0.15479, -0.08900, -0.00113);
    gf21525 = SVector4(0.01820, 0.17300, 0.32393, 0.15256);
    gf21810 = SVector4(-0.06709, 0.22941, -0.43133, 0.23721);
    gf21950 = SVector4(-0.49924, 0.19295, 0.27975, -0.03987);
    masses = {1.15299, 1.48359, 1.72923, 1.96700};
    couplings = {gf21270, gf21525, gf21810, gf21950};
    m1s = {0.1349768, 2*0.1349768, 0.493677, 0.547862};
    m2s = {0.1349768, 2*0.1349768, 0.497611, 0.547862};
    a_bkg = {
        -0.04139,
        0.00000, 0.00000,
        0.00984, 0.00000, -0.07344,
        0.01028, 0.00000, 0.05533, -0.05183};
    mat_bkg = SMatrix4Sym(a_bkg.begin(), a_bkg.end());
}

void KMatrix4_F2::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
    SMatrix4 mat_K; // Initialized as a 4x4 0-matrix
    SMatrix4 mat_C;
    // Loop over resonances
    for (int i = 0; i < 4; i++) {
        SMatrix4 temp_K; // Initialized as a 4x4 0-matrix
        SMatrix4 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_K /= denominator;
        temp_K += mat_bkg;
        // Loop over channels: 
        SVector4 B_factor;
        for (int j = 0; j < 4; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 2);
            GDouble B_den = barrierFactor(q_alpha, 2);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // Needed element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 4; j++) {
        mat_C(j, j) = KMatrix4_F2::chew(s, m1s[j], m2s[j]);
    }
    // Now mat_K should be done (ignore (s-s_0)/s_norm for now)
    SMatrix4 temp = SMatrixIdentity();
    mat_K *= mat_C;
    temp += mat_K;
    // temp *= ((s - 0.0091125) / 1); // Adler zero term only in f0
    // Now temp is the stuff that we want to invert before multiplying by the P-vector
    SMatrix4 temp_inv = KMatrix4_F2::inverse4(temp);
    // Now we cache the results
    // Note that it doesn't matter if I store it (i, j) or (j, i) as
    // long as I'm consistent later when unpacking it...
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            userVars[i * 4 + j + 2] = temp_inv(i, j).real(); // +2 because kM and kS are first in the enum
            userVars[i * 4 + j + 2 + 16] = temp_inv(i, j).imag(); // +16 to skip the real parts 
        }
    }
}



complex<GDouble> KMatrix4_F2::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector4 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(bf21270_re, bf21270_im),
        complex<GDouble>(bf21525_re, bf21525_im),
        complex<GDouble>(bf21810_re, bf21810_im),
        complex<GDouble>(bf21950_re, bf21950_im)
    };
    // Loop over resonances
    for (int i = 0; i < 4; i++) {
        SVector4 temp_P;
        SMatrix4 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_P /= denominator;
        // Loop over channels:
        SVector4 B_factor;
        for (int j = 0; j < 4; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 2);
            GDouble B_den = barrierFactor(q_alpha, 2);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        temp_P = temp_P * B_factor;
        vec_P += temp_P;
    }
    SMatrix4 temp_inv;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            temp_inv(i, j) = complex<GDouble>(userVars[i * 4 + j + 2],
                                              userVars[i * 4 + j + 2 + 16]);
        }
    }
    SVector4 res = temp_inv * vec_P;
    return res[2]; // return the KK channel contribution
}

SMatrix4 KMatrix4_F2::inverse4(SMatrix4 A) const {
    SMatrix4 M;
    SMatrix4 I = SMatrixIdentity();
    SMatrix4 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 4; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KMatrix4_F2::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KMatrix4_F2::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KMatrix4_F2::chew(double s, double m1, double m2) const {
    complex<GDouble> tot = 0;
    tot += (KMatrix4_F2::rho(s, m1, m2) / Pi()) * log((KMatrix4_F2::xi(s, m1, m2) + KMatrix4_F2::rho(s, m1, m2)) / (KMatrix4_F2::xi(s, m1, m2) - KMatrix4_F2::rho(s, m1, m2)));
    tot -= (KMatrix4_F2::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

void KMatrix4_F2::updatePar(const AmpParameter &par) {}
