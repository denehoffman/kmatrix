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
#include "KMatrix2_A0.h"


KMatrix2_A0::KMatrix2_A0(const vector<string> &args): UserAmplitude<KMatrix2_A0>(args) {
    /*
     * Usage: KMatrix2_A0 <daughter 1> <daughter 2> <Re[f_0(400)]> <Im[f_0(400)]> ...
     */
	m_daughters = pair<string, string>(args[0], args[1]);
	ba0980_re = AmpParameter(args[2]);
	ba0980_im = AmpParameter(args[3]);
    ba01450_re = AmpParameter(args[4]);
    ba01450_im = AmpParameter(args[5]);
    registerParameter(ba0980_re);
    registerParameter(ba0980_im);
    registerParameter(ba01450_re);
    registerParameter(ba01450_im);
    ga0980 = SVector2(0.43215, -0.28825);
    ga01450 = SVector2(0.19000, 0.43372);
    masses = {0.95395, 1.26767};
    couplings = {ga0980, ga01450};
    m1s = {0.1349768, 0.493677};
    m2s = {0.547862, 0.497611};
    a_bkg = {
        0.00000,
        0.00000, 0.00000};
    mat_bkg = SMatrix2Sym(a_bkg.begin(), a_bkg.end());
}

void KMatrix2_A0::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
    SMatrix2 mat_K; // Initialized as a 4x4 0-matrix
    SMatrix2 mat_C;
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SMatrix2 temp_K; // Initialized as a 4x4 0-matrix
        SMatrix2 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_K /= denominator;
        temp_K += mat_bkg;
        // Loop over channels: 
        SVector2 B_factor;
        for (int j = 0; j < 2; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 0);
            GDouble B_den = barrierFactor(q_alpha, 0);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        for (int k = 0; k < 2; k++) {
            for (int l = 0; l < 2; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // Needed element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 2; j++) {
        mat_C(j, j) = KMatrix2_A0::chew(s, m1s[j], m2s[j]);
    }
    // Now mat_K should be done (ignore (s-s_0)/s_norm for now)
    SMatrix2 temp = SMatrixIdentity();
    mat_K *= mat_C;
    temp += mat_K;
    // temp *= ((s - 0.0091125) / 1); // Adler zero term only in f0
    // Now temp is the stuff that we want to invert before multiplying by the P-vector
    SMatrix2 temp_inv = KMatrix2_A0::inverse2(temp);
    // Now we cache the results
    // Note that it doesn't matter if I store it (i, j) or (j, i) as
    // long as I'm consistent later when unpacking it...
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            userVars[i * 2 + j + 2] = temp_inv(i, j).real(); // +2 because kM and kS are first in the enum
            userVars[i * 2 + j + 2 + 4] = temp_inv(i, j).imag(); // +4 to skip the real parts 
        }
    }
}



complex<GDouble> KMatrix2_A0::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector2 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(ba0980_re, ba0980_im),
        complex<GDouble>(ba01450_re, ba01450_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector2 temp_P;
        SMatrix2 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_P /= denominator;
        // Loop over channels:
        SVector2 B_factor;
        for (int j = 0; j < 2; j++) {
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
    SMatrix2 temp_inv;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            temp_inv(i, j) = complex<GDouble>(userVars[i * 2 + j + 2],
                                              userVars[i * 2 + j + 2 + 4]);
        }
    }
    SVector2 res = temp_inv * vec_P;
    return res[1]; // return the KK channel contribution
}

SMatrix2 KMatrix2_A0::inverse2(SMatrix2 A) const {
    SMatrix2 M;
    SMatrix2 I = SMatrixIdentity();
    SMatrix2 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 2; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KMatrix2_A0::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KMatrix2_A0::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KMatrix2_A0::chew(double s, double m1, double m2) const {
    complex<GDouble> tot = 0;
    tot += (KMatrix2_A0::rho(s, m1, m2) / Pi()) * log((KMatrix2_A0::xi(s, m1, m2) + KMatrix2_A0::rho(s, m1, m2)) / (KMatrix2_A0::xi(s, m1, m2) - KMatrix2_A0::rho(s, m1, m2)));
    tot -= (KMatrix2_A0::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

void KMatrix2_A0::updatePar(const AmpParameter &par) {}
