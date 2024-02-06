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
#include "KMatrix2_Pi1.h"


KMatrix2_Pi1::KMatrix2_Pi1(const vector<string> &args): UserAmplitude<KMatrix2_Pi1>(args) {
    /*
     * Usage: KMatrix2_Pi1 <daughter 1> <daughter 2> <channel> <Re[pi_1(1600)]> <Im[pi_1(1600)]>
     */
    m_daughters = pair<string, string>(args[0], args[1]);
    channel = atoi(args[2].c_str());
    bpi11600_re = AmpParameter(args[3]);
    bpi11600_im = AmpParameter(args[4]);
    registerParameter(bpi11600_re);
    registerParameter(bpi11600_im);
    gpi11600 = SVector2(0.80564, 1.04695);
    masses = {1.38552};
    couplings = {gpi11600};
    m1s = {0.1349768, 0.1349768};
    m2s = {0.547862, 0.95778};
    a_bkg = {
        1.05000,
        0.15163, -0.24611};
    mat_bkg = SMatrix2Sym(a_bkg.begin(), a_bkg.end());
}

void KMatrix2_Pi1::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
        // if (denominator < 1E-6) {denominator += 1E-3;} // just in case
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
        mat_C(j, j) = KMatrix2_Pi1::chew(s, m1s[j], m2s[j]);
    }
    // Now mat_K should be done (ignore (s-s_0)/s_norm for now)
    SMatrix2 temp = SMatrixIdentity();
    mat_K *= mat_C;
    temp += mat_K;
    // temp *= ((s - 0.0091125) / 1); // Adler zero term only in f0
    // Now temp is the stuff that we want to invert before multiplying by the P-vector
    SMatrix2 temp_inv = KMatrix2_Pi1::inverse2(temp);
    // Now we cache the results
    for (int i = 0; i < 2; i++) {
        userVars[i + 2] = temp_inv(channel, i).real(); // +2 because kM and kS are first in the enum
        userVars[i + 2 + 2] = temp_inv(channel, i).imag(); // +2 to skip the real parts 
    }
}



complex<GDouble> KMatrix2_Pi1::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector2 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(bpi11600_re, bpi11600_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector2 temp_P;
        SMatrix2 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        // if (denominator < 1E-6) {denominator += 1E-3;} // just in case
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
    SVector2 temp_inv;
    for (int i = 0; i < 2; i++) {
        temp_inv(i) = complex<GDouble>(userVars[i + 2], userVars[i + 2 + 2]);
    }
    return Dot(temp_inv, vec_P);
}

SMatrix2 KMatrix2_Pi1::inverse2(SMatrix2 A) const {
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

complex<GDouble> KMatrix2_Pi1::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KMatrix2_Pi1::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KMatrix2_Pi1::chew(double s, double m1, double m2) const {
    complex<GDouble> tot = 0;
    tot += (KMatrix2_Pi1::rho(s, m1, m2) / Pi()) * log((KMatrix2_Pi1::xi(s, m1, m2) + KMatrix2_Pi1::rho(s, m1, m2)) / (KMatrix2_Pi1::xi(s, m1, m2) - KMatrix2_Pi1::rho(s, m1, m2)));
    tot -= (KMatrix2_Pi1::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

void KMatrix2_Pi1::updatePar(const AmpParameter &par) {}
