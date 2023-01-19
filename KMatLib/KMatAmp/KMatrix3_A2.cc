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
#include "KMatrix3_A2.h"


KMatrix3_A2::KMatrix3_A2(const vector<string> &args): UserAmplitude<KMatrix3_A2>(args) {
    /*
     * Usage: KMatrix3_A2 <daughter 1> <daughter 2> <Re[f_0(400)]> <Im[f_0(400)]> ...
     */
	m_daughters = pair<string, string>(args[0], args[1]);
	ba21320_re = AmpParameter(args[2]);
	ba21320_im = AmpParameter(args[3]);
    ba21700_re = AmpParameter(args[4]);
    ba21700_im = AmpParameter(args[5]);
    registerParameter(ba21320_re);
    registerParameter(ba21320_im);
    registerParameter(ba21700_re);
    registerParameter(ba21700_im);
    ga21320 = SVector3(0.30073, 0.21426, -0.09162);
    ga21700 = SVector3(0.68567, 0.12543, 0.00184);
    masses = {1.30080, 1.75351};
    couplings = {ga21320, ga21700};
    m1s = {0.1349768, 0.493677, 0.1349768};
    m2s = {0.547862, 0.497611, 0.95778};
    a_bkg = {
        -0.40184,
        0.00033, -0.21416,
        -0.08707, -0.06193, -0.17435};
    mat_bkg = SMatrix3Sym(a_bkg.begin(), a_bkg.end());
}

void KMatrix3_A2::calcUserVars(GDouble** pKin, GDouble* userVars) const {
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
    SMatrix3 mat_K; // Initialized as a 4x4 0-matrix
    SMatrix3 mat_C;
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SMatrix3 temp_K; // Initialized as a 4x4 0-matrix
        SMatrix3 temp_B;
        temp_K = TensorProd(couplings[i], couplings[i]);
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_K /= denominator;
        temp_K += mat_bkg;
        // Loop over channels: 
        SVector3 B_factor;
        for (int j = 0; j < 3; j++) {
            GDouble q_alpha = fabs(breakupMomentum(masses[i], m1s[j], m2s[j]));
            GDouble q = fabs(breakupMomentum(m, m1s[j], m2s[j]));
            GDouble B_num = barrierFactor(q, 2);
            GDouble B_den = barrierFactor(q_alpha, 2);
            B_factor[j] = B_num / B_den;
        }
        temp_B = TensorProd(B_factor, B_factor); 
        for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                temp_K(k, l) = temp_K(k, l) * temp_B(k, l); // Needed element-wise multiplication
            }
        }
        mat_K += temp_K;
    }
    // Loop over channels
    for (int j = 0; j < 3; j++) {
        mat_C(j, j) = KMatrix3_A2::chew(s, m1s[j], m2s[j]);
    }
    // Now mat_K should be done (ignore (s-s_0)/s_norm for now)
    SMatrix3 temp = SMatrixIdentity();
    mat_K *= mat_C;
    temp += mat_K;
    // temp *= ((s - 0.0091125) / 1); // Adler zero term only in f0
    // Now temp is the stuff that we want to invert before multiplying by the P-vector
    SMatrix3 temp_inv = KMatrix3_A2::inverse3(temp);
    // Now we cache the results
    // Note that it doesn't matter if I store it (i, j) or (j, i) as
    // long as I'm consistent later when unpacking it...
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            userVars[i * 3 + j + 2] = temp_inv(i, j).real(); // +2 because kM and kS are first in the enum
            userVars[i * 3 + j + 2 + 9] = temp_inv(i, j).imag(); // +9 to skip the real parts 
        }
    }
}



complex<GDouble> KMatrix3_A2::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM]; 
    GDouble s = userVars[kS];
    SVector3 vec_P;
    vector<complex<GDouble>> betas{
        complex<GDouble>(ba21320_re, ba21320_im),
        complex<GDouble>(ba21700_re, ba21700_im),
    };
    // Loop over resonances
    for (int i = 0; i < 2; i++) {
        SVector3 temp_P;
        SMatrix3 temp_B;
        temp_P = couplings[i];
        temp_P *= betas[i]; 
        GDouble denominator = (masses[i] * masses[i] - s);
        // if (masses[i] == m) {denominator += 1E-3;} // just in case
        if (denominator < 1E-6) {denominator += 1E-3;} // just in case
        temp_P /= denominator;
        // Loop over channels:
        SVector3 B_factor;
        for (int j = 0; j < 3; j++) {
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
    SMatrix3 temp_inv;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            temp_inv(i, j) = complex<GDouble>(userVars[i * 4 + j + 2],
                                              userVars[i * 4 + j + 2 + 9]);
        }
    }
    SVector3 res = temp_inv * vec_P;
    return res[1]; // return the KK channel contribution
}

SMatrix3 KMatrix3_A2::inverse3(SMatrix3 A) const {
    SMatrix3 M;
    SMatrix3 I = SMatrixIdentity();
    SMatrix3 temp;
    complex<GDouble> c = 1;
    for (int k = 1; k <= 3; k++) {
        M = A * M + I * c;
        temp = A * M;
        c = temp.Trace() / (-1.0 * k);
    }
    return M / (-1.0 * c);
}

complex<GDouble> KMatrix3_A2::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

complex<GDouble> KMatrix3_A2::xi(double s, double m1, double m2) const {
    return complex<GDouble>(1 - ((m1 + m2) * (m1 + m2) / s), 0);
}

complex<GDouble> KMatrix3_A2::chew(double s, double m1, double m2) const {
    complex<GDouble> tot = 0;
    tot += (KMatrix3_A2::rho(s, m1, m2) / Pi()) * log((KMatrix3_A2::xi(s, m1, m2) + KMatrix3_A2::rho(s, m1, m2)) / (KMatrix3_A2::xi(s, m1, m2) - KMatrix3_A2::rho(s, m1, m2)));
    tot -= (KMatrix3_A2::xi(s, m1, m2) / Pi()) * ((m2 - m1) / (m1 + m2)) * log(m2 / m1);
    return tot;
}

void KMatrix3_A2::updatePar(const AmpParameter &par) {}
