

#include <cassert>
#include <iostream>
#include <string>
#include <complex>
#include <cstdlib>
#include <math.h>

#include "TLorentzVector.h"

#include "IUAmpTools/Kinematics.h"
#include "KMatAmp/PhaseSpace.h"
#include "breakupMomentum.h"

PhaseSpace::PhaseSpace(const vector<string>& args): UserAmplitude<PhaseSpace>(args) {
    assert(args.size() == 3);
    m_slope = AmpParameter(args[0]);
    registerParameter(m_slope);
    m_daughters = pair<string, string>(args[1], args[2]);
}

void PhaseSpace::calcUserVars(GDouble** pKin, GDouble* userVars) const {
    TLorentzVector pTemp, pTot;
    TLorentzVector p1, p2;
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
        p1 += pTemp;
        pTot += pTemp;
    }
    for(unsigned int i=0; i < m_daughters.second.size(); ++i) {
        string num; num += m_daughters.second[i];
        int index = atoi(num.c_str());
        pTemp.SetPxPyPzE(pKin[index][1],
                         pKin[index][2],
                         pKin[index][3],
                         pKin[index][0]);
        p2 += pTemp;
        pTot += pTemp;
    }
    // Grab our hypothesis mass
    GDouble m = pTot.M();
    userVars[kM] = m;
    // Grab phase space
    complex<GDouble> rho = PhaseSpace::rho(pTot.M2(), p1.M(), p2.M());
    userVars[kRhoReal] = rho.real();
    userVars[kRhoImag] = rho.imag();
}

complex<GDouble> PhaseSpace::calcAmplitude(GDouble** pKin, GDouble* userVars) const {
    GDouble m = userVars[kM];
    complex<GDouble> rho = complex<GDouble>(userVars[kRhoReal], userVars[kRhoImag]);
    GDouble falloff = TMath::Exp(- m_slope * m * 0.0001);
    //cout << "M: " << m << " lda: " << m_slope << " rho: " << rho << " falloff: " << falloff << endl;
    //cout << "Result: " << (rho * falloff) << endl;
    return rho * falloff;
}

complex<GDouble> PhaseSpace::rho(double s, double m1, double m2) const {
    return sqrt(complex<GDouble>(((1 - ((m1 + m2) * (m1 + m2) / s)) * (1 - ((m1 - m2) * (m1 - m2) / s))), 0));
}

void PhaseSpace::updatePar(const AmpParameter &par) {}
