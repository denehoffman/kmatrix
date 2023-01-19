#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "KsKsPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

KsKsPlotGenerator::KsKsPlotGenerator(const FitResults& results) : PlotGenerator(results) {
    bookHistogram(kKsKsMass, new Histogram1D(100, 1.0, 2.0, "ksksMass", "Invariant Mass of K_{S}K_{S}"));
    bookHistogram(kKsKsMassFine, new Histogram1D(200, 1.0, 2.0, "ksksMassFine", "Invariant Mass of K_{S}K_{S}"));
    bookHistogram(kCosTheta, new Histogram1D(50, -1., 1., "cosTheta", "cos(#theta_{HX}) of Resonance Production"));
    bookHistogram(kPhi, new Histogram1D(50, -180, 180, "Phi", "#phi_{HX}"));
    bookHistogram(kt, new Histogram1D(100, 0, 2.00, "t", "-t"));
    bookHistogram(kCosThetaVMass, new Histogram2D(100, 1.0, 2.0, 50, -1., 1., "cosThetaVMass", "cos(#theta_{HX}) vs. m(K_{S}K_{S})"));
    bookHistogram(kPhiVMass, new Histogram2D(100, 1.0, 2.0, 50, -180, 180, "phiVMass", "#phi_{HX} vs. m(K_{S}K_{S})"));
}

void KsKsPlotGenerator::projectEvent(Kinematics* kin) {
    
    TLorentzVector beam_lab   = kin->particle(0); // Beam
    TLorentzVector recoil_lab = kin->particle(1); // Proton
    TLorentzVector p1_lab = kin->particle(2);     // K-Short 1
    TLorentzVector p2_lab = kin->particle(3);     // K-Short 2
    
    TLorentzVector resonance_lab = p1_lab + p2_lab;

    TLorentzVector com = -1. * (recoil_lab + p1_lab + p2_lab);
    TLorentzRotation comBoost(-com.BoostVector());

    TLorentzVector beam = comBoost * beam_lab;
    TLorentzVector recoil = comBoost * recoil_lab;
    TLorentzVector p1 = comBoost * p1_lab;
    TLorentzVector resonance = comBoost * resonance_lab;

    TLorentzRotation resRestBoost(-resonance.BoostVector());

    TLorentzVector beam_res   = resRestBoost * beam;
    TLorentzVector recoil_res = resRestBoost * recoil;
    TLorentzVector p1_res = resRestBoost * p1;

    // Helicity frame
    TVector3 z = -1. * recoil_res.Vect().Unit();
    // or GJ frame?
    // TVector3 z = beam_res.Vect().Unit();

    // normal to the production plane
    TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

    TVector3 x = y.Cross(z);

    TVector3 angles((p1_res.Vect()).Dot(x), (p1_res.Vect()).Dot(y), (p1_res.Vect()).Dot(z));
    
    // compute angles
    double cosTheta = angles.CosTheta();
    double phi = angles.Phi()*180./TMath::Pi();

    // compute invariant t
    TLorentzVector TargetP4;
    TargetP4.SetPxPyPzE(0,0,0,0.938272);
    GDouble t=(recoil-TargetP4).Mag2();
    
    fillHistogram(kKsKsMass, resonance.M());
    fillHistogram(kKsKsMassFine, resonance.M());
    fillHistogram(kCosTheta, cosTheta);
    fillHistogram(kPhi, phi);
    fillHistogram(kt, -t); // fill with -t to make positive
    fillHistogram(kCosThetaVMass, resonance.M(), cosTheta);
    fillHistogram(kPhiVMass, resonance.M(), phi);
}
