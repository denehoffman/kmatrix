#include <iostream>
#include <string>

#include "TClass.h"
#include "TApplication.h"
#include "TGClient.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TClass.h"

#include "IUAmpToolsMPI/AmpToolsInterfaceMPI.h"
#include "IUAmpTools/FitResults.h"

#include "AmpPlotter/PlotterMainWindow.h"
#include "AmpPlotter/PlotFactory.h"

#include "KMatPlot/KsKsPlotGenerator.h"
#include "KMatDataIO/ROOTDataReader.h"
#include "KMatAmp/Uniform.h"
#include "KMatAmp/Zlm.h"
#include "KMatAmp/BreitWigner.h"
#include "KMatAmp/KMatrix2_A0.h"
#include "KMatAmp/KMatrix3_A2.h"
#include "KMatAmp/KMatrix4_F2.h"
#include "KMatAmp/KMatrix5_F0.h"

typedef KsKsPlotGenerator PlotGen;

void atiSetup() {
  
    // the PlotGenerator will create an AmpToolsInterface in order
    // to create plots - this setup must happen before the
    // AmpToolsInterface is created

    AmpToolsInterface::registerAmplitude(Uniform());
    AmpToolsInterface::registerAmplitude(Zlm());
    AmpToolsInterface::registerAmplitude(BreitWigner());
    AmpToolsInterface::registerAmplitude(KMatrix2_A0());
    AmpToolsInterface::registerAmplitude(KMatrix3_A2());
    AmpToolsInterface::registerAmplitude(KMatrix4_F2());
    AmpToolsInterface::registerAmplitude(KMatrix5_F0());

    AmpToolsInterface::registerDataReader(ROOTDataReader());
}


//  THE USER SHOULD NOT HAVE TO CHANGE ANYTHING BELOW THIS LINE
// *************************************************************

using namespace std;

int main( int argc, char* argv[] ) {

        // ************************
        // usage
        // ************************

    cout << endl << " *** Viewing Results Using AmpPlotter *** " << endl << endl;

    if (argc <= 2) {
        
        cout << "Usage:" << endl << endl;
        cout << "\tampPlotter <fit results name>" << endl << endl;
        return 0;
    }

        // ************************
        // parse the command line parameters
        // ************************
    string resultsName(argv[1]);
    string outName(argv[2]);
    cout << "Loading fit results" << endl;
    FitResults results(resultsName);
    if(!results.valid()) {
        cout << "Invalid fit results in file:  " << resultsName << endl;
        exit(1);
    }

        // ************************
        // set up the plot generator
        // ************************

    atiSetup();  
    PlotGen plotGen(results);
    TFile* plotfile = new TFile(outName.c_str(), "recreate");
    TH1::AddDirectory(kFALSE);
    // enable all reactions (I think this means all the polarizations we loop over)
    vector<string> reactionList = results.reactionList();
    for (unsigned int ireact = 0; ireact < reactionList.size(); ireact++) {
        plotGen.enableReaction(reactionList[ireact]);
    }
    vector<string> amps = plotGen.uniqueAmplitudes();
    // loop over amplitude configurations
    for (unsigned int iamp = 0; iamp <= amps.size(); iamp++) {
        // ^ on the final iteration, iamp is out of bounds so
        // this just plots the sum of all amplitudes
        // turn on all amplitudes
        for (unsigned int i = 0; i < amps.size(); i++) {
            plotGen.enableAmp(i);
        }
        // turn off all amplitudes but one (or keep all on)
        if (iamp < amps.size()){
            for (unsigned int i = 0; i < amps.size(); i++){
                if (i != iamp) plotGen.disableAmp(i);
            }      
        }

        // loop over data, accMC, and genMC
        for (unsigned int iplot = 0; iplot < PlotGen::kNumTypes; iplot++) {
            if (iamp < amps.size() && iplot == PlotGen::kData) continue; // skip data plots for the individual amps

            // loop over different variables
            for (unsigned int ivar  = 0; ivar  < KsKsPlotGenerator::kNumHists; ivar++) {
                for (unsigned int ireact = 0; ireact < reactionList.size(); ireact++) {
                    // name the histogram according to all the information
                    string histname =  "h"; // sure, why not...
                    if (ivar == KsKsPlotGenerator::kKsKsMass)  histname += "KsKsMass";
                    if (ivar == KsKsPlotGenerator::kCosTheta)  break;//histname += "CosTheta";
                    if (ivar == KsKsPlotGenerator::kPhi)  break;//histname += "Phi";
                    if (ivar == KsKsPlotGenerator::kt)  break;//histname += "t";
                    if (ivar == KsKsPlotGenerator::kCosThetaVMass) histname += "CosThetaVMass";
                    if (ivar == KsKsPlotGenerator::kPhiVMass) histname += "PhiVMass";
                    if (iplot == PlotGen::kData) histname += "dat"; // only true when iamp >= amps.size() aka the total of all amplitudes
                    else if (iplot == PlotGen::kAccMC) histname += "acc";
                    else break;
                    //if (iplot == PlotGen::kGenMC) histname += "gen";
                    if (iplot == PlotGen::kGenMC) break; // maybe it won't load the generated MC?
                    if (iamp < amps.size()){
                        string ampName = amps[iamp];
                        histname += "_";
                        histname += ampName; // add the amplitude name to the fit also!
                    }

                    Histogram* hist = plotGen.projection(ivar, reactionList[ireact], iplot);
                    // set the axis labels according to the histogram information
                    string xtitle = "Mass(K_{S}K_{S}) (GeV/c^{2})";
                    string ytitle = "Events";
                    if (ivar == KsKsPlotGenerator::kCosTheta) xtitle = "cos(#theta_{HX})";
                    if (ivar == KsKsPlotGenerator::kPhi) xtitle = "#phi_{HX}";
                    if (ivar == KsKsPlotGenerator::kt) xtitle = "-t";
                    if (ivar == KsKsPlotGenerator::kCosThetaVMass) ytitle = "cos(#theta_{HX})";
                    if (ivar == KsKsPlotGenerator::kPhiVMass) ytitle = "#phi_{HX}";

                    TH1* thist = hist->toRoot(); // TH1 is the parent class for 2D histograms also
                    thist->SetName(histname.c_str());
                    thist->SetStats(0);
                    //thist->SetTitle();
                    thist->SetXTitle(xtitle.c_str());
                    thist->SetYTitle(ytitle.c_str());
                    plotfile->cd();
                    thist->Write();
                } // end ireact loop
            } // end ivar loop
        } // end iplot loop
    } // end iamp loop

  plotfile->Close();



        // ************************
        // start the GUI
        // ************************

    //cout << ">> Plot generator ready, starting GUI..." << endl;

    /*int dummy_argc = 0;
    char* dummy_argv[] = {};  
    TApplication app("app", &dummy_argc, dummy_argv);
    
    gStyle->SetFillColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetPadColor(10);
    gStyle->SetFillStyle(1001);
    gStyle->SetPalette(1);
    gStyle->SetFrameFillColor(10);
    gStyle->SetFrameFillStyle(1001);
    
    PlotFactory factory(plotGen);	
    PlotterMainWindow mainFrame(gClient->GetRoot(), factory);
        
    app.Run();
    */
    return 0;
}
