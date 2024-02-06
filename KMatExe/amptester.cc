#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"

#include "IUAmpTools/AmpToolsInterface.h"
#include "IUAmpTools/AmplitudeManager.h"
#include "IUAmpTools/ConfigFileParser.h"
#include "IUAmpTools/ConfigurationInfo.h"

#include "KMatAmp/BreitWigner.h"
#include "KMatAmp/KMatrix2_A0_KsKs.h"
#include "KMatAmp/KMatrix3_A2_KsKs.h"
#include "KMatAmp/KMatrix4_F2_KsKs.h"
#include "KMatAmp/KMatrix5_F0_KsKs.h"
#include "KMatAmp/KMatrix2_A0.h"
#include "KMatAmp/KMatrix3_A2.h"
#include "KMatAmp/KMatrix4_F2.h"
#include "KMatAmp/KMatrix5_F0.h"
#include "KMatAmp/KMatrix2_Pi1.h"
#include "KMatAmp/KMatrix5_F0b.h"
#include "KMatAmp/PhaseSpace.h"
#include "KMatAmp/Uniform.h"
#include "KMatAmp/Ylm.h"
#include "KMatAmp/Zlm.h"
#include "KMatDataIO/ROOTDataReader.h"
#include "KMatDataIO/ROOTDataReaderBootstrap.h"

using namespace std;

#include "IUAmpTools/report.h"
static const char *kModule = "printAmplitudes";

int main(int argc, char **argv) {

  // ************************
  // usage
  // ************************

  if (argc <= 1) {
    report(NOTICE, kModule) << "Usage:" << endl << endl;
    report(NOTICE, kModule) << "\tprintAmplitudes <config file name>" << endl
                            << endl;
    return 0;
  }

  report(INFO, kModule) << endl
                        << " *** Printing Amplitudes *** " << endl
                        << endl;

  // ************************
  // parse the command line parameters
  // ************************

  string cfgname(argv[1]);

  report(INFO, kModule) << "Config file name = " << cfgname << endl << endl;

  // ************************
  // parse the config file
  // ************************

  ConfigFileParser parser(cfgname);
  ConfigurationInfo *cfgInfo = parser.getConfigurationInfo();
  cfgInfo->display();

  // ************************
  // AmpToolsInterface
  // ************************

  AmpToolsInterface::registerAmplitude(BreitWigner());
  AmpToolsInterface::registerAmplitude(Zlm());
  AmpToolsInterface::registerAmplitude(Ylm());
  AmpToolsInterface::registerAmplitude(KMatrix2_A0_KsKs());
  AmpToolsInterface::registerAmplitude(KMatrix3_A2_KsKs());
  AmpToolsInterface::registerAmplitude(KMatrix4_F2_KsKs());
  AmpToolsInterface::registerAmplitude(KMatrix5_F0_KsKs());
  AmpToolsInterface::registerAmplitude(KMatrix2_A0());
  AmpToolsInterface::registerAmplitude(KMatrix3_A2());
  AmpToolsInterface::registerAmplitude(KMatrix4_F2());
  AmpToolsInterface::registerAmplitude(KMatrix5_F0());
  AmpToolsInterface::registerAmplitude(KMatrix2_Pi1());
  AmpToolsInterface::registerAmplitude(KMatrix5_F0b());
  AmpToolsInterface::registerAmplitude(Uniform());
  AmpToolsInterface::registerAmplitude(PhaseSpace());
  AmpToolsInterface::registerDataReader(ROOTDataReader());
  AmpToolsInterface::registerDataReader(ROOTDataReaderBootstrap());

  AmpToolsInterface ATI(cfgInfo, AmpToolsInterface::kPlotGeneration);

  DataReader *dataReader =
      ATI.genMCReader(cfgInfo->reactionList()[0]->reactionName());
  for (int i = 0; i < 10; i++) {
    Kinematics *kin = dataReader->getEvent();
    ATI.printEventDetails(cfgInfo->reactionList()[0]->reactionName(), kin);
    delete kin;
  }
}
