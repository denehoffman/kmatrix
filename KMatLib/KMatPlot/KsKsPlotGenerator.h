#if !(defined KSKSPLOTGENERATOR)
#define KSKSPLOTGENERATOR

#include <string>
#include <vector>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class KsKsPlotGenerator : public PlotGenerator {

public:
  // create an index for different histograms
  enum {
    kKsKsMass = 0,
    kKsKsMassFine,
    kKsKsMassCourse,
    kCosTheta,
    kPhi,
    kt,
    kCosThetaVMass,
    kPhiVMass,
    kNumHists
  };
  KsKsPlotGenerator(const FitResults &results);

private:
  void projectEvent(Kinematics *kin);
};

#endif
