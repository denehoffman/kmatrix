#if !(defined KSKSPLOTGENERATOR)
#define KSKSPLOTGENERATOR

#include <vector>
#include <string>

#include "IUAmpTools/PlotGenerator.h"

using namespace std;

class FitResults;
class Kinematics;

class KsKsPlotGenerator : public PlotGenerator {
    
public:
  
    // create an index for different histograms
    enum{kKsKsMass = 0, kCosTheta = 1, kPhi = 2, kt = 3, kCosThetaVMass = 4, kPhiVMass = 5, kNumHists = 6};
    KsKsPlotGenerator(const FitResults& results);
    
private:
        
    void projectEvent(Kinematics* kin);
  
};

#endif
