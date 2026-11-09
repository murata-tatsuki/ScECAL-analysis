#include <iostream>
#include <stdexcept>
#include "MCDigi.h"
#include "Extract.h"
#include "TemperaturePhysical.hh"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 8 && argc != 9) {
        cerr << "Usage: MCAnalysis input.root pedestal.root hl.root mip.root sipm.root threshold.root output.root [threshold]" << endl;
        return 1;
    }
    try {
        MCDigi MCD;
        Extract EXTR;
        int thre = argc == 9 ? stoi(argv[8]) : 0;
        cout << "threshold file " << argv[6] << endl;
        cout << "threshold " << thre << endl;
        cout << "Process : Simulation Digi" << endl;
        EXTR.PedestalExtract(argv[2]);
        EXTR.HLGainRatioExtract(argv[3]);
        EXTR.MIPPeakValueExtract(argv[4]);
        EXTR.SiPMGainExtract(argv[5]);
 TemperaturePhysical::initialize(EXTR,argv[5]);
        EXTR.TresholdExtract(argv[6], thre);
        MCD.SimuDigi(argv[1], argv[7], EXTR);
        cout << "Finish  : Simulation Digi" << endl;
    } catch (const std::exception& error) {
        cerr << "Digitization failed: " << error.what() << endl;
        return 1;
    }
    return 0;
}
