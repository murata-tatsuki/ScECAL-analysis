#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include "Calibration.h"
#include "TROOT.h" //https://root-forum.cern.ch/t/trees-with-vectors/8276/2

using namespace std;

int main(int argc, char* argv[])
{
	double start = clock();
	gROOT->ProcessLine("#include <vector>"); //https://root-forum.cern.ch/t/trees-with-vectors/8276/2
	cout << "start of main program" << endl;

// 	treeCalibration t;
// 	t.Calibration(argv[1], argv[2], argv[3]);

    if (argc != 6) {
        cerr << "Usage: calib input.root pedestal.root hl.root mip.root output.root" << endl;
        return 1;
    }
    try {
        processCali cali;
        cali.Calibration(argv[1],argv[2],argv[3],argv[4],argv[5]);
    } catch (const std::exception& error) {
        cerr << "Calibration failed: " << error.what() << endl;
        return 1;
    }

	cout << "end of main program" << endl;
	double end = clock();
	cout << "end of RawToRoot : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
	return 0;
}
