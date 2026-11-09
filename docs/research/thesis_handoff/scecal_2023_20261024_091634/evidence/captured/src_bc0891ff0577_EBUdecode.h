#ifndef EbuDecode_H_
#define EbuDecode_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// scintillator strips wrt. 6 SPIROC2E chips * 36 channels
double* EBUdecode(int LayerIDs, int ChipIDs, int ChannelIDs);
double tempReconstruction(int LayerID, double* position, vector<double> tempSensor);

#endif
