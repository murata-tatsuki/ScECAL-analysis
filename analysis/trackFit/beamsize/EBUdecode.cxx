#include "EBUdecode.h"
#include <fstream>
#include <iostream>

using namespace std;

// scintillator strips wrt. 6 SPIROC2E chips * 36 channels
double* EBUdecode(int LayerIDs, int ChipIDs, int ChannelIDs)
{
	// static const int layerNu = 32;
	static const int chipNu = 6;
	static const int chnNu = 36;
	int decodeID[chipNu][chnNu] = { 0, 42, 1, 43, 2, 44, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 54, 13, 55, 14, 56, 15, 57, 16, 58, 17, 59, 18, 60, 19, 61, 20, 62, 21, 22, 23,
																	24, 66, 25, 67, 26, 68, 27, 69, 28, 70, 29, 71, 30, 72, 31, 73, 32, 74, 33, 75, 34, 76, 35, 77, 36, 78, 37, 79, 38, 80, 39, 81, 40, 82, 41, 83,
																	149, 148, 147, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 63, 64, 65, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
																	150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197, 156, 198, 157, 199, 158, 200, 159, 201, 160, 202, 161, 203, 162, 204, 163, 205, 164, 206, 165, 207, 166, 208, 167, 209,
																	191, 190, 189, 188, 146, 187, 145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181, 139, 180, 138, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
																	137, 136, 135, 134, 133, 132, 131, 130, 129, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 45, 46, 47, 48, 49, 50, 51, 52, 53, 210, 210, 210, 210, 210, 210 };

	int ScintillatorIDs = decodeID[ChipIDs][ChannelIDs];
	double layerZ;

	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	int _yID = ScintillatorIDs / rowNu;
	int _xID = ScintillatorIDs % rowNu;
	static double _position[3] = { 0 };
	double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;
	double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;

	// for prototype test
	if (LayerIDs % 2 == 0) {
		_position[0] = -y0;
		_position[1] = -x0;
	}
	if (LayerIDs % 2 == 1) {
		_position[0] = -x0;
		_position[1] = -y0;
	}
	if (LayerIDs % 2 == 0)
		layerZ = 1 + LayerIDs / 2 * 19.9;
	//else layerZ = 12.95+(LayerIDs-1)/2*19.9;
	else
		layerZ = 12.2 + (LayerIDs - 1) / 2 * 19.9;
	_position[2] = layerZ;
	// for grouped CR test
	// if(LayerIDs%4==0) {
	//     _position[0] = -y0;
	//     _position[1] = -x0;
	// }
	// if(LayerIDs%4==1) {
	//     _position[0] = x0;
	//     _position[1] = y0;
	// }
	// if(LayerIDs%4==2) {
	//     _position[0] = y0;
	//     _position[1] = x0;
	// }
	// if(LayerIDs%4==3) {
	//     _position[0] = -x0;
	//     _position[1] = -y0;
	// }
	// _position[2] = layerZ[LayerIDs];

	return _position;
}

double tempReconstruction(int LayerID, double* position, vector<double> tempSensor)
{
	double localTemperature;
	const int sensorNo = 16;
	double posX[sensorNo] = { 39.5986, -106.197, 89.9922, -22.1234, -103.886, -0.508, -60.2488, 59.1312, 28.2702, -50.3428, 100.99, -106.68, -106.299, 105.918, -47.0916, 53.2638 };
	double posY[sensorNo] = { -34.3916, -34.0868, -11.1506, -13.1064, 8.0772, 7.366, 34.5186, 33.7312, 54.0766, 75.819, 101.473, 101.498, -101.752, -102.032, -70.6882, -70.3326 };
	double sum = 0;
	double length = 0;
	double x, y;
	for (int g = 0; g < sensorNo; g++) {
		if (LayerID % 2 == 0) {
			x = -posY[g];
			y = -posX[g];
		}
		if (LayerID % 2 == 1) {
			x = -posX[g];
			y = -posY[g];
		}
		length = std::sqrt(std::pow((*position - x), 2) + std::pow((*(position + 1) - y), 2));
		sum += 1 / length;
	}
	for (int g = 0; g < sensorNo; g++) {
		if (LayerID % 2 == 0) {
			x = -posY[g];
			y = -posX[g];
		}
		if (LayerID % 2 == 1) {
			x = -posX[g];
			y = -posY[g];
		}
		//if(LayerID%2==0) {x = posY[g]; y = -posX[g];}
		//if(LayerID%2==1) {x = -posX[g]; y = posY[g];}
		length = std::sqrt(std::pow((*position - x), 2) + std::pow((*(position + 1) - y), 2));
		double weight = 1 / length / sum;
		localTemperature += weight * tempSensor[g];
	}
	return localTemperature;
}
