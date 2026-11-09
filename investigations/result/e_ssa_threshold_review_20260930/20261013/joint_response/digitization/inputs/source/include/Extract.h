#ifndef _EXTRACT_H
#define _EXTRACT_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
// #include "MCDigi.h"

using namespace std;

class Extract
{
	public:
		static constexpr const double MIPEnergy = 0.305;
		static const int LAYERNO = 30;
		static const int CHIPNO = 6;
		static const int CHANNO = 36;
		static const int CELLNO = 6300;
		void PedestalExtract(string inPedFileName);
		void HLGainRatioExtract(string inHLGainRatioFileName);
		void MIPPeakValueExtract(string inMIPFileName);
		void SiPMGainExtract(string SiPMGainFileName);
		void TresholdExtract(string inThresholdFileName, int thre);
		void TresholdExtract(int thre);
	
		double pedMeanHigh[LAYERNO][CHIPNO][CHANNO]; 
		double pedSigmaHigh[LAYERNO][CHIPNO][CHANNO]; 
		double pedMeanLow[LAYERNO][CHIPNO][CHANNO]; 
		double pedSigmaLow[LAYERNO][CHIPNO][CHANNO]; 

		double HLRatio[LAYERNO][CHIPNO][CHANNO];
		double HLIntercept[LAYERNO][CHIPNO][CHANNO]; // LG ADC intercept = -calibration HG offset / HLRatio
		double HighADCSwitch[LAYERNO][CHIPNO][CHANNO];
		
		double MIPPeakValue[LAYERNO][CHIPNO][CHANNO]; 
		double SiPMGain[LAYERNO][CHIPNO][CHANNO]; 
		double SiPMGainError[LAYERNO][CHIPNO][CHANNO]; 

		double thresholdValue[LAYERNO][CHIPNO][CHANNO];
		double thresholdSigmaValue[LAYERNO][CHIPNO][CHANNO];
		double thresholdValue_with_pedestal[LAYERNO][CHIPNO][CHANNO];
		double thresholdValue_chip[LAYERNO][CHIPNO];
		double thresholdSigmaValue_chip[LAYERNO][CHIPNO];
		double thresholdValue_chipCount[LAYERNO][CHIPNO];
};

#endif
