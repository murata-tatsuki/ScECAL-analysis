#include <iostream>
//#include "TROOT.h"
#include <stdio.h>
#include <stdlib.h>
#include "RawtoRoot.h"
//#include "EventAlign.h"

using namespace std;

Int_t main(int argc, char* argv[])
{
	double start = clock();
	raw2Root tw;
	//tw.RMFelixTag(argv[1],argv[2]);
	//tw.BER(argv[2]);
	tw.Decode(argv[1], argv[2]);
	//tw.AnaSPS(argv[3],argv[4]);
	//tw.LEDCalib(argv[2],argv[2],8,1800,2500);
	//tw.AnaGainRatio(argv[1],argv[2]);
	//tw.AnaPed(argv[3],argv[2]);
	//tw.AnaCos(argv[3],argv[4]);
	//eventAlign tr(argv[1],argv[2],argv[4]);

	//return 0;
	double end = clock();
	cout << "end of RawToRoot : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
}
