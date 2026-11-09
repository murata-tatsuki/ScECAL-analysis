#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2D.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <cmath>
#include "TLegend.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "EBUdecode.h"
#include "EBUdecode.cxx"
#include "langaus.C"        // langaus(chargeH, &fitFunc, &peakP, &peakPError);
// #include "../RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest_SPS2022/analysis/ECAL_Analysis/include/EBUdecode.h"

using namespace std;


// エネルギー分布を求めるマクロ
// calibを使う

const bool pedestal_figure_saving = false;


double* EBUdecode_itr(int LayerIDs, int ChipIDs, int ChannelIDs){
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

// 解析結果をまとめて保持する構造体
struct ThresholdFitResult {
    bool isSuccess;        // フィット成功フラグ
    double Ath;            // 閾値 [ADC]
    double AthErr;         // 閾値のエラー [ADC]
    double sigma;          // なまり幅 (ノイズ幅) [ADC]
    double sigmaErr;       // なまり幅のエラー [ADC]
    double p0;             // 指数関数の振幅
    double p1;             // 指数関数の減衰率
    double p2;             // オフセット (定数ノイズ)
    double chi2ndf;        // Chi2 / NDF (フィット精度の指標)
    TF1* f_draw = nullptr; // 描画用TF1のポインタ
};

// 1. フィット関数の定義 ( Exp * erf )
Double_t ThresholdExpFitFunc(Double_t *x, Double_t *par) {
    Double_t A     = x[0];
    Double_t p0    = par[0]; // 指数振幅
    Double_t p1    = par[1]; // 減衰率
    Double_t p2    = par[2]; // オフセット
    Double_t Ath   = par[3]; // 閾値 (ADC)
    Double_t sigma = par[4]; // なまり幅 (ADC)

    if (sigma <= 0.0) return 0.0;

    // バックグラウンド関数（本来の入力スペクトル）
    Double_t bg  = p0 * TMath::Exp(-p1 * A) + p2;
    // エレクトロニクスの検出効率 (Sカーブ)
    Double_t eff = 0.5 * (1.0 + TMath::Erf((A - Ath) / (TMath::Sqrt(2.0) * sigma)));

    return bg * eff;
}

// 2. メインの解析関数
ThresholdFitResult FitChannelThreshold(TH1F* h_adc_in, int rebinFactor = 4, double fitRangeWidthLow = 15.0, double fitRangeWidthHigh = 40.0, double fitMin = -1., double fitMax = -1.) {
  ThresholdFitResult res;
  res.isSuccess = false;

  // 入力チェック（データが空の場合）
  if (!h_adc_in || h_adc_in->GetEntries() < 100) {
    std::cerr << "[Warning] Histogram is empty or has too few entries." << std::endl;
    return res;
  }

  // 元のヒストグラムを変更しないようクローンを作成しRebin
  TH1F* h_fit = (TH1F*)h_adc_in->Clone(Form("%s_rebinned", h_adc_in->GetName()));
  if (rebinFactor > 1) h_fit->Rebin(rebinFactor);

  // --- 初期値の自動推定 (Heuristics) ---
  double maxContent = h_fit->GetMaximum();
  int maxBin = h_fit->GetMaximumBin();
  double maxX = h_fit->GetBinCenter(maxBin);

  // 最大カウントの 15% を超える最初のビンを「閾値の初期値」とする
  double Ath_init = maxX;
  if (fitMin == -1 || fitMax == -1) {
    for (int i = 1; i <= h_fit->GetNbinsX(); ++i) {
      if (h_fit->GetBinContent(i) > maxContent * 0.15) {
        Ath_init = h_fit->GetBinCenter(i);
        break;
      }
    }

    // フィット範囲の自動設定
    fitMin = Ath_init - fitRangeWidthLow;
    fitMax = Ath_init + fitRangeWidthHigh;
    if (fitMin < h_fit->GetXaxis()->GetXmin()) fitMin = h_fit->GetXaxis()->GetXmin();

    // TF1 オブジェクトの生成

  } else {
    h_fit->GetXaxis()->SetRangeUser(fitMin, fitMax); // 一時的に範囲制限
    double maxContent = h_fit->GetMaximum();
    int maxBin = h_fit->GetMaximumBin();
    double Ath_init = h_fit->GetBinCenter(maxBin); 
    h_fit->GetXaxis()->UnZoom();                     // 制限を解除
  }
  TF1 *f_fit = new TF1("f_fit_exp_erf", ThresholdExpFitFunc, fitMin, fitMax, 5);
  f_fit->SetParNames("p0", "p1", "p2", "Ath", "sigma");

  // 初期値と探索制限（ParLimits）の設定
  // f_fit->SetParameters(maxContent * 2.0, 0.02, 0.0, Ath_init, 3.0);
  // f_fit->SetParLimits(0, 0.0, maxContent * 10.0);   // p0 > 0
  f_fit->SetParameters(maxContent * 2.0, 0.01, 0.0, Ath_init, 5.0); // sigmaの初期値を5.0 ADCに
  f_fit->SetParLimits(1, 0.0001, 0.5);            // p1 > 0 (右肩下がり)
  f_fit->SetParLimits(2, 0.0, maxContent * 0.5);   // p2 >= 0
  f_fit->SetParLimits(3, Ath_init - 15.0, Ath_init + 15.0); // Ath
  // f_fit->SetParLimits(4, 0.3, 12.0);              // sigma (ノイズ幅の現実的な範囲)
  // f_fit->SetParLimits(4, 1.0, 15.0);
  f_fit->SetParLimits(4, 0.5, 30.0);

  if(!(fitMin == -1 || fitMax == -1)){
    f_fit->SetParameters(maxContent * 2.0, 0.01, 0.0, (fitMin + fitMax) / 2.0, 5.0);
    f_fit->SetParLimits(3, fitMin, fitMax); // ★ Athがペデスタル側（fitMinより左）に行かないようロック
  }

  // フィット実行 (Q: Quiet, S: SaveResult, R: Range)
  TFitResultPtr fitRes = h_fit->Fit(f_fit, "Q S R N");

    // 結果の格納
  if (fitRes.Get() && fitRes->IsValid()) {
    res.isSuccess = true;
    res.p0       = f_fit->GetParameter(0);
    res.p1       = f_fit->GetParameter(1);
    res.p2       = f_fit->GetParameter(2);
    res.Ath      = f_fit->GetParameter(3);
    res.AthErr   = f_fit->GetParError(3);
    res.sigma    = f_fit->GetParameter(4);
    res.sigmaErr = f_fit->GetParError(4);

    TF1* f_draw = new TF1(Form("f_draw_%s", h_adc_in->GetName()), ThresholdExpFitFunc, fitMin, fitMax, 5);
    f_draw->SetParameters(res.p0 / (double)rebinFactor, // 高さを元に戻す
                          res.p1, 
                          res.p2 / (double)rebinFactor, // 高さを元に戻す
                          res.Ath, 
                          res.sigma);
    f_draw->SetLineColor(kRed); // 赤色に設定
    f_draw->SetLineWidth(2);
    res.f_draw = f_draw;
    h_adc_in->GetListOfFunctions()->Add(f_draw);
        
    if (fitRes->Ndf() > 0) res.chi2ndf = fitRes->Chi2() / fitRes->Ndf();
    else res.chi2ndf = -1.0;
  } else {
    res.isSuccess = false;
  }

  // 一時オブジェクトの破棄
  delete h_fit;
  delete f_fit;

  return res;
}


int main(int argc, char* argv[])
{
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./between_files  output.root Nenegry E1...En E1Nfiles...EnNfiles input1.root input2.root  ...  figure_path " << endl;
        return 1;
    }
  gROOT->SetBatch(kTRUE);

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;

  int Nenegry = atoi(argv[2]);
  double energy[Nenegry];
  int energy_files[Nenegry];
  for(int i=0;i<Nenegry;i++){
    energy[i] = atoi(argv[3+i]);
    if(energy[i]<1) energy[i] = 0.5;
    energy_files[i] = atoi(argv[3+Nenegry+i]);
  }


  int rawfilenum = argc - 2 - 2*Nenegry - 1;
  int irawfilenum = 3 + 2*Nenegry;


  double CycleID, TriggerID;
  vector<int> *CellID = nullptr;
  vector<int> *BCID = nullptr;
  vector<int> *HitTag = nullptr;
  vector<int> *GainTag = nullptr;
  vector<double> *HG_Charge = nullptr;
  vector<double> *LG_Charge = nullptr;
  vector<vector<double>> *Temperature = nullptr;

  // double cycleID_ssa, triggerID_ssa;
  // vector<int> *cellIDs_ssa = nullptr;
  // vector<int> *BCIDs_ssa = nullptr;
  // vector<int> *hitTags_ssa = nullptr;
  // vector<int> *gainTags_ssa = nullptr;
  // vector<double> *charges_ssa = nullptr;
  // vector<double> *times_ssa = nullptr;
  // vector<double> *temp_ssa = nullptr;
  // vector<double> *posX_ssa = nullptr;
  // vector<double> *posY_ssa = nullptr;
  // vector<double> *posZ_ssa = nullptr;
  // vector<int> *ssaTag_ssa = nullptr;

  // double EnergyDep;
  // vector<int> *cellID = nullptr;
  // vector<double> *Hit_Energy = nullptr;
  // vector<double> *Hit_X = nullptr;
  // vector<double> *Hit_Y = nullptr;
  // vector<double> *Hit_Z = nullptr;



  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int doublelayerNu = 2;
  const double _xInterval = 5.3;  // 300 um gap in width direction
  const double _yInterval = 45.4; // 400 um gap in length direction
  const int rowNu = 42;
  const int columnNu = 5;
  const double _zInitialPosition = 1.;
  const double _zIntervalEven = 12.7; // mm
  const double _zIntervalOdd = 19.9; // mm
  double _xMin = -rowNu*_xInterval/2.;
  double _xMax = rowNu*_xInterval/2.;
  double _yMin = -columnNu*_yInterval/2.;
  double _yMax = columnNu*_yInterval/2.;
  double _zMin = -10; double _zMax = 320; int _zBin = 165;
  double coormin = -(rowNu-1)*_xInterval/2.;
  double accurateBinY[10] = {-113.3,-68.3, -67.9,-22.9, -22.5,22.5, 22.9,67.9, 68.3,113.3};
  double accurateBinX[84] = {0};
  for(int i=0;i<rowNu;i++){
      accurateBinX[2*i] = -111.15 + i*_xInterval;
      accurateBinX[2*i+1] = accurateBinX[2*i] + 5.0;
  }

  const int range_maximum = 12000;
  const int range_maximum_PS = 500;
  const double binWidth_PS = 0.1;




  TFile fileout(argv[1],"RECREATE");
  // gain_histo->Write();


  cout<<"=========== PedestalExtract ============"<<endl;
    TFile* pedFile = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/pedestal2023_SPS.root");
    TTree* pedTree = (TTree*)pedFile->Get("ChnLevel"); 
    if(!pedTree)  {cout<<"!!! GET PED FILE FAILED !!!"<<endl;}
    vector<int>* _pedCellID = nullptr;
    vector<double>* _pedMeanHigh = nullptr;
    vector<double>* _pedSigmaHigh = nullptr;
    vector<double>* _pedMeanLow = nullptr;
    vector<double>* _pedSigmaLow = nullptr;
    pedTree ->SetBranchAddress("CellID",&_pedCellID);
    pedTree ->SetBranchAddress("PedHighMean",&_pedMeanHigh);
    pedTree ->SetBranchAddress("PedHighSig",&_pedSigmaHigh);
    pedTree ->SetBranchAddress("PedLowMean",&_pedMeanLow);
    pedTree ->SetBranchAddress("PedLowSig",&_pedSigmaLow);

    double meanPedHigh15  = 0 ; double meanPedHigh10  = 0;
    double sigmaPedHigh15 = 0 ; double sigmaPedHigh10 = 0;
    double chnNoPedHigh15 = 0 ; double chnNoPedHigh10 = 0;
    double meanPedLow15   = 0 ; double meanPedLow10   = 0;
    double sigmaPedLow15  = 0 ; double sigmaPedLow10  = 0;
    double chnNoPedLow15  = 0 ; double chnNoPedLow10  = 0;
    double pedMeanHigh[layerNu][chipNu][channelNu];
    double pedSigmaHigh[layerNu][chipNu][channelNu];
    double pedMeanLow[layerNu][chipNu][channelNu];
    double pedSigmaLow[layerNu][chipNu][channelNu];
    TH2F *ped_map[layerNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          pedMeanHigh[i_layer][i_chip][i_channel] = 1;
          pedSigmaHigh[i_layer][i_chip][i_channel] = 1;
          pedMeanLow[i_layer][i_chip][i_channel] = 1;
          pedSigmaLow[i_layer][i_chip][i_channel] = 1;
        }
      }
      double xNu_ = i_layer%2==0 ? columnNu : rowNu;
      double yNu_ = i_layer%2==0 ? rowNu : columnNu;
      double xMax_ = i_layer%2==0 ? _yMax : _xMax;
      double yMax_ = i_layer%2==0 ? _xMax : _yMax;
      // double accuratexNu_ = xNu_ * 2 - 1;
      // double accurateyNu_ = yNu_ * 2 - 1;
      // double* accurateXBin = i_layer%2==0 ? accurateBinY : accurateBinX;
      // double* accurateYBin = i_layer%2==0 ? accurateBinX : accurateBinY;
      ped_map[i_layer] = new TH2F(Form("ped_map_%d",i_layer),Form("layer %d",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
    }
    TH1F *pedestals[2][2];
    TH1F *pedestals_sigma[2][2];
    for(int pix=0;pix<2;pix++){
      for(int hl=0;hl<2;hl++){
        int pixel = pix == 0 ? 10 : 15;
        string HL = hl == 0 ? "high" : "low";
        pedestals[pix][hl] = new TH1F(Form("ped_%s_%d",HL.c_str(),pixel), Form("%s gain pedestal (%d um)",HL.c_str(),pixel), 300,300,600);
        pedestals_sigma[pix][hl] = new TH1F(Form("ped_sigma_%s_%d",HL.c_str(),pixel), Form("%s gain pedestal sigma (%d um)",HL.c_str(),pixel), 1000,0,10);
      }
    }
    for(int entry=0; entry!=pedTree->GetEntries(); ++entry){
      pedTree->GetEntry(entry);
      //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
      for(size_t ihit=0; ihit!=_pedCellID->size(); ++ihit)
      {
        int layerID = _pedCellID->at(ihit)/1e5;
        int chipID = _pedCellID->at(ihit)%100000/1e4;
        int channelID = _pedCellID->at(ihit)%100;        

        if (layerID>3 && layerID<28) {
          // if(_ChiSqr/(double)_NDF>2) continue;
          chnNoPedHigh10++;
          chnNoPedLow10++;
          meanPedHigh10 += _pedMeanHigh->at(ihit);
          meanPedLow10  += _pedMeanLow ->at(ihit);
          sigmaPedHigh10 += _pedSigmaHigh->at(ihit);
          sigmaPedLow10  += _pedSigmaLow ->at(ihit);
        } else {
           // if(_ChiSqr/(double)_NDF>1.7) continue;
          chnNoPedHigh15++;
          chnNoPedLow15++;
          meanPedHigh15 += _pedMeanHigh->at(ihit);
          meanPedLow15  += _pedMeanLow ->at(ihit);
          sigmaPedHigh15 += _pedSigmaHigh->at(ihit);
          sigmaPedLow15  += _pedSigmaLow ->at(ihit);
        }
        pedMeanHigh[layerID][chipID][channelID] = _pedMeanHigh->at(ihit);
        pedSigmaHigh[layerID][chipID][channelID] = _pedSigmaHigh->at(ihit);
        pedMeanLow[layerID][chipID][channelID] = _pedMeanLow->at(ihit);
        pedSigmaLow[layerID][chipID][channelID] = _pedSigmaLow->at(ihit);
      }
    }
    meanPedHigh15 /= chnNoPedHigh15;
    meanPedLow15  /= chnNoPedLow15;
    meanPedHigh10 /= chnNoPedHigh10;
    meanPedLow10  /= chnNoPedLow10;
    sigmaPedHigh15 /= chnNoPedHigh15;
    sigmaPedLow15  /= chnNoPedLow15;
    sigmaPedHigh10 /= chnNoPedHigh10;
    sigmaPedLow10  /= chnNoPedLow10;
    cout << " HIGH GAIN PEDESTAL ... ..."  << endl;
    cout << " channel count(15um:10um) : " << chnNoPedHigh15 << " , " << chnNoPedHigh10 << endl;
    cout << " mean of pedestal High : "    << meanPedHigh15  << " , " << meanPedHigh10  << endl;
    cout << " sigma of pedestal High: "    << sigmaPedHigh15 << " , " << sigmaPedHigh10 << endl;
    cout << " LOW GAIN PEDESTAL ... ..."   << endl;
    cout << " channel count(15um:10um) : " << chnNoPedLow15  << " , " << chnNoPedLow10  << endl;
    cout << " mean of pedestal Low : "     << meanPedLow15   << " , " << meanPedLow10   << endl;
    cout << " sigma of pedestal Low: "     << sigmaPedLow15  << " , " << sigmaPedLow10  << endl;

    for (int il=0; il<layerNu; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if(pedMeanHigh[il][ip][ic]==1) cout << il << ", " << ip << ", " << ic << endl;
          if (il>3 && il<28) {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh10 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow10 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh10 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow10 : pedSigmaLow[il][ip][ic];

            pedestals[0][0]->Fill(pedMeanHigh[il][ip][ic]);
            pedestals[0][1]->Fill(pedMeanLow[il][ip][ic]);
            pedestals_sigma[0][0]->Fill(pedSigmaHigh[il][ip][ic]);
            pedestals_sigma[0][1]->Fill(pedSigmaLow[il][ip][ic]);
          } else {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh15 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow15 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh15 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow15 : pedSigmaLow[il][ip][ic];

            if(il<30){
              pedestals[1][0]->Fill(pedMeanHigh[il][ip][ic]);
              pedestals[1][1]->Fill(pedMeanLow[il][ip][ic]);
              pedestals_sigma[1][0]->Fill(pedSigmaHigh[il][ip][ic]);
              pedestals_sigma[1][1]->Fill(pedSigmaLow[il][ip][ic]);
            }
          }

          if(ip==5 && ic>29) continue;
          double* position = EBUdecode_itr(il,ip,ic);
          double int_x = il%2==0 ? _yInterval : _xInterval;
          double int_y = il%2==0 ? _xInterval : _yInterval;
          double xMax_ = il%2==0 ? _yMax : _xMax;
          double yMax_ = il%2==0 ? _xMax : _yMax;
          int ix_ = (position[0]+xMax_+0.5)/int_x;
          int iy_ = (position[1]+yMax_+0.5)/int_y;
          // int x_ = il%2==0 ? ix_ : iy_;
          // int y_ = il%2==0 ? iy_ : ix_;

          ped_map[il]->SetBinContent(ix_+1,iy_+1, pedMeanHigh[il][ip][ic]);
        }
      }
    }
    if(pedestal_figure_saving){
      TFile fileout_pedestal("../result/pedestal/pedestal.root","RECREATE");
      fileout_pedestal.cd();
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_ped_map_per_layer%d", i_layer), Form("ped map Layer%d", i_layer), 2560, 2560);
        ped_map[i_layer]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("../result/pedestal/figures/pedestal_map_Layer%d.png",i_layer));
        ped_map[i_layer]->Write(Form("pedestal_map_Layer%d", i_layer));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_ped_map_layer%d", i), Form("ped map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          ped_map[index]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("../result/pedestal/figures/ped_map%d.png",i));
        c1->Write(Form("ped_map%d",i));
        delete c1;
      }
      for(int i=0;i<2;i++){
        int pixel = i == 0 ? 10 : 15;
        TCanvas* c1 = new TCanvas(Form("c_ped_%d", i), Form("ped %d", i), 2560, 1440);
        c1->Divide(2,2);
        c1->cd();
        for(int ii=0;ii<2;ii++){
          c1->cd(ii+1);
          gStyle->SetOptStat(1);
          pedestals[i][ii]->Draw();
          c1->cd(ii+3);
          pedestals_sigma[i][ii]->Draw();
        }
        c1->Update();
        c1->SaveAs(Form("../result/pedestal/figures/ped_%dum.png",pixel));
        c1->Write(Form("ped_%dum.png",pixel));
        delete c1;
      }
    }
  pedFile->Close();
  fileout.cd();




  // dataを入れるもの
  

  cout << "made save TTree" << endl;
  TTree *treeout = new TTree("threshold", "Threshold Fit Results");
  int layer_out, chip_out, channel_out, CellID_out;
  double Ath, AthErr, sigma, sigmaErr, chi2ndf, p0, p1, p2;
  treeout->Branch("CellID",       &CellID_out,  "CellID/I");
  treeout->Branch("layer",        &layer_out,   "layer/I");
  treeout->Branch("chip",         &chip_out,    "chip/I");
  treeout->Branch("channel",      &channel_out, "channel/I");
  treeout->Branch("threshold",    &Ath,         "Ath/D");
  treeout->Branch("thresholdErr", &AthErr,      "AthErr/D");
  treeout->Branch("sigma",        &sigma,       "sigma/D");
  treeout->Branch("sigmaErr",     &sigmaErr,    "sigmaErr/D");
  treeout->Branch("chi2ndf",      &chi2ndf,     "chi2ndf/D");
  treeout->Branch("p0",           &p0,          "p0/D");
  treeout->Branch("p1",           &p1,          "p1/D");
  treeout->Branch("p2",           &p2,          "p2/D");


  TH1F* all_adc_channel[layerNu][chipNu][channelNu];
  TH1F* all_adc_channel_raw[layerNu][chipNu][channelNu];
  for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
    for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
      for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
        all_adc_channel[i_layer][i_chip][i_channel] = new TH1F(Form("all_adc_channel_%d_%d_%d",i_layer, i_chip, i_channel),Form("layer%d chip%d channel%d", i_layer, i_chip, i_channel),500,20,520);
        all_adc_channel[i_layer][i_chip][i_channel]->SetXTitle("[ADC]");
        all_adc_channel[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.075);
        all_adc_channel[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.075);
        all_adc_channel_raw[i_layer][i_chip][i_channel] = new TH1F(Form("all_adc_channel_raw_%d_%d_%d",i_layer, i_chip, i_channel),Form("layer%d chip%d channel%d", i_layer, i_chip, i_channel),500,300,800);
        all_adc_channel_raw[i_layer][i_chip][i_channel]->SetXTitle("[ADC]");
        all_adc_channel_raw[i_layer][i_chip][i_channel]->GetXaxis()->SetLabelSize(0.075);
        all_adc_channel_raw[i_layer][i_chip][i_channel]->GetYaxis()->SetLabelSize(0.075);
      }
    }
  }



  // data をとってきてる
  for(int ienergy=0; ienergy<Nenegry; ienergy++){
    if(Nenegry>1) cout << energy[ienergy] << " GeV" << endl;

    TH1F* adc_channel[layerNu][chipNu][channelNu];
    TH1F* adc_channel_raw[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          adc_channel[i_layer][i_chip][i_channel] = new TH1F(Form("adc_channel_%d_%d_%d",i_layer, i_chip, i_channel),Form("%g GeV layer%d chip%d channel%d",energy[ienergy], i_layer, i_chip, i_channel),500,-50,450);
          adc_channel_raw[i_layer][i_chip][i_channel] = new TH1F(Form("adc_channel_raw_%d_%d_%d",i_layer, i_chip, i_channel),Form("%g GeV layer%d chip%d channel%d",energy[ienergy], i_layer, i_chip, i_channel),500,300,800);
        }
      }
    }

    TFile *filein[energy_files[ienergy]];
    TTree *tree[energy_files[ienergy]];
    int entry_max[energy_files[ienergy]];
    for(int i=0; i<energy_files[ienergy]; i++){
      cout << argv[irawfilenum+i] << endl;
      filein[i] = new TFile(argv[irawfilenum+i]);
      tree[i] = (TTree*) filein[i]->Get("Raw_Hit");
      entry_max[i] = tree[i]->GetEntries();
    }

    for(int irawfile=0; irawfile<energy_files[ienergy]; irawfile++){
      cout << irawfile << "/" << energy_files[ienergy] << endl;

      tree[irawfile]->SetBranchAddress("CellID", &CellID);
      tree[irawfile]->SetBranchAddress("HG_Charge", &HG_Charge);
      tree[irawfile]->SetBranchAddress("LG_Charge", &LG_Charge);
      tree[irawfile]->SetBranchAddress("HitTag", &HitTag);

      for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
        tree[irawfile]->GetEntry(ientry);

        for(int ihit=0;ihit<CellID->size();ihit++){
          if(HitTag->at(ihit)==0) continue;
          int _layerID = CellID->at(ihit)/1e5;
          int _chipID  = (CellID->at(ihit)%100000) /1e4;
          int _chanID  = CellID->at(ihit)%100;
          if(_layerID>=30) continue;
          if( (_chipID<0 || _chipID>5) || (_chanID<0 || _chanID>35) ) cout << "error in cellId " << CellID->at(ihit) << endl;

          adc_channel[_layerID][_chipID][_chanID]->Fill(HG_Charge->at(ihit) - pedMeanHigh[_layerID][_chipID][_chanID]);
          adc_channel_raw[_layerID][_chipID][_chanID]->Fill(HG_Charge->at(ihit));
          all_adc_channel[_layerID][_chipID][_chanID]->Fill(HG_Charge->at(ihit) - pedMeanHigh[_layerID][_chipID][_chanID]);
          all_adc_channel_raw[_layerID][_chipID][_chanID]->Fill(HG_Charge->at(ihit));
        }
      }

      CellID = nullptr;
      HG_Charge = nullptr;
      LG_Charge = nullptr;
      HitTag = nullptr;
      delete filein[irawfile];
    }

    fileout.cd();
    if(Nenegry == 1){
      cout << energy[ienergy] << "GeV" << endl;
      fileout.mkdir(Form("%dGeV", (int)energy[ienergy]));
      fileout.cd(Form("%dGeV", (int)energy[ienergy]));
      std::ofstream outfile(Form("txt/threshold_results_%gGeV.txt",energy[ienergy]));
      outfile << "channel/I:Ath/D:AthErr/D:sigma/D:sigmaErr/D:chi2ndf/D" << std::endl;
      TH2F *thre_map[layerNu][2];
      TH1F *thresholds[2];
      TH1F *thresholds_sigma[2];
      for(int pix=0;pix<2;pix++){
        int pixel = pix == 0 ? 10 : 15;
        thresholds[pix] = new TH1F(Form("thre_%d",pixel), Form("threshold (%d um)",pixel), 300,0,300);
        thresholds_sigma[pix] = new TH1F(Form("thre_sigma_%d",pixel), Form("threshold sigma (%d um)",pixel), 3000,0,30);
      }
      for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
        double xNu_ = i_layer%2==0 ? columnNu : rowNu;
        double yNu_ = i_layer%2==0 ? rowNu : columnNu;
        double xMax_ = i_layer%2==0 ? _yMax : _xMax;
        double yMax_ = i_layer%2==0 ? _xMax : _yMax;
        thre_map[i_layer][0] = new TH2F(Form("thre_map_%d",i_layer),Form("layer %d threshold",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        thre_map[i_layer][1] = new TH2F(Form("threSigma_map_%d",i_layer),Form("layer %d threshold sigma",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
        if(i_layer>=30) continue;

        for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
          TCanvas* c1 = new TCanvas(Form("adc_Layer%d_Chip%d", i_layer, i_chip), Form("Layer%d Chip%d", i_layer, i_chip), 2560, 1440);
          c1->Divide(6, 6);
          gStyle->SetOptFit(1111);
          gStyle->SetStatFontSize(0.075);
          for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
            if(i_chip==5 && i_channel>29) continue;
            c1->cd();
            c1->cd(i_channel+1);
            adc_channel[i_layer][i_chip][i_channel]->Draw();
            ThresholdFitResult result = FitChannelThreshold(adc_channel[i_layer][i_chip][i_channel], 10, 15, 40, pedSigmaHigh[i_layer][i_chip][i_channel]*5, pedSigmaHigh[i_layer][i_chip][i_channel]*5+300);
            if (result.isSuccess) {
              if (result.f_draw) result.f_draw->Draw("same"); // ★ 赤線を重ねて描画
              std::cout << Form("Layer %2d Chip %2d Ch %2d: Ath = %.2f +- %.2f ADC, sigma = %.2f +- %.2f ADC (chi2/ndf = %.2f)", i_layer, i_chip, i_channel, result.Ath, result.AthErr, result.sigma, result.sigmaErr, result.chi2ndf) << std::endl;
              outfile << i_layer <<"," << i_chip << "," << i_channel << " " << result.Ath << " " << result.AthErr << " " << result.sigma << " " << result.sigmaErr << " " << result.chi2ndf << std::endl;
              
              CellID_out    = i_layer * 1e5  +  i_chip * 1e4  + i_channel;
              layer_out    = i_layer;
              chip_out     = i_chip;
              channel_out  = i_channel;
              Ath      = result.Ath;
              AthErr   = result.AthErr;
              sigma    = result.sigma;
              sigmaErr = result.sigmaErr;
              chi2ndf  = result.chi2ndf;
              p0       = result.p0;
              p1       = result.p1;
              p2       = result.p2;
              treeout->Fill();

              double* position = EBUdecode_itr(i_layer,i_chip,i_channel);
              double int_x = i_layer%2==0 ? _yInterval : _xInterval;
              double int_y = i_layer%2==0 ? _xInterval : _yInterval;
              double xMax_ = i_layer%2==0 ? _yMax : _xMax;
              double yMax_ = i_layer%2==0 ? _xMax : _yMax;
              int ix_ = (position[0]+xMax_+0.5)/int_x;
              int iy_ = (position[1]+yMax_+0.5)/int_y;
              thre_map[i_layer][0]->SetBinContent(ix_+1,iy_+1, result.Ath);
              thre_map[i_layer][1]->SetBinContent(ix_+1,iy_+1, result.sigma);

              int pixel = (i_layer>3 && i_layer<28) ? 0 : 1;
              thresholds[pixel]->Fill(result.Ath);
              thresholds_sigma[pixel]->Fill(result.sigma);
            } else {
              std::cerr << Form("Layer %2d Chip %2d Ch %2d: Fit Failed!", i_layer, i_chip, i_channel) << std::endl;
            }
          }
          c1->Update();
          c1->Write(Form("adc_Layer%d_Chip%d", i_layer, i_chip));
          c1->SaveAs(Form("%s/%gGeV/%gGeV_layer%dchip%d.png",argv[argc-1],energy[ienergy],energy[ienergy],i_layer,i_chip));
          delete c1;
        }
      }
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
          TCanvas* c1 = new TCanvas(Form("raw_adc_Layer%d_Chip%d", i_layer, i_chip), Form("Layer%d Chip%d", i_layer, i_chip), 2560, 1440);
          c1->Divide(6, 6);
          for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
            c1->cd();
            c1->cd(i_channel+1);
            adc_channel_raw[i_layer][i_chip][i_channel]->Draw();
          }
          c1->Update();
          c1->Write(Form("raw_adc_Layer%d_Chip%d", i_layer, i_chip));
          c1->SaveAs(Form("%s/%gGeV/%gGeV_raw_layer%dchip%d.png",argv[argc-1],energy[ienergy],energy[ienergy],i_layer,i_chip));
          delete c1;
        }
      }
      for(int meansigma=0;meansigma<2;meansigma++){
        string MS = meansigma==0 ? "": "Sigma";
        string MSname = meansigma==0 ? "": "_sigma";
        for (int i_layer = 0; i_layer < 30; ++i_layer) {
          TCanvas* c1 = new TCanvas(Form("c_thre_map_per_layer%d", i_layer), Form("thre map Layer%d", i_layer), 2560, 2560);
          thre_map[i_layer][meansigma]->Draw("colz");
          gStyle->SetOptStat(0);
          c1->Update();
          c1->SaveAs(Form("%s/%gGeV/%gGeV_thre%sMap_layer%d.png",argv[argc-1],energy[ienergy],energy[ienergy],MS.c_str(),i_layer));
          thre_map[i_layer][meansigma]->Write(Form("threshold_map%s_Layer%d", MSname.c_str(),i_layer));
          delete c1;
        }
        for(int i=0;i<2;i++){
          TCanvas* c1 = new TCanvas(Form("c_thre_map_layer%d", i), Form("thre map Layer%d", i), 2560, 2560);
          c1->Divide(4,4);
          c1->cd();
          for(int _layer=0;_layer<16;_layer++){
            c1->cd(_layer+1);
            int index = i*16 + _layer;
            gStyle->SetOptStat(0);
            thre_map[index][meansigma]->Draw("colz"); 
          }
          c1->Update();
          c1->SaveAs(Form("%s/%gGeV/%gGeV_thre%s_map%d.png",argv[argc-1],energy[ienergy],energy[ienergy],MSname.c_str(),i));
          c1->Write(Form("thre%s_map%d",MSname.c_str(),i));
          delete c1;
        }
      }
      if(1){
        TCanvas* c1 = new TCanvas(Form("c_thre"), Form("thre"), 2560, 1440);
        c1->Divide(2,2);
        for(int i=0;i<2;i++){
          int pixel = i == 0 ? 10 : 15;
          c1->cd();
          c1->cd(i+1);
          gStyle->SetOptStat(1);
          thresholds[i]->Draw();
          c1->cd(i+3);
          thresholds_sigma[i]->Draw();
        }
        c1->Update();
        c1->SaveAs(Form("%s/%gGeV/%gGeV_thre_um.png",argv[argc-1],energy[ienergy],energy[ienergy]));
        c1->Write(Form("thre_um.png"));
        delete c1;
      }
      
    }
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          delete adc_channel[i_layer][i_chip][i_channel];
          delete adc_channel_raw[i_layer][i_chip][i_channel];
        }
      }
    }
    fileout.cd();
    irawfilenum += energy_files[ienergy];

    
  }

  if(Nenegry > 1){
    fileout.cd();
    std::ofstream outfile(Form("txt/threshold_results.txt"));
    outfile << "channel/I:Ath/D:AthErr/D:sigma/D:sigmaErr/D:chi2ndf/D" << std::endl;
    TH2F *thre_map[layerNu][2];
    TH1F *thresholds[2];
    TH1F *thresholds_sigma[2];
    for(int pix=0;pix<2;pix++){
      int pixel = pix == 0 ? 10 : 15;
      thresholds[pix] = new TH1F(Form("ped_%d",pixel), Form("threshold (%d um)",pixel), 300,0,300);
      thresholds_sigma[pix] = new TH1F(Form("ped_sigma_%d",pixel), Form("threshold sigma (%d um)",pixel), 3000,0,30);
    }
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      double xNu_ = i_layer%2==0 ? columnNu : rowNu;
      double yNu_ = i_layer%2==0 ? rowNu : columnNu;
      double xMax_ = i_layer%2==0 ? _yMax : _xMax;
      double yMax_ = i_layer%2==0 ? _xMax : _yMax;
      thre_map[i_layer][0] = new TH2F(Form("thre_map_%d",i_layer),Form("layer %d threshold",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
      thre_map[i_layer][1] = new TH2F(Form("threSigma_map_%d",i_layer),Form("layer %d threshold sigma",i_layer), xNu_,-xMax_,xMax_, yNu_,-yMax_,yMax_);
      if(i_layer>=30) continue;

      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        TCanvas* c1 = new TCanvas(Form("adc_Layer%d_Chip%d", i_layer, i_chip), Form("Layer%d Chip%d", i_layer, i_chip), 2560, 1440);
        c1->Divide(6, 6);
        gStyle->SetOptFit(1111);
        gStyle->SetStatFontSize(0.075);
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          if(i_chip==5 && i_channel>29) continue;
          c1->cd();
          c1->cd(i_channel+1);
          all_adc_channel[i_layer][i_chip][i_channel]->Draw();
          ThresholdFitResult result = FitChannelThreshold(all_adc_channel[i_layer][i_chip][i_channel], 10, 15, 40, pedSigmaHigh[i_layer][i_chip][i_channel]*5, pedSigmaHigh[i_layer][i_chip][i_channel]*5+300);
          if (result.isSuccess) {
            if (result.f_draw) result.f_draw->Draw("same"); // ★ 赤線を重ねて描画
            std::cout << Form("Layer %2d Chip %2d Ch %2d: Ath = %.2f +- %.2f ADC, sigma = %.2f +- %.2f ADC (chi2/ndf = %.2f)", i_layer, i_chip, i_channel, result.Ath, result.AthErr, result.sigma, result.sigmaErr, result.chi2ndf) << std::endl;
            outfile << i_layer <<"," << i_chip << "," << i_channel << " " << result.Ath << " " << result.AthErr << " " << result.sigma << " " << result.sigmaErr << " " << result.chi2ndf << std::endl;
            
            CellID_out    = i_layer * 1e5  +  i_chip * 1e4  + i_channel;
            layer_out    = i_layer;
            chip_out     = i_chip;
            channel_out  = i_channel;
            Ath      = result.Ath;
            AthErr   = result.AthErr;
            sigma    = result.sigma;
            sigmaErr = result.sigmaErr;
            chi2ndf  = result.chi2ndf;
            p0       = result.p0;
            p1       = result.p1;
            p2       = result.p2;
            treeout->Fill();

            double* position = EBUdecode_itr(i_layer,i_chip,i_channel);
            double int_x = i_layer%2==0 ? _yInterval : _xInterval;
            double int_y = i_layer%2==0 ? _xInterval : _yInterval;
            double xMax_ = i_layer%2==0 ? _yMax : _xMax;
            double yMax_ = i_layer%2==0 ? _xMax : _yMax;
            int ix_ = (position[0]+xMax_+0.5)/int_x;
            int iy_ = (position[1]+yMax_+0.5)/int_y;
            thre_map[i_layer][0]->SetBinContent(ix_+1,iy_+1, result.Ath);
            thre_map[i_layer][1]->SetBinContent(ix_+1,iy_+1, result.sigma);

            int pixel = (i_layer>3 && i_layer<28) ? 0 : 1;
            thresholds[pixel]->Fill(result.Ath);
            thresholds_sigma[pixel]->Fill(result.sigma);
          } else {
            std::cerr << Form("Layer %2d Chip %2d Ch %2d: Fit Failed!", i_layer, i_chip, i_channel) << std::endl;
          }
        }
        c1->Update();
        c1->Write(Form("adc_Layer%d_Chip%d", i_layer, i_chip));
        c1->SaveAs(Form("%s/layer%dchip%d.png",argv[argc-1],i_layer,i_chip));
        delete c1;
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          delete all_adc_channel[i_layer][i_chip][i_channel];
        }
      }
    }
    for (int i_layer = 0; i_layer < 30; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        TCanvas* c1 = new TCanvas(Form("raw_adc_Layer%d_Chip%d", i_layer, i_chip), Form("Layer%d Chip%d", i_layer, i_chip), 2560, 1440);
        c1->Divide(6, 6);
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          c1->cd();
          c1->cd(i_channel+1);
          all_adc_channel_raw[i_layer][i_chip][i_channel]->Draw();
        }
        c1->Update();
        c1->Write(Form("raw_adc_Layer%d_Chip%d", i_layer, i_chip));
        c1->SaveAs(Form("%s/raw_layer%dchip%d.png",argv[argc-1],i_layer,i_chip));
        delete c1;
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          delete all_adc_channel_raw[i_layer][i_chip][i_channel];
        }
      }
    }
    for(int meansigma=0;meansigma<2;meansigma++){
      string MS = meansigma==0 ? "": "Sigma";
      string MSname = meansigma==0 ? "": "_sigma";
      for (int i_layer = 0; i_layer < 30; ++i_layer) {
        TCanvas* c1 = new TCanvas(Form("c_thre_map_per_layer%d", i_layer), Form("thre map Layer%d", i_layer), 2560, 2560);
        thre_map[i_layer][meansigma]->Draw("colz");
        gStyle->SetOptStat(0);
        c1->Update();
        c1->SaveAs(Form("%s/thre%sMap_layer%d.png",argv[argc-1],MS.c_str(),i_layer));
        thre_map[i_layer][meansigma]->Write(Form("threshold_map%s_Layer%d", MSname.c_str(),i_layer));
        delete c1;
      }
      for(int i=0;i<2;i++){
        TCanvas* c1 = new TCanvas(Form("c_thre_map_layer%d", i), Form("thre map Layer%d", i), 2560, 2560);
        c1->Divide(4,4);
        c1->cd();
        for(int _layer=0;_layer<16;_layer++){
          c1->cd(_layer+1);
          int index = i*16 + _layer;
          gStyle->SetOptStat(0);
          thre_map[index][meansigma]->Draw("colz"); 
        }
        c1->Update();
        c1->SaveAs(Form("%s/thre%s_map%d.png",argv[argc-1],MSname.c_str(),i));
        c1->Write(Form("thre%s_map%d",MSname.c_str(),i));
        delete c1;
      }
    }
    if(1){
      TCanvas* c1 = new TCanvas(Form("c_thre"), Form("thre"), 2560, 1440);
      c1->Divide(2,2);
      for(int i=0;i<2;i++){
        int pixel = i == 0 ? 10 : 15;
        c1->cd();
        c1->cd(i+1);
        gStyle->SetOptStat(1);
        thresholds[i]->Draw();
        c1->cd(i+3);
        thresholds_sigma[i]->Draw();
      }
      c1->Update();
      c1->SaveAs(Form("%s/thre_um.png",argv[argc-1]));
      c1->Write(Form("thre_um.png"));
      delete c1;
    }
  }

  fileout.cd();
  treeout->Write();









}