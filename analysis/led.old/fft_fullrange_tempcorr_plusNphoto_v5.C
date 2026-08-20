// fft_fullrange_tempcorr_plusNphoto_v5.C
//  - SingleHitのMPV（MIP_FIT.LandauMPV）を使用
//  - 旧MIP_Fit(Layer/Chip/Chn)にもフォールバック対応
//  - 壊れたROOT/欠損ファイルを安全にスキップ（zlib/TBasketの大量メッセージは抑制）
//  - 低統計/ヒスト欠落のチャンネルはサマリのみ表示
//  - Chi2/NDF（MPV fit）による品質判定追加
//  - pedSigma<50 のスプリット出力
//  - Nphoto も 2D マップ保存（レイヤ別）

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <memory>
#include <cfloat>
#include <cstring>

#include "TROOT.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "TTree.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TError.h"   // gErrorIgnoreLevel

namespace cfg {

  // ---- pedestal（平均のみ使用。無ければフォールバック） ----
  const char* PEDFILE = "/megraid01/users/takatsu_t/beamtest/ECAL_crstk_simulation/Analysis/share/pedestal2023_SPS.root";

  // ---- 入出力 ----
  const char* INDIR_HIST   = "/megraid01/users/data_beamtest/disk-BeamTestData1/Analyzed/mip/hist/mipScan/";
  const char* INDIR_DECODE = "/megraid01/users/data_beamtest/disk-BeamTestData1/Analyzed/decode/mipScan/";
  const char* OUTDIR       = "/megraid01/users/takatsu_t/Analysis/ScECAL_KEK/Results_gainMaps_tempCorr";

  // MPV (MIPFit) の場所（各 Run ごと）
  const char* MIPFIT_BASE  = "/megraid01/users/takatsu_t/Analysis/ScECAL_KEK/Results/MIPFit";

  // ---- ped sharp（緩和）----
  const char* PED_SHARP_FILE =
    "/megraid01/users/data_beamtest/ECAL_data/analysed/2025/kek/ped/roots_0717/hist/Pedestal/"
    "layer_except_10_11_14_15_20_21_26_27/ECAL_Run170_20250605_155655.root";
  constexpr double PED_SIGMA_MAX      = 200.0;
  constexpr double PED_REL_SIGMA_MAX  = 0.50;

  // ---- 幾何 ----
  constexpr int LMIN = 0;
  constexpr int LMAX = 31;
  constexpr int NCHIP_PER_LAYER = 6; // 0..5
  constexpr int CHIP_MIN = 0;
  constexpr int CHIP_MAX = NCHIP_PER_LAYER-1;
  constexpr int CH0  = 0;
  constexpr int CH1  = 35;
  constexpr int NCHIP_AXIS = NCHIP_PER_LAYER;

  // ---- しきい ----
  constexpr int    MIN_ENT        = 2000;    // 入力ヒスト最低統計（SPE 側）
  constexpr int    MIP_MIN_ENT    = 2000;    // MPV 側の統計しきい（SingleHitでも同値）
  constexpr double SNR_MIN_REF    = 5.0;     // 参考用
  constexpr double RELERR_MAX_REF = 0.15;    // 参考用
  constexpr double MPV_RCHI2_MAX  = 8.0;     // MPVの reduced-chi2 カット（Chi2/NDF）

  // ---- FFT ----
  constexpr bool   PAD_TO_POW2 = true;
  constexpr int    PAD_MAX     = 1<<15;
  constexpr double PED_FALLBACK = 400.0;

  // ---- 温度補正 ----
  constexpr double TREF_FIXED_C       = 25.0;
  constexpr double ALPHA_L15 = -0.020; // 15µm (L0-3)
  constexpr double ALPHA_L10 = -0.017; // 10µm (L4-31)

  // ---- 出力上書き ----
  constexpr bool   OVERWRITE    = true;

  // 入力ヒストの名前
  const char* HDIR_SINGLE = "SingleHit/layer%d/chip%d";
  const char* HDIR_ALL    = "AllHits/layer%d/chip%d";
  const char* HNAME_SINGLE = "Layer%d_Chip%d_Chn%02d_1hit";
  const char* HNAME_ALL    = "Layer%d_Chip%d_Chn%02d_all";

  // 期待ゲイン帯域（ADC/p.e.）
  constexpr double GMIN_L15 = 5.0, GMAX_L15 = 30.0; // L0-3
  constexpr double GMIN_L10 = 5.0, GMAX_L10 = 30.0; // L4-31
  constexpr double GMIN_DEF =  3.0, GMAX_DEF = 60.0;

  // 可視化
  constexpr bool   WRITE_GAIN_HIST   = false;

  // ---- ログ制御 ----
  constexpr bool VERBOSE_MISSING_HIST = false;   // 個別missingは抑制
  constexpr bool VERBOSE_LOWSTAT      = false;   // 個別low statは抑制
}

// ---------- ユーティリティ ----------
inline int nextPow2(int n){ int p=1; while(p<n) p<<=1; return p; }
inline bool isFinitePos(double x){ return std::isfinite(x) && (x>0.0); }
inline double gainErr(double Npad,double k0){ const double sigk=0.5/std::sqrt(12.0); return (k0>0?(Npad*sigk)/(k0*k0):-1); }

template<class V> double median(V v){
  if(v.empty()) return 0.0;
  size_t n=v.size()/2; std::nth_element(v.begin(),v.begin()+n,v.end());
  double m=v[n];
  if(v.size()%2==0){ std::nth_element(v.begin(),v.begin()+n-1,v.end()); m=0.5*(m+v[n-1]); }
  return m;
}

inline void setHistStyle(TH1* h){ if(!h) return; h->SetOption("HIST"); h->SetFillStyle(0); h->SetLineWidth(2); }

inline void bandByLayer(int L,double &gmin,double &gmax){
  if(L>=0 && L<=3){ gmin=cfg::GMIN_L15; gmax=cfg::GMAX_L15; return; }
  if(L>=4 && L<=31){ gmin=cfg::GMIN_L10; gmax=cfg::GMAX_L10; return; }
  gmin=cfg::GMIN_DEF; gmax=cfg::GMAX_DEF;
}

// ---- スコープでROOTの冗長エラー抑制 ----
struct ErrorSilencer {
  int old;
  ErrorSilencer(int newLevel): old(gErrorIgnoreLevel) { gErrorIgnoreLevel = newLevel; }
  ~ErrorSilencer(){ gErrorIgnoreLevel = old; }
};

// ---------- FFT & ピークピック ----------
TH1D* doFFT(const std::vector<double>& vin,int Npad,const char* name){
  std::vector<double> v(Npad,0.0);
  for(size_t i=0;i<vin.size() && i<(size_t)Npad;++i) v[i]=vin[i];
  Int_t n=Npad;
  TVirtualFFT* fft=TVirtualFFT::FFT(1,&n,"R2C ES K");
  fft->SetPoints(v.data()); fft->Transform();
  std::vector<double> re(Npad/2+1), im(Npad/2+1);
  fft->GetPointsComplex(re.data(), im.data());
  delete fft;

  TH1D* h=new TH1D(name,name,Npad/2,0,Npad/2);
  for(int k=0;k<Npad/2;++k) h->SetBinContent(k+1,std::hypot(re[k],im[k]));
  setHistStyle(h);
  return h;
}

struct KPick{ double k0=-1,gain=-1,err=-1; int kp=-1; double snr=-1; };

KPick pickK0(const TH1D* h,int Npad,int L){
  KPick r; if(!h) return r;
  double gmin,gmax; bandByLayer(L,gmin,gmax);
  int kmin=std::max(1,(int)std::floor(Npad/gmax));
  int kmax=std::min(h->GetNbinsX(),(int)std::ceil(Npad/gmin));
  if(kmax<kmin) std::swap(kmax,kmin);
  int kp=kmin; double pk=h->GetBinContent(kp);
  for(int k=kmin+1;k<=kmax;++k){ double a=h->GetBinContent(k); if(a>pk){pk=a; kp=k;} }
  r.kp=kp;
  double w=0, kw=0;
  for(int dk=-1; dk<=+1; ++dk){
    int kk=kp+dk; if(kk<1||kk>h->GetNbinsX()) continue;
    double a=h->GetBinContent(kk);
    w+=a; kw+=kk*a;
  }
  r.k0  =(w>0?kw/w:kp);
  r.gain=(r.k0>0?Npad/r.k0:-1);
  r.err =(r.k0>0?gainErr(Npad,r.k0):-1);

  std::vector<double> noise; noise.reserve(kmax-kmin+1);
  for(int k=kmin;k<=kmax;++k){ if(std::abs(k-kp)<=2) continue; noise.push_back(h->GetBinContent(k)); }
  double med=median(noise);
  r.snr=(med>0?pk/med:-1);
  return r;
}

// ============ pedestal（平均のみ使用。無ければフォールバック） ============
struct PedDB{
  TH2D *hM=nullptr,*hS=nullptr;
  bool load(const char* fn){
    TFile f(fn,"READ"); if(f.IsZombie()) return false;
    hM=(TH2D*)f.Get("high_gain_mean");
    hS=(TH2D*)f.Get("high_gain_rms");
    if(!hM||!hS) return false;
    hM->SetDirectory(0); hS->SetDirectory(0);
    return true;
  }
  bool get(int L,int chip,int ch,double &mean,double &rms) const{
    if(!hM||!hS) { mean=cfg::PED_FALLBACK; rms=0; return true; }
    int x=L*cfg::NCHIP_PER_LAYER+chip+1;
    int y=ch+1;
    if(x<1||x>hM->GetNbinsX()||y<1||y>hM->GetNbinsY()) { mean=cfg::PED_FALLBACK; rms=0; return true; }
    mean=hM->GetBinContent(x,y);
    rms =hS->GetBinContent(x,y);
    if(mean<=0){ mean=cfg::PED_FALLBACK; }
    return true;
  }
}; static PedDB gPed;

// ---- pedestal sharpness（緩和版）----
struct PedSharpReader {
  TFile* f=nullptr; TDirectory* dir=nullptr; bool ready=false;
  bool open(const char* path){
    if(ready) return true;
    f=TFile::Open(path,"READ");
    if(!f || f->IsZombie()){ std::cerr<<"Error in <Ped>: cannot open "<<path<<"\n"; return false; }
    dir = dynamic_cast<TDirectory*>(f->Get("ChannelLevel"));
    if(!dir) dir=f;
    ready=true; return true;
  }
  bool getHG(int L,int chip,int ch,double &mean,double &sigma){
    if(!ready) return false;
    TString hname = Form("HighGain_layer%d_chip%d_chan%d", L, chip, ch);
    TH1* h = dynamic_cast<TH1*>(dir->Get(hname));
    if(!h) return false;
    TF1* fit = nullptr;
    if(h->GetListOfFunctions()){
      TIter nx(h->GetListOfFunctions()); TObject* o;
      while((o=nx())){
        TF1* f1 = dynamic_cast<TF1*>(o);
        if(!f1) continue;
        if(TString(f1->GetName()).Contains("gaus",TString::kIgnoreCase)){ fit=f1; break; }
        if(!fit) fit=f1;
      }
    }
    if(fit && fit->GetNpar()>=3){
      mean = fit->GetParameter(1);
      sigma= std::abs(fit->GetParameter(2));
    }else{
      mean = h->GetMean();
      sigma= h->GetRMS();
    }
    return true;
  }
}; static PedSharpReader gPedSharp;

inline bool pedIsSharp_relaxed(int L,int chip,int ch, double &mean,double &sigma){
  mean=0; sigma=1e9;
  if(!gPedSharp.open(cfg::PED_SHARP_FILE)) return false;
  if(!gPedSharp.getHG(L, chip, ch, mean, sigma)) return false;
  bool ok=true;
  if(sigma<=0) ok=false;
  if(cfg::PED_SIGMA_MAX>0 && sigma>cfg::PED_SIGMA_MAX) ok=false;
  if(mean>0 && cfg::PED_REL_SIGMA_MAX>0 && (sigma/mean)>cfg::PED_REL_SIGMA_MAX) ok=false;
  return ok;
}

// --- 再帰 mkdir（出力ファイル内）---
TDirectory* ensureDir(TFile* fout, const TString& fullpath){
  if(!fout) return nullptr;
  TDirectory* cur = fout;
  std::unique_ptr<TObjArray> toks(fullpath.Tokenize("/"));
  for (int i=0;i<toks->GetEntries();++i){
    TString name = ((TObjString*)toks->At(i))->GetString();
    if (name.IsNull()) continue;
    TDirectory* next = dynamic_cast<TDirectory*>(cur->Get(name));
    if (!next) next = cur->mkdir(name);
    cur = next;
  }
  return cur;
}

enum ESource { kSingle, kAll };

// ---- 1 run 内の集計カウンタ ----
struct RunCounters {
  long missing_single=0, missing_all=0;
  long lowstat_single=0,  lowstat_all=0;
};

// 入力ヒスト取得（chip ディレクトリ対応）
TH1D* getRawHist(TFile* fin, int L, int chip, int ch, ESource src){
  if(!fin) return nullptr;
  const char* dpat = (src==kSingle)? cfg::HDIR_SINGLE : cfg::HDIR_ALL;
  const char* hpat1= (src==kSingle)? cfg::HNAME_SINGLE: cfg::HNAME_ALL;

  TString dpath = Form(dpat, L, chip);
  TDirectory* d = dynamic_cast<TDirectory*>(fin->Get(dpath));
  if(!d) return nullptr;

  TString h1 = Form(hpat1, L, chip, ch);
  TH1D* h = dynamic_cast<TH1D*>(d->Get(h1));
  if(h) return h;

  TString h2 = Form("Layer%d_Chip%d_Chn%02d", L, chip, ch);
  h = dynamic_cast<TH1D*>(d->Get(h2));
  return h;
}

inline void zRangeByLayer(int L, double &zmin, double &zmax){
  if (L>=0 && L<=3) { zmin = 15.0; zmax = 27.0; }
  else            { zmin =  6.0; zmax = 15.0; }
}

struct MapPack {
  TH2D* h2[ (cfg::LMAX - cfg::LMIN + 1) ] = {nullptr};
  TString tag; // "SingleHit" / "AllHits"
  TDirectory* baseDir=nullptr; // maps or maps_AllHits
};

MapPack* createGainMaps(TFile* fout, const char* tag){
  auto mp = new MapPack();
  mp->tag = tag;
  TString dname = (tag && TString(tag)=="AllHits") ? "maps_AllHits" : "maps";
  mp->baseDir = ensureDir(fout, dname);

  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    auto* h = new TH2D(Form("gainMap_%s_L%d", tag, L),
                       Form("Layer %d Gain (%s) [ADC/p.e.];chip;channel",L,tag),
                       cfg::NCHIP_AXIS, -0.5, cfg::NCHIP_AXIS-0.5,
                       (cfg::CH1-cfg::CH0+1), cfg::CH0-0.5, cfg::CH1+0.5);
    h->SetStats(false);
    double zmin,zmax; zRangeByLayer(L,zmin,zmax);
    h->SetMinimum(zmin); h->SetMaximum(zmax);
    h->SetDirectory(nullptr);
    mp->h2[L - cfg::LMIN] = h;
  }
  return mp;
}

void writeGainMapsAndCanvases(TFile* fout, MapPack* mp){
  if(!mp) return;
  TDirectory* d = mp->baseDir ? mp->baseDir : ensureDir(fout, (mp->tag=="AllHits") ? "maps_AllHits" : "maps");
  TDirectory::TContext cx(d);

  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    TH2D* h = mp->h2[L - cfg::LMIN];
    if(!h) continue;
    h->Write("", TObject::kOverwrite);

    TCanvas* c = new TCanvas(Form("c_gainMap_%s_L%d", mp->tag.Data(), L),
                             Form("Layer %d Gain Map (%s)",L, mp->tag.Data()), 900, 600);
    c->cd(); h->Draw("COLZ"); c->Write(); delete c;
  }
}

// ---------- Nphoto マップ（2D, レイヤ別） ----------
struct NphotoMaps {
  // tag: "SingleHit" or "AllHits"; cat: "main", "ref", "refPed50"
  std::map<TString, TH2D*> h2_byLayer; // key="nphotoMap_<tag>_<cat>_L<layer>"

  static TString key(int L, const char* tag, const char* cat){
    return TString::Format("nphotoMap_%s_%s_L%d", tag, cat, L);
  }
};

std::unique_ptr<NphotoMaps> createNphotoMaps(TFile* fout, const char* tag, const char* cat){
  auto nm = std::make_unique<NphotoMaps>();
  TDirectory* d = ensureDir(fout, TString::Format("NphotoMaps/%s/%s", tag, cat));
  TDirectory::TContext cx(d);
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    auto *h = new TH2D(NphotoMaps::key(L,tag,cat), TString::Format("Nphoto (%s,%s) L%d;chip;channel",tag,cat,L),
                       cfg::NCHIP_AXIS,-0.5,cfg::NCHIP_AXIS-0.5,
                       (cfg::CH1-cfg::CH0+1), cfg::CH0-0.5, cfg::CH1+0.5);
    h->SetDirectory(d);
    nm->h2_byLayer[h->GetName()] = h;
  }
  return nm;
}

void fillNphotoMap(NphotoMaps* nm, int L, int chip, int ch, double nphoto, const char* tag, const char* cat){
  if(!nm) return;
  TString k = NphotoMaps::key(L,tag,cat);
  auto it = nm->h2_byLayer.find(k);
  if(it==nm->h2_byLayer.end()) return;
  it->second->SetBinContent(chip+1, ch+1, nphoto);
}

// ---------- 位相/可視化補助 ----------
double estimatePhase(const TH1* h,double period){
  if(period<=0) return 0.0;
  int nb = h->GetNbinsX();
  double step = period/400.0;
  double bestPhi=0, bestS=-1e300;
  for(double phi=0; phi<period; phi+=step){
    double s=0;
    for(int b=1;b<=nb;++b){
      double x=h->GetBinCenter(b), y=h->GetBinContent(b);
      s += y*std::cos(2*TMath::Pi()*(x-phi)/period);
    }
    if(s>bestS){ bestS=s; bestPhi=phi; }
  }
  return bestPhi;
}

TF1* makeCos(const TH1* h,double period,double phase,const char* name){
  double x1=h->GetXaxis()->GetXmin(), x2=h->GetXaxis()->GetXmax();
  TF1* f=new TF1(name,"[3] + [0]*cos(2*TMath::Pi()*(x-[2])/[1])",x1,x2);
  double ymax=h->GetMaximum();
  f->SetParameters(0.6*ymax, period, phase, 0.3*ymax);
  f->SetNpx(5000);
  f->SetLineColor(kRed); f->SetLineWidth(2);
  return f;
}

struct ChResult{
  int layer=-1, chip=-1, ch=-1;
  double gain=-1, err=-1, snr=-1;
  int good_main=0;   // 主解析: 統計≥2000 & 異常値除外
  int good_ref =0;   // 参考: 主解析 + (SNR, relerr, temp ok, ped-sharp)
  double pedMean=0, pedSigma=0; int pedOK=0;
};

// 温度（run平均）取得（壊れた Raw_Hit の zlib/TBasket 警告は抑制して試読）
bool getRunMeanTemp(const TString& runFile, double &meanT, int &nevUsed){
  meanT = NAN; nevUsed=0;
  TString fdec = TString(cfg::INDIR_DECODE) + runFile;

  if (gSystem->AccessPathName(fdec)) {
    Warning("getRunTemps","decode file missing: %s", fdec.Data());
    return false;
  }

  std::unique_ptr<TFile> fd(TFile::Open(fdec,"READ"));
  if(!fd || fd->IsZombie()){
    Warning("getRunTemps","cannot open decode file: %s", fdec.Data());
    return false;
  }
  TTree *t = (TTree*)fd->Get("Raw_Hit");
  if(!t){
    Warning("getRunTemps","no TTree Raw_Hit in %s", fdec.Data());
    return false;
  }

  { ErrorSilencer sl(kError);
    Int_t okCnt=0, tryCnt=0;
    Long64_t nTry = std::min<Long64_t>(t->GetEntries(), 50);
    for(Long64_t i=0;i<nTry;++i){ auto br = t->GetEntry(i); tryCnt++; if(br>0) okCnt++; }
    if(tryCnt>0 && okCnt==0){
      Warning("getRunTemps","decode seems corrupted (no entries readable): %s", fdec.Data());
      return false;
    }
  }

  std::vector<std::vector<double>> *Temp=nullptr;
  t->SetBranchAddress("Temperature",&Temp);

  Long64_t n=t->GetEntries();
  const Long64_t MAXE = std::min<Long64_t>(n, 50000);
  long double sum=0; long long cnt=0;

  { ErrorSilencer sl(kError);
    for(Long64_t i=0;i<MAXE;++i){
      if(t->GetEntry(i)<=0) continue;
      if(!Temp) continue;
      for(const auto& arr:*Temp){
        for(double v:arr){ if(std::isfinite(v)){ sum+=v; cnt++; } }
      }
    }
  }

  if(cnt==0){
    Warning("getRunTemps","no finite temperature in %s", fdec.Data());
    return false;
  }
  meanT = (double)(sum/cnt);
  nevUsed = (int)MAXE;
  Info("Temp","%s meanT=%.3fC (events=%d) src=decode/mipScan", runFile.Data(), meanT, nevUsed);
  return true;
}

// =============================================================
// FFTを用いてTH1Dからgainを求める関数
// =============================================================
double GetGainFromFFT(const TH1D* hIn, int layer = 0, double pedMean = 0.0)
{
  if (!hIn) return -1.0;

  // pedestal補正
  const int nb = hIn->GetNbinsX();
  std::vector<double> v(nb, 0.0);
  double mean = 0.0;
  for (int i = 0; i < nb; ++i) {
    double y = hIn->GetBinContent(i + 1);
    v[i] = y;
    mean += y;
  }
  mean /= nb;
  for (int i = 0; i < nb; ++i) v[i] -= mean;  // 平均除去

  // パディング（2の冪長に拡張）
  int Npad = cfg::PAD_TO_POW2 ? nextPow2(nb) : nb;
  if (cfg::PAD_TO_POW2 && Npad > cfg::PAD_MAX) Npad = cfg::PAD_MAX;

  // FFT実行
  std::unique_ptr<TH1D> hFFT(doFFT(v, Npad, "hFFT_temp"));

  // FFTスペクトルからピークを抽出
  KPick pk = pickK0(hFFT.get(), Npad, layer);

  // 結果を返す（異常値なら負の値）
  if (!isFinitePos(pk.gain)) return -1.0;
  return pk.gain;
}


// 1 チャンネル処理（Single/All 共通）
void processChannel_common(TFile* fin,TFile* fout,int L,int chip,int ch,
                           TH1D* hGain, std::vector<ChResult>& out,
                           MapPack* mp, ESource src, RunCounters* rc){
  TH1D* hRaw = getRawHist(fin, L, chip, ch, src);
  if(!hRaw){
    if(rc){ if(src==kSingle) rc->missing_single++; else rc->missing_all++; }
    if(cfg::VERBOSE_MISSING_HIST)
      Warning("process","missing hist L%d chip%d ch%d (%s)",L,chip,ch,(src==kSingle)?"SingleHit":"AllHits");
    return;
  }
  if(hRaw->GetEntries()<cfg::MIN_ENT){
    if(rc){ if(src==kSingle) rc->lowstat_single++; else rc->lowstat_all++; }
    if(cfg::VERBOSE_LOWSTAT)
      Warning("process","low stat hist L%d chip%d ch%d (%s) ent=%lld",
              L,chip,ch,(src==kSingle)?"SingleHit":"AllHits",(Long64_t)hRaw->GetEntries());
    return;
  }

  double mean=cfg::PED_FALLBACK, rms=0;
  gPed.get(L,chip,ch,mean,rms);

  const int nb = hRaw->GetNbinsX();
  const double xmin = hRaw->GetXaxis()->GetXmin() - mean;
  const double xmax = hRaw->GetXaxis()->GetXmax() - mean;

  TH1D* hCorr = new TH1D("tmp_hCorr","ped sub;ADC-ped (ADC);counts", nb, xmin, xmax);
  hCorr->SetDirectory(nullptr);
  for(int b=1;b<=nb;++b){
    double xraw=hRaw->GetBinCenter(b);
    double y   =hRaw->GetBinContent(b);
    hCorr->Fill(xraw-mean,y);
  }
  setHistStyle(hCorr);

  // 平均除去
  std::vector<double> v(nb,0.0); double mu=0;
  for(int i=0;i<nb;++i){ double y=hCorr->GetBinContent(i+1); v[i]=y; mu+=y; }
  if(nb>0) mu/=nb;
  for(int i=0;i<nb;++i) v[i]-=mu;

  // FFT
  int Npad=cfg::PAD_TO_POW2?nextPow2(nb):nb;
  if(cfg::PAD_TO_POW2 && Npad>cfg::PAD_MAX) Npad=cfg::PAD_MAX;
  TH1D* hFFT=doFFT(v,Npad,"tmp_hFFT");
  hFFT->SetDirectory(nullptr);

  // 基本周波数 → ゲイン
  auto pk=pickK0(hFFT,Npad,L);

  // ped-sharp 指標（緩和）
  double pM=0,pS=0; int pOK = pedIsSharp_relaxed(L, chip, ch, pM, pS) ? 1 : 0;

  ChResult cr; cr.layer=L; cr.chip=chip; cr.ch=ch;
  cr.gain=pk.gain; cr.err=pk.err; cr.snr=pk.snr;
  cr.pedMean=pM; cr.pedSigma=pS; cr.pedOK=pOK;

  // 主解析: 異常値除外のみ（統計はここに来る前に通過済み）
  cr.good_main = (isFinitePos(cr.gain) && isFinitePos(cr.err) && std::isfinite(cr.snr)) ? 1 : 0;

  // 参考: 主解析 + (SNR, relerr, ped-sharp)
  bool fitOK = (cr.snr>=cfg::SNR_MIN_REF) && (cr.err>0) && (cr.gain>0) && ((cr.err/cr.gain)<=cfg::RELERR_MAX_REF);
  bool pedOK = (pOK!=0);
  cr.good_ref = (cr.good_main && fitOK && pedOK) ? 1 : 0;

  // 出力保持
  if(pk.gain>0){
    if(hGain){ hGain->SetBinContent(ch+1,pk.gain); hGain->SetBinError(ch+1,pk.err); }
    if(mp && mp->h2[L - cfg::LMIN]) mp->h2[L - cfg::LMIN]->SetBinContent(chip+1, ch+1, pk.gain);
  }
  out.push_back(cr);

  // 個別オーバレイ保存
  TString base = (src==kSingle) ? "layer" : "AllHits/layer";
  TString dpath = Form("%s%d/chip%d/ch%d", base.Data(), L, chip, ch);
  TDirectory* d = ensureDir(fout, dpath);
  {
    TDirectory::TContext cx(d);
    TH1D* hRawOut = (TH1D*)hRaw->Clone((src==kSingle) ? "hRaw" : "hRaw_All");
    hRawOut->SetDirectory(d); hRawOut->Write("", TObject::kOverwrite);

    TH1D* hCorrOut = (TH1D*)hCorr->Clone( (src==kSingle) ? "hCorr" : "hCorr_All" );
    hCorrOut->SetDirectory(d); hCorrOut->Write("", TObject::kOverwrite);

    TH1D* hFFTOut  = (TH1D*)hFFT ->Clone( (src==kSingle) ? "hFFT_full" : "hFFT_full_All" );
    hFFTOut->SetDirectory(d);  hFFTOut->Write("", TObject::kOverwrite);

    TCanvas *c=new TCanvas("c_overlay","overlay",800,600);
    gPad->SetLogy(false);
    hCorrOut->Draw("HIST");
    if(pk.gain>0){
      double phase = estimatePhase(hCorrOut,pk.gain);
      TF1* f=makeCos(hCorrOut,pk.gain,phase,"fWave");
      double xmin=hCorrOut->GetXaxis()->GetXmin(), xmax=hCorrOut->GetXaxis()->GetXmax();
      double ymin=f->GetMinimum(xmin,xmax), ymax=std::max(hCorrOut->GetMaximum(), f->GetMaximum(xmin,xmax));
      hCorrOut->SetMinimum(std::min(0.0,ymin*1.1));
      hCorrOut->SetMaximum(ymax*1.1);
      hCorrOut->Draw("HIST");
      f->Draw("SAME");
    }
    c->Write(); delete c;

    delete hRawOut; delete hCorrOut; delete hFFTOut;
  }
  delete hCorr; delete hFFT;
}

// 1 run 処理
void analyzeOneRun_full(const TString& runFile){
  TString inHist  = TString(cfg::INDIR_HIST)+runFile;
  if (gSystem->AccessPathName(inHist)) {
    Warning("RUN","hist file missing: %s", inHist.Data());
    return;
  }

  TString out     = TString(cfg::OUTDIR)+"/"+runFile;
  out.ReplaceAll(".root","_FFTped_full.root");

  if(!cfg::OVERWRITE && !gSystem->AccessPathName(out)){
    Info("RUN","[skip] %s (exists)", runFile.Data());
    return;
  }

  std::unique_ptr<TFile> fin(TFile::Open(inHist,"READ"));
  if(!fin || fin->IsZombie()){
    Error("RUN","cannot open hist file: %s", inHist.Data());
    return;
  }

  gSystem->mkdir(cfg::OUTDIR,kTRUE);
  std::unique_ptr<TFile> fout(TFile::Open(out,"RECREATE"));
  if(!fout || fout->IsZombie()){
    Error("RUN","cannot create output: %s", out.Data());
    return;
  }
  Info("RUN","%s -> %s", runFile.Data(), out.Data());

  TH1::AddDirectory(kFALSE);
  gStyle->SetEndErrorSize(8);
  gStyle->SetErrorX(0);

  // run 平均温度
  double meanT= NAN; int neT=0;
  bool haveT = getRunMeanTemp(runFile, meanT, neT);
  if(!haveT){
    Warning("Temp","No temperature info for %s -> use scale=1", runFile.Data());
  }

  RunCounters rc{}; // 低統計/欠落の集計

  // ---------- SingleHit ----------
  MapPack* mapsS = createGainMaps(fout.get(), "SingleHit");
  TH1D* hGain_1D_S[ (cfg::LMAX-cfg::LMIN+1) ] = {nullptr};
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    hGain_1D_S[L-cfg::LMIN] = new TH1D(Form("gain_full_L%d",L),"",
                                      cfg::CH1-cfg::CH0+1, cfg::CH0-0.5, cfg::CH1+0.5);
    hGain_1D_S[L-cfg::LMIN]->SetDirectory(nullptr);
  }
  std::vector<ChResult> resS; resS.reserve((cfg::LMAX-cfg::LMIN+1)*cfg::NCHIP_PER_LAYER*(cfg::CH1-cfg::CH0+1));
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L)
    for(int chip=cfg::CHIP_MIN; chip<=cfg::CHIP_MAX; ++chip)
      for(int ch=cfg::CH0; ch<=cfg::CH1; ++ch)
        processChannel_common(fin.get(), fout.get(), L, chip, ch,
                              hGain_1D_S[L-cfg::LMIN], resS, mapsS, kSingle, &rc);

  if(cfg::WRITE_GAIN_HIST){
    TDirectory::TContext cx(fout.get());
    for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){ setHistStyle(hGain_1D_S[L-cfg::LMIN]); hGain_1D_S[L-cfg::LMIN]->Write(); }
  }
  writeGainMapsAndCanvases(fout.get(), mapsS);
  delete mapsS;

  // per-channel 結果ツリー（SingleHit）
  {
    TDirectory::TContext cx(fout.get());
    int layer, chip, ch; double gain, err, snr; int good_main, good_ref; double pedMean,pedSigma; int pedOK;
    TTree* tRes=new TTree("GoodChannels","SPE pick result (SingleHit) with main/ref flags");
    tRes->SetDirectory(fout.get());
    tRes->Branch("layer",&layer,"layer/I");
    tRes->Branch("chip",&chip,"chip/I");
    tRes->Branch("ch",&ch,"ch/I");
    tRes->Branch("gain",&gain,"gain/D");
    tRes->Branch("err",&err,"err/D");
    tRes->Branch("snr",&snr,"snr/D");
    tRes->Branch("good_main",&good_main,"good_main/I");
    tRes->Branch("good_ref",&good_ref,"good_ref/I");
    tRes->Branch("pedMean",&pedMean,"pedMean/D");
    tRes->Branch("pedSigma",&pedSigma,"pedSigma/D");
    tRes->Branch("pedOK",&pedOK,"pedOK/I");
    for(const auto&r:resS){
      layer=r.layer; chip=r.chip; ch=r.ch; gain=r.gain; err=r.err; snr=r.snr;
      good_main=r.good_main; good_ref=r.good_ref; pedMean=r.pedMean; pedSigma=r.pedSigma; pedOK=r.pedOK;
      tRes->Fill();
    }
    tRes->Write();
  }
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L) delete hGain_1D_S[L-cfg::LMIN];

  // ---------- AllHits ----------
  MapPack* mapsA = createGainMaps(fout.get(), "AllHits");
  TH1D* hGain_1D_A[ (cfg::LMAX-cfg::LMIN+1) ] = {nullptr};
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    hGain_1D_A[L-cfg::LMIN] = new TH1D(Form("gain_full_All_L%d",L),"",
                                      cfg::CH1-cfg::CH0+1, cfg::CH0-0.5, cfg::CH1+0.5);
    hGain_1D_A[L-cfg::LMIN]->SetDirectory(nullptr);
  }
  std::vector<ChResult> resA; resA.reserve((cfg::LMAX-cfg::LMIN+1)*cfg::NCHIP_PER_LAYER*(cfg::CH1-cfg::CH0+1));
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L)
    for(int chip=cfg::CHIP_MIN; chip<=cfg::CHIP_MAX; ++chip)
      for(int ch=cfg::CH0; ch<=cfg::CH1; ++ch)
        processChannel_common(fin.get(), fout.get(), L, chip, ch,
                              hGain_1D_A[L-cfg::LMIN], resA, mapsA, kAll, &rc);

  if(cfg::WRITE_GAIN_HIST){
    TDirectory::TContext cx(fout.get());
    for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){ setHistStyle(hGain_1D_A[L-cfg::LMIN]); hGain_1D_A[L-cfg::LMIN]->Write(); }
  }
  writeGainMapsAndCanvases(fout.get(), mapsA);
  delete mapsA;

  // per-channel 結果ツリー（AllHits）
  {
    TDirectory::TContext cx(fout.get());
    int layer, chip, ch; double gain, err, snr; int good_main, good_ref; double pedMean,pedSigma; int pedOK;
    TTree* tResA=new TTree("GoodChannels_AllHits","SPE pick result (AllHits) with main/ref flags");
    tResA->SetDirectory(fout.get());
    tResA->Branch("layer",&layer,"layer/I");
    tResA->Branch("chip",&chip,"chip/I");
    tResA->Branch("ch",&ch,"ch/I");
    tResA->Branch("gain",&gain,"gain/D");
    tResA->Branch("err",&err,"err/D");
    tResA->Branch("snr",&snr,"snr/D");
    tResA->Branch("good_main",&good_main,"good_main/I");
    tResA->Branch("good_ref",&good_ref,"good_ref/I");
    tResA->Branch("pedMean",&pedMean,"pedMean/D");
    tResA->Branch("pedSigma",&pedSigma,"pedSigma/D");
    tResA->Branch("pedOK",&pedOK,"pedOK/I");
    for(const auto&r:resA){
      layer=r.layer; chip=r.chip; ch=r.ch; gain=r.gain; err=r.err; snr=r.snr;
      good_main=r.good_main; good_ref=r.good_ref; pedMean=r.pedMean; pedSigma=r.pedSigma; pedOK=r.pedOK;
      tResA->Fill();
    }
    tResA->Write();
  }
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L) delete hGain_1D_A[L-cfg::LMIN];

  // ---------- 温度情報 ----------
  {
    TDirectory::TContext cx(fout.get());
    TDirectory* d=ensureDir(fout.get(),"runInfo");
    TDirectory::TContext cy(d);
    TParameter<double>("RunMeanTempC", haveT?meanT:NAN).Write("RunMeanTempC");
    TParameter<int>("RunTempNevents", neT).Write("RunTempNevents");
  }

  // ---------- Nphoto 計算（SingleHit MPV を使用） ----------
  // MPV 読み出し（key = L*100000 + chip*1000 + ch ではなく、ここでは map のキーは (L,chip,ch) タプル化）
  struct Key { int L, C, ch; bool operator<(const Key& o) const {
    if(L!=o.L) return L<o.L; if(C!=o.C) return C<o.C; return ch<o.ch; } };
  std::map<Key,double> mpvMap;
  std::map<Key,int>    mpvEntries;

  auto loadMPV = [&](const TString& runName)->void{
    TString mipDir  = Form("%s/%s", cfg::MIPFIT_BASE, runName.Data());
    TString mipFile = Form("%s/mip_%s.root", mipDir.Data(), runName.Data());
    if (gSystem->AccessPathName(mipFile)) {
      TString alt = Form("%s/%s.root", mipDir.Data(), runName.Data());
      if (!gSystem->AccessPathName(alt)) mipFile = alt;
    }

    if (gSystem->AccessPathName(mipFile)) {
      Warning("MIP","cannot find MIP file: %s (dir=%s)", mipFile.Data(), mipDir.Data());
      return;
    }

    std::unique_ptr<TFile> fm(TFile::Open(mipFile, "READ"));
    if (!fm || fm->IsZombie()) {
      Warning("MIP","cannot open MIP file: %s", mipFile.Data());
      return;
    }

    // 新: MIP_FIT（SingleHitのみ）を優先
    TTree* tNew = (TTree*)fm->Get("MIP_FIT");
    if (tNew){
      Int_t    ChnEntries=0, CellID=0, NDF=0;
      Double_t LandauMPV=0, ChiSquare=0;
      tNew->SetBranchAddress("ChnEntries",&ChnEntries);
      tNew->SetBranchAddress("CellID",&CellID);
      tNew->SetBranchAddress("LandauMPV",&LandauMPV);
      if (tNew->GetBranch("ChiSquare")) tNew->SetBranchAddress("ChiSquare",&ChiSquare);
      if (tNew->GetBranch("NDF"))       tNew->SetBranchAddress("NDF",&NDF);

      Long64_t n=tNew->GetEntries(), kept=0, skipStat=0, skipBad=0, skipRChi2=0, skipOOB=0;
      for(Long64_t i=0;i<n;++i){
        tNew->GetEntry(i);
        // decode CellID = L*100000 + C*10000 + ch
        int L = CellID/100000;
        int C = (CellID%100000)/10000;
        int ch= (CellID%10000);
        if(L<cfg::LMIN || L>cfg::LMAX || C<0 || C>=cfg::NCHIP_PER_LAYER || ch<cfg::CH0 || ch>cfg::CH1){ ++skipOOB; continue; }
        if(ChnEntries < cfg::MIP_MIN_ENT){ ++skipStat; continue; }
        if(!(LandauMPV>0) || !std::isfinite(LandauMPV) || LandauMPV>1e6){ ++skipBad; continue; }
        double rchi2 = (NDF>0)? (ChiSquare/NDF) : 1e9;
        if(!(rchi2 >= 0.0) || rchi2 > cfg::MPV_RCHI2_MAX){ ++skipRChi2; continue; }

        Key k{L,C,ch};
        mpvMap[k]     = LandauMPV;
        mpvEntries[k] = ChnEntries;
        kept++;
      }
      Info("MIP","MPV(SingleHit) loaded from %s : rows=%lld, kept=%lld, skipped=[oob=%lld, stat<%d=%lld, bad=%lld, rchi2>%g=%lld]",
           mipFile.Data(), n, kept, skipOOB, cfg::MIP_MIN_ENT, skipStat, skipBad, cfg::MPV_RCHI2_MAX, skipRChi2);
      return;
    }

    // 旧: MIP_Fit（互換）
    TTree* tOld = (TTree*)fm->Get("MIP_Fit");
    if (tOld){
      Int_t Layer=0, Chip=0, Chn=0, Entries=0, NDF=0;
      Double_t MPV=0, Chi2=0;
      tOld->SetBranchAddress("Layer",&Layer);
      tOld->SetBranchAddress("Chip",&Chip);
      tOld->SetBranchAddress("Chn",&Chn);
      tOld->SetBranchAddress("Entries",&Entries);
      tOld->SetBranchAddress("MPV",&MPV);
      if (tOld->GetBranch("Chi2")) tOld->SetBranchAddress("Chi2",&Chi2);
      if (tOld->GetBranch("NDF"))  tOld->SetBranchAddress("NDF",&NDF);

      Long64_t n=tOld->GetEntries(), kept=0, skipStat=0, skipBad=0, skipRChi2=0, skipOOB=0;
      for(Long64_t i=0;i<n;++i){
        tOld->GetEntry(i);
        if(Layer<cfg::LMIN || Layer>cfg::LMAX || Chip<0 || Chip>=cfg::NCHIP_PER_LAYER || Chn<cfg::CH0 || Chn>cfg::CH1){ ++skipOOB; continue; }
        if(Entries < cfg::MIP_MIN_ENT){ ++skipStat; continue; }
        if(!(MPV>0) || !std::isfinite(MPV) || MPV>1e6){ ++skipBad; continue; }
        double rchi2 = (NDF>0)? (Chi2/NDF) : 1e9;
        if(!(rchi2 >= 0.0) || rchi2 > cfg::MPV_RCHI2_MAX){ ++skipRChi2; continue; }

        Key k{Layer,Chip,Chn};
        mpvMap[k]     = MPV;
        mpvEntries[k] = Entries;
        kept++;
      }
      Info("MIP","MPV(compat:old) loaded from %s : rows=%lld, kept=%lld, skipped=[oob=%lld, stat<%d=%lld, bad=%lld, rchi2>%g=%lld]",
           mipFile.Data(), n, kept, skipOOB, cfg::MIP_MIN_ENT, skipStat, skipBad, cfg::MPV_RCHI2_MAX, skipRChi2);
      return;
    }

    Warning("MIP","tree missing: MIP_FIT & MIP_Fit (file=%s)", mipFile.Data());
  };

  // runName の抽出と MPV ロード
  TString runName = gSystem->BaseName(runFile);
  if (runName.EndsWith(".root")) runName.ReplaceAll(".root","");
  mpvMap.clear(); mpvEntries.clear();
  loadMPV(runName);

  auto getMPV = [&](int L,int chip,int ch,double &mpv,int &ent)->bool{
    Key k{L,chip,ch};
    auto it = mpvMap.find(k);
    if(it==mpvMap.end()) return false;
    mpv = it->second;
    ent = mpvEntries[k];
    return true;
  };

  // ---------- Nphoto 出力（ツリー & マップ） ----------
  auto writeNphoto = [&](const char* tagSPE, const std::vector<ChResult>& rr, const char* treeNameMain, const char* treeNameRef, const char* treeNameRefPed50){
    TDirectory::TContext cx(fout.get());
    TDirectory* d = ensureDir(fout.get(), TString::Format("Nphoto/%s", tagSPE));
    TDirectory::TContext cy(d);

    // ツリー
    Int_t L,chip,ch; Double_t mpv,gain,nphoto; Int_t mpvEnt;
    Double_t pedSigma;
    TTree* tMain = new TTree(treeNameMain, "Nphoto main (entries>=2000 & finite)");
    tMain->Branch("L",&L,"L/I"); tMain->Branch("chip",&chip,"chip/I"); tMain->Branch("ch",&ch,"ch/I");
    tMain->Branch("mpv",&mpv,"mpv/D"); tMain->Branch("gain",&gain,"gain/D");
    tMain->Branch("nphoto",&nphoto,"nphoto/D"); tMain->Branch("mpvEntries",&mpvEnt,"mpvEntries/I");
    tMain->Branch("pedSigma",&pedSigma,"pedSigma/D");

    TTree* tRef  = new TTree(treeNameRef , "Nphoto reference (main + fit & ped-sharp)");
    tRef ->Branch("L",&L,"L/I"); tRef ->Branch("chip",&chip,"chip/I"); tRef ->Branch("ch",&ch,"ch/I");
    tRef ->Branch("mpv",&mpv,"mpv/D"); tRef ->Branch("gain",&gain,"gain/D");
    tRef ->Branch("nphoto",&nphoto,"nphoto/D"); tRef ->Branch("mpvEntries",&mpvEnt,"mpvEntries/I");
    tRef ->Branch("pedSigma",&pedSigma,"pedSigma/D");

    TTree* tRef50  = new TTree(treeNameRefPed50 , "Nphoto reference (ref & pedSigma<50)");
    tRef50 ->Branch("L",&L,"L/I"); tRef50 ->Branch("chip",&chip,"chip/I"); tRef50 ->Branch("ch",&ch,"ch/I");
    tRef50 ->Branch("mpv",&mpv,"mpv/D"); tRef50 ->Branch("gain",&gain,"gain/D");
    tRef50 ->Branch("nphoto",&nphoto,"nphoto/D"); tRef50 ->Branch("mpvEntries",&mpvEnt,"mpvEntries/I");
    tRef50 ->Branch("pedSigma",&pedSigma,"pedSigma/D");

    // マップ（2D, レイヤ別）
    auto mapMain  = createNphotoMaps(fout.get(), tagSPE, "main");
    auto mapRef   = createNphotoMaps(fout.get(), tagSPE, "ref");
    auto mapRef50 = createNphotoMaps(fout.get(), tagSPE, "refPed50");

    int filledMain=0, filledRef=0, filledRef50=0;

    for(const auto &r: rr){
      L=r.layer; chip=r.chip; ch=r.ch;
      double m; int entM=0;
      if(!getMPV(L,chip,ch,m,entM)) continue;
      if(entM < cfg::MIP_MIN_ENT) continue;

      if(!isFinitePos(r.gain) || !isFinitePos(m)) continue;

      mpv=m; mpvEnt=entM; gain=r.gain; nphoto = mpv/gain;
      pedSigma = r.pedSigma;

      if(isFinitePos(nphoto)){
        tMain->Fill(); filledMain++;
        fillNphotoMap(mapMain.get(), L, chip, ch, nphoto, tagSPE, "main");
      }

      bool fitOK = (r.snr>=cfg::SNR_MIN_REF) && (r.err>0) && (r.gain>0) && ((r.err/r.gain)<=cfg::RELERR_MAX_REF);
      bool pedOK = (r.pedOK!=0);
      if(fitOK && pedOK && isFinitePos(nphoto)){
        tRef->Fill(); filledRef++;
        fillNphotoMap(mapRef.get(), L, chip, ch, nphoto, tagSPE, "ref");
        if (pedSigma < 50.0){
          tRef50->Fill(); filledRef50++;
          fillNphotoMap(mapRef50.get(), L, chip, ch, nphoto, tagSPE, "refPed50");
        }
      }
    }

    tMain->Write("", TObject::kOverwrite);
    tRef ->Write("", TObject::kOverwrite);
    tRef50->Write("", TObject::kOverwrite);

    Info("Nphoto","%s: main=%d, ref=%d, refPed50=%d", tagSPE, filledMain, filledRef, filledRef50);
  };

  writeNphoto("SingleHit", resS, "Nphoto_main_SingleHit", "Nphoto_ref_SingleHit", "Nphoto_refPed50_SingleHit");
  writeNphoto("AllHits"  , resA, "Nphoto_main_AllHits"  , "Nphoto_ref_AllHits"  , "Nphoto_refPed50_AllHits"  );

  // 低統計/欠落のサマリのみ表示
  Info("summary","Run=%s: missing(Single/All)=%ld/%ld, lowstat(Single/All)=%ld/%ld",
       runFile.Data(), rc.missing_single, rc.missing_all, rc.lowstat_single, rc.lowstat_all);

  Info("RUN","Done run: %s", runFile.Data());
}

// ==================== メイン（引数あり/なし両対応） ====================
// 引数: only_run_path="" → 全run。非空ならそのファイル名の run だけ実行
int fft_fullrange_tempcorr_plusNphoto_v5(const char* only_run_path /*= ""*/){
  gROOT->SetBatch(kTRUE);

  // pedestal（平均用）
  if(!gPed.load(cfg::PEDFILE)){
    Warning("Ped","failed to load pedestal: %s (will fallback %.1f)", cfg::PEDFILE, cfg::PED_FALLBACK);
  }

  // 実行対象の run リスト
  std::vector<TString> runFiles;

  if(only_run_path && std::strlen(only_run_path)>0){
    TString base = gSystem->BaseName(only_run_path);
    if(!base.EndsWith(".root")){
      Error("MAIN","ONLY-RUN must be a .root file: %s", only_run_path);
      return 1;
    }
    TString inHist = TString(cfg::INDIR_HIST)+base;
    if (gSystem->AccessPathName(inHist)) {
      Error("MAIN","hist file not found: %s", inHist.Data());
      return 2;
    }
    runFiles.push_back(base);
    Info("MAIN","[single-run mode] %s", base.Data());
  }else{
    TSystemDirectory dir("indir", cfg::INDIR_HIST);
    std::unique_ptr<TList> list(dir.GetListOfFiles());
    if(!list){
      Error("MAIN","cannot list %s", cfg::INDIR_HIST);
      return 3;
    }
    TIter next(list.get()); TSystemFile* fobj=nullptr;
    while((fobj=(TSystemFile*)next())){
      TString name=fobj->GetName();
      if(fobj->IsDirectory()) continue;
      if(!name.EndsWith(".root")) continue;
      if(name.Contains("_FFTped_full.root")) continue; // 出力はスキップ
      runFiles.push_back(name);
    }
  }

  // 実行
  int nRun=0;
  for(const auto& rf: runFiles){
    analyzeOneRun_full(rf);
    nRun++;
  }
  Info("MAIN","Done. processed=%d", nRun);
  return 0;
}

// 互換：引数なし版（全run）
void fft_fullrange_tempcorr_plusNphoto_v5(){
  fft_fullrange_tempcorr_plusNphoto_v5("");
}

// ショートカット：Run171だけ
int fft_v5_run171_only(){
  return fft_fullrange_tempcorr_plusNphoto_v5(
    "/megraid01/users/data_beamtest/disk-BeamTestData1/Analyzed/decode/mipScan/ECAL_Run171_20250605_162329.root"
  );
}





void analyzeOneRun_full_(const TString& runFile){
  TString inHist  = TString(cfg::INDIR_HIST)+runFile;
  if (gSystem->AccessPathName(inHist)) {
    Warning("RUN","hist file missing: %s", inHist.Data());
    return;
  }

  TString out     = TString(cfg::OUTDIR)+"/"+runFile;
  out.ReplaceAll(".root","_FFTped_full.root");

  if(!cfg::OVERWRITE && !gSystem->AccessPathName(out)){
    Info("RUN","[skip] %s (exists)", runFile.Data());
    return;
  }

  std::unique_ptr<TFile> fin(TFile::Open(inHist,"READ"));
  if(!fin || fin->IsZombie()){
    Error("RUN","cannot open hist file: %s", inHist.Data());
    return;
  }

  gSystem->mkdir(cfg::OUTDIR,kTRUE);
  std::unique_ptr<TFile> fout(TFile::Open(out,"RECREATE"));
  if(!fout || fout->IsZombie()){
    Error("RUN","cannot create output: %s", out.Data());
    return;
  }
  Info("RUN","%s -> %s", runFile.Data(), out.Data());

  TH1::AddDirectory(kFALSE);
  gStyle->SetEndErrorSize(8);
  gStyle->SetErrorX(0);

  // run 平均温度
  double meanT= NAN; int neT=0;
  bool haveT = getRunMeanTemp(runFile, meanT, neT);
  if(!haveT){
    Warning("Temp","No temperature info for %s -> use scale=1", runFile.Data());
  }

  RunCounters rc{}; // 低統計/欠落の集計

  // ---------- SingleHit ----------
  MapPack* mapsS = createGainMaps(fout.get(), "SingleHit");
  TH1D* hGain_1D_S[ (cfg::LMAX-cfg::LMIN+1) ] = {nullptr};
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    hGain_1D_S[L-cfg::LMIN] = new TH1D(Form("gain_full_L%d",L),"",
                                      cfg::CH1-cfg::CH0+1, cfg::CH0-0.5, cfg::CH1+0.5);
    hGain_1D_S[L-cfg::LMIN]->SetDirectory(nullptr);
  }
  std::vector<ChResult> resS; resS.reserve((cfg::LMAX-cfg::LMIN+1)*cfg::NCHIP_PER_LAYER*(cfg::CH1-cfg::CH0+1));
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L)
    for(int chip=cfg::CHIP_MIN; chip<=cfg::CHIP_MAX; ++chip)
      for(int ch=cfg::CH0; ch<=cfg::CH1; ++ch)
        processChannel_common(fin.get(), fout.get(), L, chip, ch,
                              hGain_1D_S[L-cfg::LMIN], resS, mapsS, kSingle, &rc);

  if(cfg::WRITE_GAIN_HIST){
    TDirectory::TContext cx(fout.get());
    for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){ setHistStyle(hGain_1D_S[L-cfg::LMIN]); hGain_1D_S[L-cfg::LMIN]->Write(); }
  }
  writeGainMapsAndCanvases(fout.get(), mapsS);
  delete mapsS;

  // per-channel 結果ツリー（SingleHit）
  {
    TDirectory::TContext cx(fout.get());
    int layer, chip, ch; double gain, err, snr; int good_main, good_ref; double pedMean,pedSigma; int pedOK;
    TTree* tRes=new TTree("GoodChannels","SPE pick result (SingleHit) with main/ref flags");
    tRes->SetDirectory(fout.get());
    tRes->Branch("layer",&layer,"layer/I");
    tRes->Branch("chip",&chip,"chip/I");
    tRes->Branch("ch",&ch,"ch/I");
    tRes->Branch("gain",&gain,"gain/D");
    tRes->Branch("err",&err,"err/D");
    tRes->Branch("snr",&snr,"snr/D");
    tRes->Branch("good_main",&good_main,"good_main/I");
    tRes->Branch("good_ref",&good_ref,"good_ref/I");
    tRes->Branch("pedMean",&pedMean,"pedMean/D");
    tRes->Branch("pedSigma",&pedSigma,"pedSigma/D");
    tRes->Branch("pedOK",&pedOK,"pedOK/I");
    for(const auto&r:resS){
      layer=r.layer; chip=r.chip; ch=r.ch; gain=r.gain; err=r.err; snr=r.snr;
      good_main=r.good_main; good_ref=r.good_ref; pedMean=r.pedMean; pedSigma=r.pedSigma; pedOK=r.pedOK;
      tRes->Fill();
    }
    tRes->Write();
  }
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L) delete hGain_1D_S[L-cfg::LMIN];

  // ---------- AllHits ----------
  MapPack* mapsA = createGainMaps(fout.get(), "AllHits");
  TH1D* hGain_1D_A[ (cfg::LMAX-cfg::LMIN+1) ] = {nullptr};
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){
    hGain_1D_A[L-cfg::LMIN] = new TH1D(Form("gain_full_All_L%d",L),"",
                                      cfg::CH1-cfg::CH0+1, cfg::CH0-0.5, cfg::CH1+0.5);
    hGain_1D_A[L-cfg::LMIN]->SetDirectory(nullptr);
  }
  std::vector<ChResult> resA; resA.reserve((cfg::LMAX-cfg::LMIN+1)*cfg::NCHIP_PER_LAYER*(cfg::CH1-cfg::CH0+1));
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L)
    for(int chip=cfg::CHIP_MIN; chip<=cfg::CHIP_MAX; ++chip)
      for(int ch=cfg::CH0; ch<=cfg::CH1; ++ch)
        processChannel_common(fin.get(), fout.get(), L, chip, ch,
                              hGain_1D_A[L-cfg::LMIN], resA, mapsA, kAll, &rc);

  if(cfg::WRITE_GAIN_HIST){
    TDirectory::TContext cx(fout.get());
    for(int L=cfg::LMIN; L<=cfg::LMAX; ++L){ setHistStyle(hGain_1D_A[L-cfg::LMIN]); hGain_1D_A[L-cfg::LMIN]->Write(); }
  }
  writeGainMapsAndCanvases(fout.get(), mapsA);
  delete mapsA;

  // per-channel 結果ツリー（AllHits）
  {
    TDirectory::TContext cx(fout.get());
    int layer, chip, ch; double gain, err, snr; int good_main, good_ref; double pedMean,pedSigma; int pedOK;
    TTree* tResA=new TTree("GoodChannels_AllHits","SPE pick result (AllHits) with main/ref flags");
    tResA->SetDirectory(fout.get());
    tResA->Branch("layer",&layer,"layer/I");
    tResA->Branch("chip",&chip,"chip/I");
    tResA->Branch("ch",&ch,"ch/I");
    tResA->Branch("gain",&gain,"gain/D");
    tResA->Branch("err",&err,"err/D");
    tResA->Branch("snr",&snr,"snr/D");
    tResA->Branch("good_main",&good_main,"good_main/I");
    tResA->Branch("good_ref",&good_ref,"good_ref/I");
    tResA->Branch("pedMean",&pedMean,"pedMean/D");
    tResA->Branch("pedSigma",&pedSigma,"pedSigma/D");
    tResA->Branch("pedOK",&pedOK,"pedOK/I");
    for(const auto&r:resA){
      layer=r.layer; chip=r.chip; ch=r.ch; gain=r.gain; err=r.err; snr=r.snr;
      good_main=r.good_main; good_ref=r.good_ref; pedMean=r.pedMean; pedSigma=r.pedSigma; pedOK=r.pedOK;
      tResA->Fill();
    }
    tResA->Write();
  }
  for(int L=cfg::LMIN; L<=cfg::LMAX; ++L) delete hGain_1D_A[L-cfg::LMIN];

  // ---------- 温度情報 ----------
  {
    TDirectory::TContext cx(fout.get());
    TDirectory* d=ensureDir(fout.get(),"runInfo");
    TDirectory::TContext cy(d);
    TParameter<double>("RunMeanTempC", haveT?meanT:NAN).Write("RunMeanTempC");
    TParameter<int>("RunTempNevents", neT).Write("RunTempNevents");
  }

  // ---------- Nphoto 計算（SingleHit MPV を使用） ----------
  // MPV 読み出し（key = L*100000 + chip*1000 + ch ではなく、ここでは map のキーは (L,chip,ch) タプル化）
  struct Key { int L, C, ch; bool operator<(const Key& o) const {
    if(L!=o.L) return L<o.L; if(C!=o.C) return C<o.C; return ch<o.ch; } };
  std::map<Key,double> mpvMap;
  std::map<Key,int>    mpvEntries;

  auto loadMPV = [&](const TString& runName)->void{
    TString mipDir  = Form("%s/%s", cfg::MIPFIT_BASE, runName.Data());
    TString mipFile = Form("%s/mip_%s.root", mipDir.Data(), runName.Data());
    if (gSystem->AccessPathName(mipFile)) {
      TString alt = Form("%s/%s.root", mipDir.Data(), runName.Data());
      if (!gSystem->AccessPathName(alt)) mipFile = alt;
    }

    if (gSystem->AccessPathName(mipFile)) {
      Warning("MIP","cannot find MIP file: %s (dir=%s)", mipFile.Data(), mipDir.Data());
      return;
    }

    std::unique_ptr<TFile> fm(TFile::Open(mipFile, "READ"));
    if (!fm || fm->IsZombie()) {
      Warning("MIP","cannot open MIP file: %s", mipFile.Data());
      return;
    }

    // 新: MIP_FIT（SingleHitのみ）を優先
    TTree* tNew = (TTree*)fm->Get("MIP_FIT");
    if (tNew){
      Int_t    ChnEntries=0, CellID=0, NDF=0;
      Double_t LandauMPV=0, ChiSquare=0;
      tNew->SetBranchAddress("ChnEntries",&ChnEntries);
      tNew->SetBranchAddress("CellID",&CellID);
      tNew->SetBranchAddress("LandauMPV",&LandauMPV);
      if (tNew->GetBranch("ChiSquare")) tNew->SetBranchAddress("ChiSquare",&ChiSquare);
      if (tNew->GetBranch("NDF"))       tNew->SetBranchAddress("NDF",&NDF);

      Long64_t n=tNew->GetEntries(), kept=0, skipStat=0, skipBad=0, skipRChi2=0, skipOOB=0;
      for(Long64_t i=0;i<n;++i){
        tNew->GetEntry(i);
        // decode CellID = L*100000 + C*10000 + ch
        int L = CellID/100000;
        int C = (CellID%100000)/10000;
        int ch= (CellID%10000);
        if(L<cfg::LMIN || L>cfg::LMAX || C<0 || C>=cfg::NCHIP_PER_LAYER || ch<cfg::CH0 || ch>cfg::CH1){ ++skipOOB; continue; }
        if(ChnEntries < cfg::MIP_MIN_ENT){ ++skipStat; continue; }
        if(!(LandauMPV>0) || !std::isfinite(LandauMPV) || LandauMPV>1e6){ ++skipBad; continue; }
        double rchi2 = (NDF>0)? (ChiSquare/NDF) : 1e9;
        if(!(rchi2 >= 0.0) || rchi2 > cfg::MPV_RCHI2_MAX){ ++skipRChi2; continue; }

        Key k{L,C,ch};
        mpvMap[k]     = LandauMPV;
        mpvEntries[k] = ChnEntries;
        kept++;
      }
      Info("MIP","MPV(SingleHit) loaded from %s : rows=%lld, kept=%lld, skipped=[oob=%lld, stat<%d=%lld, bad=%lld, rchi2>%g=%lld]",
           mipFile.Data(), n, kept, skipOOB, cfg::MIP_MIN_ENT, skipStat, skipBad, cfg::MPV_RCHI2_MAX, skipRChi2);
      return;
    }

    // 旧: MIP_Fit（互換）
    TTree* tOld = (TTree*)fm->Get("MIP_Fit");
    if (tOld){
      Int_t Layer=0, Chip=0, Chn=0, Entries=0, NDF=0;
      Double_t MPV=0, Chi2=0;
      tOld->SetBranchAddress("Layer",&Layer);
      tOld->SetBranchAddress("Chip",&Chip);
      tOld->SetBranchAddress("Chn",&Chn);
      tOld->SetBranchAddress("Entries",&Entries);
      tOld->SetBranchAddress("MPV",&MPV);
      if (tOld->GetBranch("Chi2")) tOld->SetBranchAddress("Chi2",&Chi2);
      if (tOld->GetBranch("NDF"))  tOld->SetBranchAddress("NDF",&NDF);

      Long64_t n=tOld->GetEntries(), kept=0, skipStat=0, skipBad=0, skipRChi2=0, skipOOB=0;
      for(Long64_t i=0;i<n;++i){
        tOld->GetEntry(i);
        if(Layer<cfg::LMIN || Layer>cfg::LMAX || Chip<0 || Chip>=cfg::NCHIP_PER_LAYER || Chn<cfg::CH0 || Chn>cfg::CH1){ ++skipOOB; continue; }
        if(Entries < cfg::MIP_MIN_ENT){ ++skipStat; continue; }
        if(!(MPV>0) || !std::isfinite(MPV) || MPV>1e6){ ++skipBad; continue; }
        double rchi2 = (NDF>0)? (Chi2/NDF) : 1e9;
        if(!(rchi2 >= 0.0) || rchi2 > cfg::MPV_RCHI2_MAX){ ++skipRChi2; continue; }

        Key k{Layer,Chip,Chn};
        mpvMap[k]     = MPV;
        mpvEntries[k] = Entries;
        kept++;
      }
      Info("MIP","MPV(compat:old) loaded from %s : rows=%lld, kept=%lld, skipped=[oob=%lld, stat<%d=%lld, bad=%lld, rchi2>%g=%lld]",
           mipFile.Data(), n, kept, skipOOB, cfg::MIP_MIN_ENT, skipStat, skipBad, cfg::MPV_RCHI2_MAX, skipRChi2);
      return;
    }

    Warning("MIP","tree missing: MIP_FIT & MIP_Fit (file=%s)", mipFile.Data());
  };

  // runName の抽出と MPV ロード
  TString runName = gSystem->BaseName(runFile);
  if (runName.EndsWith(".root")) runName.ReplaceAll(".root","");
  mpvMap.clear(); mpvEntries.clear();
  loadMPV(runName);

  auto getMPV = [&](int L,int chip,int ch,double &mpv,int &ent)->bool{
    Key k{L,chip,ch};
    auto it = mpvMap.find(k);
    if(it==mpvMap.end()) return false;
    mpv = it->second;
    ent = mpvEntries[k];
    return true;
  };

  // ---------- Nphoto 出力（ツリー & マップ） ----------
  auto writeNphoto = [&](const char* tagSPE, const std::vector<ChResult>& rr, const char* treeNameMain, const char* treeNameRef, const char* treeNameRefPed50){
    TDirectory::TContext cx(fout.get());
    TDirectory* d = ensureDir(fout.get(), TString::Format("Nphoto/%s", tagSPE));
    TDirectory::TContext cy(d);

    // ツリー
    Int_t L,chip,ch; Double_t mpv,gain,nphoto; Int_t mpvEnt;
    Double_t pedSigma;
    TTree* tMain = new TTree(treeNameMain, "Nphoto main (entries>=2000 & finite)");
    tMain->Branch("L",&L,"L/I"); tMain->Branch("chip",&chip,"chip/I"); tMain->Branch("ch",&ch,"ch/I");
    tMain->Branch("mpv",&mpv,"mpv/D"); tMain->Branch("gain",&gain,"gain/D");
    tMain->Branch("nphoto",&nphoto,"nphoto/D"); tMain->Branch("mpvEntries",&mpvEnt,"mpvEntries/I");
    tMain->Branch("pedSigma",&pedSigma,"pedSigma/D");

    TTree* tRef  = new TTree(treeNameRef , "Nphoto reference (main + fit & ped-sharp)");
    tRef ->Branch("L",&L,"L/I"); tRef ->Branch("chip",&chip,"chip/I"); tRef ->Branch("ch",&ch,"ch/I");
    tRef ->Branch("mpv",&mpv,"mpv/D"); tRef ->Branch("gain",&gain,"gain/D");
    tRef ->Branch("nphoto",&nphoto,"nphoto/D"); tRef ->Branch("mpvEntries",&mpvEnt,"mpvEntries/I");
    tRef ->Branch("pedSigma",&pedSigma,"pedSigma/D");

    TTree* tRef50  = new TTree(treeNameRefPed50 , "Nphoto reference (ref & pedSigma<50)");
    tRef50 ->Branch("L",&L,"L/I"); tRef50 ->Branch("chip",&chip,"chip/I"); tRef50 ->Branch("ch",&ch,"ch/I");
    tRef50 ->Branch("mpv",&mpv,"mpv/D"); tRef50 ->Branch("gain",&gain,"gain/D");
    tRef50 ->Branch("nphoto",&nphoto,"nphoto/D"); tRef50 ->Branch("mpvEntries",&mpvEnt,"mpvEntries/I");
    tRef50 ->Branch("pedSigma",&pedSigma,"pedSigma/D");

    // マップ（2D, レイヤ別）
    auto mapMain  = createNphotoMaps(fout.get(), tagSPE, "main");
    auto mapRef   = createNphotoMaps(fout.get(), tagSPE, "ref");
    auto mapRef50 = createNphotoMaps(fout.get(), tagSPE, "refPed50");

    int filledMain=0, filledRef=0, filledRef50=0;

    for(const auto &r: rr){
      L=r.layer; chip=r.chip; ch=r.ch;
      double m; int entM=0;
      if(!getMPV(L,chip,ch,m,entM)) continue;
      if(entM < cfg::MIP_MIN_ENT) continue;

      if(!isFinitePos(r.gain) || !isFinitePos(m)) continue;

      mpv=m; mpvEnt=entM; gain=r.gain; nphoto = mpv/gain;
      pedSigma = r.pedSigma;

      if(isFinitePos(nphoto)){
        tMain->Fill(); filledMain++;
        fillNphotoMap(mapMain.get(), L, chip, ch, nphoto, tagSPE, "main");
      }

      bool fitOK = (r.snr>=cfg::SNR_MIN_REF) && (r.err>0) && (r.gain>0) && ((r.err/r.gain)<=cfg::RELERR_MAX_REF);
      bool pedOK = (r.pedOK!=0);
      if(fitOK && pedOK && isFinitePos(nphoto)){
        tRef->Fill(); filledRef++;
        fillNphotoMap(mapRef.get(), L, chip, ch, nphoto, tagSPE, "ref");
        if (pedSigma < 50.0){
          tRef50->Fill(); filledRef50++;
          fillNphotoMap(mapRef50.get(), L, chip, ch, nphoto, tagSPE, "refPed50");
        }
      }
    }

    tMain->Write("", TObject::kOverwrite);
    tRef ->Write("", TObject::kOverwrite);
    tRef50->Write("", TObject::kOverwrite);

    Info("Nphoto","%s: main=%d, ref=%d, refPed50=%d", tagSPE, filledMain, filledRef, filledRef50);
  };

  writeNphoto("SingleHit", resS, "Nphoto_main_SingleHit", "Nphoto_ref_SingleHit", "Nphoto_refPed50_SingleHit");
  writeNphoto("AllHits"  , resA, "Nphoto_main_AllHits"  , "Nphoto_ref_AllHits"  , "Nphoto_refPed50_AllHits"  );

  // 低統計/欠落のサマリのみ表示
  Info("summary","Run=%s: missing(Single/All)=%ld/%ld, lowstat(Single/All)=%ld/%ld",
       runFile.Data(), rc.missing_single, rc.missing_all, rc.lowstat_single, rc.lowstat_all);

  Info("RUN","Done run: %s", runFile.Data());
}