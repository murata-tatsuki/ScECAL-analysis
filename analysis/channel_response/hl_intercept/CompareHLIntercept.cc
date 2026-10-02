// Compare saved ChannelResponse histograms; do not modify production calibration.
#include <TFile.h>
#include <TTree.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLeaf.h>
#include <array>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace hl {
    namespace fs=std::filesystem;
    constexpr double NaN=std::numeric_limits<double>::quiet_NaN();
    const std::array<std::string,3> labels= {
        "data","before","after"
    };
    const std::array<int,3> colors= {
        kBlack,kOrange+7,kAzure+1
    };
    void require(bool ok,const std::string& m) {
        if(!ok)throw std::runtime_error(m);
    }
    struct Options {
        std::array<std::string,3> inputs;
        std::string dataBefore,reference,output,energy="100";
        double low=800,high=2200,minBin=20,minSpan=400,minHits=200,pivot=1500,switchMargin=600,gainSwitch=2600;
        int minBins=8,top=20;
        std::vector<int> channels;
    };
    std::unique_ptr<TFile> open(const std::string& path) {
        auto f=std::unique_ptr<TFile>(TFile::Open(path.c_str(),"READ"));
        require(f&&!f->IsZombie(),"Cannot open "+path);
        return f;
    }
    template<class T>T* get(TDirectory* d,const std::string& name) {
        auto*o=dynamic_cast<T*>(d->Get(name.c_str()));
        require(o,"Missing/wrong object "+std::string(d->GetPath())+"/"+name);
        return o;
    }
    template<class T>std::unique_ptr<T> hist(TDirectory*d,const std::string& name) {
        auto*p=get<T>(d,name);
        p->SetDirectory(nullptr);
        return std::unique_ptr<T>(p);
    }
    std::string named(TFile& f,const char*name) {
        return get<TNamed>(&f,name)->GetTitle();
    }
    void sameAxis(const TAxis&a,const TAxis&b) {
        require(a.GetNbins()==b.GetNbins(),"Different histogram bin counts");
        for(int i=1;i<=a.GetNbins()+1;++i)require(a.GetBinLowEdge(i)==b.GetBinLowEdge(i),"Different histogram bin edges");
    }
    void sameData(const TH1&a,const TH1&b) {
        sameAxis(*a.GetXaxis(),*b.GetXaxis());
        sameAxis(*a.GetYaxis(),*b.GetYaxis());
        require(a.GetNcells()==b.GetNcells()&&a.GetEntries()==b.GetEntries(),"Data changed between before/after");
        for(int i=0;i<a.GetNcells();++i)require(a.GetBinContent(i)==b.GetBinContent(i)&&a.GetBinError(i)==b.GetBinError(i),"Data bin content/error changed between before/after");
        auto*p=dynamic_cast<const TProfile*>(&a);
        auto*q=dynamic_cast<const TProfile*>(&b);
        if(p&&q)for(int i=0;i<a.GetNcells();++i)require(p->GetBinEntries(i)==q->GetBinEntries(i),"Data profile statistics changed");
        double x[13]= {
        },y[13]= {
        };
        a.GetStats(x);
        b.GetStats(y);
        for(int i=0;i<13;++i)require(x[i]==y[i],"Data unbinned statistics changed");
    }
    struct Ref {
        double slope=NaN,b=NaN,xmax=NaN;
    };
    std::map<int,Ref> readReference(const std::string& path) {
        auto f=open(path);
        auto*t=get<TTree>(f.get(),"InterCalib");
        std::map<std::string,std::string> types= {
            {
                "CellID","Int_t"
            }, {
                "Slope","Double_t"
            }, {
                "Intercept","Double_t"
            }, {
                "SlopeError","Double_t"
            }, {
                "InterceptError","Double_t"
            }, {
                "ChiSquare","Double_t"
            }, {
                "NDF","Double_t"
            }, {
                "XMax","Double_t"
            }, {
                "Statistics","Int_t"
            }
        };
        require(t->GetListOfBranches()->GetEntries()==9,"Reference InterCalib must have exactly 9 branches");
        for(auto&[n,ty]:types) {
            auto*l=t->GetLeaf(n.c_str());
            require(l&&ty==l->GetTypeName(),"Reference branch schema mismatch: "+n);
        }
        int id;
        Ref r;
        t->SetBranchAddress("CellID",&id);
        t->SetBranchAddress("Slope",&r.slope);
        t->SetBranchAddress("Intercept",&r.b);
        t->SetBranchAddress("XMax",&r.xmax);
        std::map<int,Ref> a;
        for(Long64_t i=0;i<t->GetEntries();++i) {
            t->GetEntry(i);
            require(!a.count(id),"Duplicate reference CellID");
            a[id]=r;
        }
        return a;
    }
    struct Fit {
        int status=1,bins=0;
        double s=NaN,b=NaN,se=NaN,be=NaN,cov=NaN,chi=NaN,ndf=0,n=0,xlo=NaN,xhi=NaN,atPivot=NaN,pivotError=NaN;
    };
    // Status: 0=valid, 1=missing channel, 2=insufficient common bins, 3=insufficient
    // span, 4=insufficient hits, 5=singular/nonfinite/nonpositive slope, 6=no reference.
    Fit fitProfile(const TProfile&p,const std::vector<int>&bins,const Options&o) {
        Fit f;
        f.bins=bins.size();
        f.status=2;
        if(bins.empty())return f;
        f.xlo=p.GetBinCenter(bins.front());
        f.xhi=p.GetBinCenter(bins.back());
        for(int i:bins)f.n+=p.GetBinEntries(i);
        if(f.bins<o.minBins)return f;
        f.status=3;
        if(f.xhi-f.xlo<o.minSpan)return f;
        f.status=4;
        if(f.n<o.minHits)return f;
        double sw=0,sx=0,sy=0;
        for(int i:bins) {
            double w=1/std::pow(p.GetBinError(i),2);
            sw+=w;
            sx+=w*p.GetBinCenter(i);
            sy+=w*p.GetBinContent(i);
        }
        double mx=sx/sw,my=sy/sw,sxx=0,sxy=0;
        for(int i:bins) {
            double w=1/std::pow(p.GetBinError(i),2),dx=p.GetBinCenter(i)-mx;
            sxx+=w*dx*dx;
            sxy+=w*dx*(p.GetBinContent(i)-my);
        }
        f.status=5;
        if(!(sxx>0)||!std::isfinite(sxx))return f;
        f.s=sxy/sxx;
        f.b=my-f.s*mx;
        f.se=std::sqrt(1/sxx);
        f.be=std::sqrt(1/sw+mx*mx/sxx);
        f.cov=-mx/sxx;
        f.chi=0;
        f.ndf=f.bins-2;
        for(int i:bins)f.chi+=std::pow((p.GetBinContent(i)-f.s*p.GetBinCenter(i)-f.b)/p.GetBinError(i),2);
        f.atPivot=f.s*o.pivot+f.b;
        f.pivotError=std::sqrt(1/sw+std::pow(o.pivot-mx,2)/sxx);
        if(std::isfinite(f.s)&&f.s>0&&std::isfinite(f.b)&&std::isfinite(f.chi))f.status=0;
        return f;
    }
    std::vector<int> commonBins(const std::array<std::unique_ptr<TProfile>,3>&ps,double lo,double hi,const Options&o) {
        std::vector<int> bins;
        for(int i=1;i<=ps[0]->GetNbinsX();++i) {
            if(ps[0]->GetXaxis()->GetBinLowEdge(i)<lo||ps[0]->GetXaxis()->GetBinUpEdge(i)>hi)continue;
            bool valid=true;
            for(auto&p:ps)valid=valid&&p->GetBinEntries(i)>=o.minBin&&p->GetBinError(i)>0&&std::isfinite(p->GetBinError(i))&&std::isfinite(p->GetBinContent(i));
            if(valid)bins.push_back(i);
        }
        return bins;
    }
    struct Distance {
        double bias=NaN,rms=NaN,weight=0;
    };
    std::array<Distance,2> distances(const std::array<std::unique_ptr<TProfile>,3>&ps,const std::vector<int>&bins) {
        std::array<Distance,2>d;
        double sw=0;
        std::array<double,2>sy= {
        },sy2= {
        };
        for(int i:bins) {
            double w=std::min( {
                ps[0]->GetBinEntries(i),ps[1]->GetBinEntries(i),ps[2]->GetBinEntries(i)
            });
            sw+=w;
            for(int k=0;k<2;++k) {
                double delta=ps[k+1]->GetBinContent(i)-ps[0]->GetBinContent(i);
                sy[k]+=w*delta;
                sy2[k]+=w*delta*delta;
            }
        }
        for(int k=0;k<2;++k) {
            d[k].weight=sw;
            if(sw) {
                d[k].bias=sy[k]/sw;
                d[k].rms=std::sqrt(sy2[k]/sw);
            }
        }
        return d;
    }
    struct Residual {
        double n=0,flow=0,mean=NaN,rms=NaN,q16=NaN,median=NaN,q84=NaN,width=NaN,tail=NaN;
        std::vector<std::pair<double,double>> points;
    };
    double quantile(const std::vector<std::pair<double,double>>&v,double n,double q) {
        double sum=0;
        for(auto [x,w]:v) {
            sum+=w;
            if(sum>=q*n)return x;
        }
        return NaN;
    }
    Residual residual(const TH2D&h,const Fit&f,const std::vector<int>&bins) {
        Residual r;
        if(f.status)return r;
        double sum=0,sum2=0;
        for(int i:bins) {
            r.flow+=h.GetBinContent(i,0)+h.GetBinContent(i,h.GetNbinsY()+1);
            for(int j=1;j<=h.GetNbinsY();++j) {
                double n=h.GetBinContent(i,j);
                if(!n)continue;
                double v=h.GetYaxis()->GetBinCenter(j)-f.s*h.GetXaxis()->GetBinCenter(i)-f.b;
                r.points.emplace_back(v,n);
                r.n+=n;
                sum+=n*v;
                sum2+=n*v*v;
            }
        }
        if(r.n) {
            r.mean=sum/r.n;
            r.rms=std::sqrt(std::max(0.,sum2/r.n-r.mean*r.mean));
            std::sort(r.points.begin(),r.points.end());
            r.q16=quantile(r.points,r.n,.16);
            r.median=quantile(r.points,r.n,.5);
            r.q84=quantile(r.points,r.n,.84);
            r.width=(r.q84-r.q16)/2;
        }
        return r;
    }
    struct Stats {
        double hits=NaN,perEvent=NaN,hgMean=NaN,lgMean=NaN,outside=NaN,energyMean=NaN,energyPerEvent=NaN,energyOutside=NaN;
    };
    Stats stats(const TH2D&raw,const TH2D&energy,double events) {
        Stats s;
        double x[7],y[7];
        raw.GetStats(x);
        energy.GetStats(y);
        s.hits=raw.GetEntries();
        s.perEvent=s.hits/events;
        s.outside=s.hits?1-x[0]/s.hits:NaN;
        s.hgMean=x[0]?x[2]/x[0]:NaN;
        s.lgMean=x[0]?x[4]/x[0]:NaN;
        s.energyMean=y[0]?y[4]/y[0]:NaN;
        s.energyPerEvent=y[4]/events;
        s.energyOutside=energy.GetEntries()?1-y[0]/energy.GetEntries():NaN;
        return s;
    }
    struct SwitchCheck {
        int leftStatus=1,rightStatus=1;
        double cut=NaN,jump=NaN,error=NaN;
    };
    SwitchCheck switchCheck(const TProfile& p,double pedestal,const Options&o) {
        SwitchCheck result;
        result.cut=o.gainSwitch-pedestal;
        Options local=o;
        local.minBins=3;
        local.minSpan=60;
        local.minHits=60;
        local.pivot=result.cut;
        std::array<std::vector<int>,2> bins;
        for(int i=1;i<=p.GetNbinsX();++i) {
            if(p.GetBinEntries(i)<o.minBin||!(p.GetBinError(i)>0)||!std::isfinite(p.GetBinError(i)))continue;
            double lo=p.GetXaxis()->GetBinLowEdge(i),hi=p.GetXaxis()->GetBinUpEdge(i);
            if(lo>=result.cut-250&&hi<=result.cut)bins[0].push_back(i);
            if(lo>=result.cut&&hi<=result.cut+250)bins[1].push_back(i);
        }
        auto left=fitProfile(p,bins[0],local),right=fitProfile(p,bins[1],local);
        result.leftStatus=left.status;
        result.rightStatus=right.status;
        if(!left.status&&!right.status) {
            result.jump=right.atPivot-left.atPivot;
            result.error=std::hypot(left.pivotError,right.pivotError);
        }
        return result;
    }
    // Total variation between unit-area projections, including both axes' flow bins.
    double projectionTV(const TH2D&a,const TH2D&b,bool hg) {
        sameAxis(*a.GetXaxis(),*b.GetXaxis());
        sameAxis(*a.GetYaxis(),*b.GetYaxis());
        if(!(a.GetEntries()>0&&b.GetEntries()>0))return NaN;
        int n=hg?a.GetNbinsX()+2:a.GetNbinsY()+2,m=hg?a.GetNbinsY()+2:a.GetNbinsX()+2;
        double tv=0;
        for(int i=0;i<n;++i) {
            double x=0,y=0;
            for(int j=0;j<m;++j) {
                x+=hg?a.GetBinContent(i,j):a.GetBinContent(j,i);
                y+=hg?b.GetBinContent(i,j):b.GetBinContent(j,i);
            }
            tv+=std::abs(x/a.GetEntries()-y/b.GetEntries());
        }
        return tv/2;
    }
    struct Row {
        int id=0;
        Ref ref;
        double upper=NaN;
        std::array<Fit,3>fit,inner;
        std::array<Stats,3>stat;
        std::array<Residual,3>res;
        std::array<Distance,2>rawDiff,energyDiff;
        double hgTV=NaN,lgTV=NaN;
        std::array<SwitchCheck,3>sw;
    };
    struct FitWriter {
        std::unique_ptr<TFile> file;
        TTree *tree,*diag;
        int id=0,statistics=0,status=0,bins=0;
        double s,b,se,be,chi,ndf,xmax,cov,lo,hi,n,atPivot,pivotError;
        FitWriter(const fs::path& path,const Options&o,int sample) {
            file=std::make_unique<TFile>(path.c_str(),"RECREATE");
            require(!file->IsZombie(),"Cannot create fit ROOT");
            tree=new TTree("InterCalib","HG-LG profile linear fit; see provenance and FitDiagnostics");
            tree->Branch("CellID",&id,"CellID/I");
            tree->Branch("Slope",&s,"Slope/D");
            tree->Branch("Intercept",&b,"Intercept/D");
            tree->Branch("SlopeError",&se,"SlopeError/D");
            tree->Branch("InterceptError",&be,"InterceptError/D");
            tree->Branch("ChiSquare",&chi,"ChiSquare/D");
            tree->Branch("NDF",&ndf,"NDF/D");
            tree->Branch("XMax",&xmax,"XMax/D");
            tree->Branch("Statistics",&statistics,"Statistics/I");
            diag=new TTree("FitDiagnostics","All requested physical channels; InterCalib contains successful fits only");
            diag->Branch("CellID",&id,"CellID/I");
            diag->Branch("Status",&status,"Status/I");
            diag->Branch("Bins",&bins,"Bins/I");
            diag->Branch("Statistics",&n,"Statistics/D");
            diag->Branch("SlopeInterceptCovariance",&cov,"SlopeInterceptCovariance/D");
            diag->Branch("FitXMin",&lo,"FitXMin/D");
            diag->Branch("FitXMax",&hi,"FitXMax/D");
            diag->Branch("LGAtPivot",&atPivot,"LGAtPivot/D");
            diag->Branch("LGAtPivotError",&pivotError,"LGAtPivotError/D");
            TNamed("sample",labels[sample].c_str()).Write();
            TNamed("source",o.inputs[sample].c_str()).Write();
            TNamed("reference",o.reference.c_str()).Write();
            TNamed("method","Weighted least squares of mean LG against HG bin centers; common bins among data/before/after, inverse SEM squared weights. Errors are formal, not scaled by chi2. Not the original multigaus fit.").Write();
            TNamed("XMax_policy","Copied from reference InterCalib, NOT measured here. Actual fit support is FitDiagnostics/FitXMin,FitXMax. Statistics is hits in accepted common fit bins.").Write();
            TNamed("status_codes","0=valid;1=missing_channel;2=insufficient_common_bins;3=insufficient_HG_span;4=insufficient_hits;5=invalid_fit;6=missing_reference").Write();
            TNamed("limitations","Profile bin-center approximation; no per-hit mean HG in each bin. SEM errors exclude binning/range/model systematics. No direct replacement of production calibration is performed.").Write();
            TParameter<double>("pivot_HG_ADC",o.pivot).Write();
        }
        void fill(const Row&r,int k) {
            auto&f=r.fit[k];
            id=r.id;
            status=f.status;
            bins=f.bins;
            n=f.n;
            cov=f.cov;
            lo=f.xlo;
            hi=f.xhi;
            atPivot=f.atPivot;
            pivotError=f.pivotError;
            diag->Fill();
            if(status)return;
            s=f.s;
            b=f.b;
            se=f.se;
            be=f.be;
            chi=f.chi;
            ndf=f.ndf;
            xmax=r.ref.xmax;
            require(n<=std::numeric_limits<int>::max(),"Statistics exceeds reference Int_t capacity");
            statistics=std::lround(n);
            tree->Fill();
        }
        void finish() {
            file->cd();
            tree->Write();
            diag->Write();
            file->Close();
        }
    };
    std::string statusName(int s) {
        static const char*n[]= {
            "valid","missing_channel","insufficient_common_bins","insufficient_HG_span","insufficient_hits","invalid_fit","missing_reference"
        };
        return n[s];
    }
    void writeTables(const std::vector<Row>&rows,const fs::path&dir,const Options&o) {
        std::ofstream f(dir/"channels.tsv"),c(dir/"comparison.tsv");
        f<<std::setprecision(12);
        c<<std::setprecision(12);
        f<<"energy_GeV\tcellid\tlayer\tchip\tchannel\tsample\tstatus\tslope\tintercept_LG_ADC\tslope_error\tintercept_error\tcov_slope_intercept\tchi2\tndf\tchi2_ndf\tfit_bins\tfit_hits\tfit_HG_min\tfit_HG_max\tupper_limit\tLG_at_pivot\tLG_at_pivot_error\tinner_status\tinner_slope\tinner_intercept\tinner_minus_base_at_pivot\tref_slope\tref_intercept\tref_XMax\thits\thits_per_event\tHG_mean_inrange\tLG_mean_inrange\traw_outside_fraction\tselected_energy_mean_inrange_MeV\tselected_energy_sum_inrange_per_event_MeV\tenergy_outside_fraction\tbinned_residual_hits\tbinned_residual_flow\tbinned_residual_mean_ADC\tbinned_residual_RMS_ADC\tbinned_residual_median_ADC\tbinned_residual_width68_ADC\tbinned_tail_fraction_data5width68\tswitch_HG_pedsub_ADC\tswitch_left_status\tswitch_right_status\tswitch_jump_proxy_MeV\tswitch_jump_formal_error_MeV\n";
        c<<"energy_GeV\tcellid\tall_fits_valid\tcommon_weight\tbias_before_LG_ADC\tbias_after_LG_ADC\tRMS_before_LG_ADC\tRMS_after_LG_ADC\tRMS_improvement_LG_ADC\tdelta_slope\tdelta_intercept\tinput_intercept\tdelta_intercept_minus_input\tbefore_intercept\tafter_intercept_minus_input\tafter_slope_minus_input\tbefore_minus_data_slope\tafter_minus_data_slope\tbefore_minus_data_at_pivot\tafter_minus_data_at_pivot\tenergy_bias_before_MeV\tenergy_bias_after_MeV\tenergy_RMS_before_MeV\tenergy_RMS_after_MeV\tHG_projection_TV_before_after\tLG_projection_TV_before_after\n";
        for(auto&r:rows) {
            for(int k=0;k<3;++k) {
                auto&a=r.fit[k];
                auto&b=r.inner[k];
                auto&s=r.stat[k];
                auto&t=r.res[k];
                f<<o.energy<<'\t'<<r.id<<'\t'<<r.id/100000<<'\t'<<r.id%100000/10000<<'\t'<<r.id%100<<'\t'<<labels[k]<<'\t'<<statusName(a.status)<<'\t'<<a.s<<'\t'<<a.b<<'\t'<<a.se<<'\t'<<a.be<<'\t'<<a.cov<<'\t'<<a.chi<<'\t'<<a.ndf<<'\t'<<(a.ndf?a.chi/a.ndf:NaN)<<'\t'<<a.bins<<'\t'<<a.n<<'\t'<<a.xlo<<'\t'<<a.xhi<<'\t'<<r.upper<<'\t'<<a.atPivot<<'\t'<<a.pivotError<<'\t'<<statusName(b.status)<<'\t'<<b.s<<'\t'<<b.b<<'\t'<<b.atPivot-a.atPivot<<'\t'<<r.ref.slope<<'\t'<<r.ref.b<<'\t'<<r.ref.xmax<<'\t'<<s.hits<<'\t'<<s.perEvent<<'\t'<<s.hgMean<<'\t'<<s.lgMean<<'\t'<<s.outside<<'\t'<<s.energyMean<<'\t'<<s.energyPerEvent<<'\t'<<s.energyOutside<<'\t'<<t.n<<'\t'<<t.flow<<'\t'<<t.mean<<'\t'<<t.rms<<'\t'<<t.median<<'\t'<<t.width<<'\t'<<t.tail<<'\t'<<r.sw[k].cut<<'\t'<<statusName(r.sw[k].leftStatus)<<'\t'<<statusName(r.sw[k].rightStatus)<<'\t'<<r.sw[k].jump<<'\t'<<r.sw[k].error<<'\n';
            }
            auto&a=r.fit[0];
            auto&b=r.fit[1];
            auto&z=r.fit[2];
            c<<o.energy<<'\t'<<r.id<<'\t'<<(!a.status&&!b.status&&!z.status)<<'\t'<<r.rawDiff[0].weight<<'\t'<<r.rawDiff[0].bias<<'\t'<<r.rawDiff[1].bias<<'\t'<<r.rawDiff[0].rms<<'\t'<<r.rawDiff[1].rms<<'\t'<<r.rawDiff[0].rms-r.rawDiff[1].rms<<'\t'<<z.s-b.s<<'\t'<<z.b-b.b<<'\t'<<r.ref.b<<'\t'<<z.b-b.b-r.ref.b<<'\t'<<b.b<<'\t'<<z.b-r.ref.b<<'\t'<<z.s-r.ref.slope<<'\t'<<b.s-a.s<<'\t'<<z.s-a.s<<'\t'<<b.atPivot-a.atPivot<<'\t'<<z.atPivot-a.atPivot<<'\t'<<r.energyDiff[0].bias<<'\t'<<r.energyDiff[1].bias<<'\t'<<r.energyDiff[0].rms<<'\t'<<r.energyDiff[1].rms<<'\t'<<r.hgTV<<'\t'<<r.lgTV<<'\n';
        }
        require(bool(f)&&bool(c),"Cannot write tables");
    }
    void overview(const std::vector<Row>&rows,const fs::path&dir) {
        gStyle->SetTitleAlign(23);
        gStyle->SetTitleX(.5);
        gStyle->SetTitleY(.98);
        gStyle->SetTitleH(.05);
        gStyle->SetTitleFontSize(.035);
        gStyle->SetPadTopMargin(.12);
        gStyle->SetPadBottomMargin(.14);
        gStyle->SetPadLeftMargin(.14);
        TFile out((dir/"summary.root").c_str(),"RECREATE");
        require(!out.IsZombie(),"Cannot create summary ROOT");
        TH2D improvement("RMS_improvement","RMS(before-data) - RMS(after-data);Layer;Chip*36 + channel;ADC",30,-.5,29.5,216,-.5,215.5),closure("intercept_closure","Delta intercept - input intercept;Layer;Chip*36 + channel;ADC",30,-.5,29.5,216,-.5,215.5),mask("valid_common_fit","Three valid fits (1=yes, 0=no);Layer;Chip*36 + channel",30,-.5,29.5,216,-.5,215.5);
        for(int i=0;i<improvement.GetNcells();++i) {
            improvement.SetBinContent(i,NaN);
            closure.SetBinContent(i,NaN);
        }
        TGraphErrors delta,slope,before,after;
        delta.SetName("delta_intercept_vs_input");
        slope.SetName("delta_slope_vs_cellid");
        before.SetName("bias_before_vs_cellid");
        after.SetName("bias_after_vs_cellid");
        TH1D improvements("RMS_improvement_distribution","Channel count;RMS improvement [ADC];Channels",200,-20,20),closureDist("intercept_closure_distribution","Channel count;Delta intercept - input [ADC];Channels",200,-10,10);
        for(auto&r:rows) {
            int x=r.id/100000+1,y=r.id%100000/10000*36+r.id%100+1;
            bool valid=!r.fit[0].status&&!r.fit[1].status&&!r.fit[2].status;
            mask.SetBinContent(x,y,valid);
            improvement.SetBinContent(x,y,NaN);
            closure.SetBinContent(x,y,NaN);
            if(std::isfinite(r.rawDiff[0].rms)) {
                double v=r.rawDiff[0].rms-r.rawDiff[1].rms;
                improvement.SetBinContent(x,y,v);
                improvements.Fill(v);
                int n=before.GetN();
                before.SetPoint(n,r.id,r.rawDiff[0].bias);
                after.SetPoint(n,r.id,r.rawDiff[1].bias);
            }
            if(!valid)continue;
            auto&b=r.fit[1];
            auto&a=r.fit[2];
            double d=a.b-b.b-r.ref.b;
            closure.SetBinContent(x,y,d);
            closureDist.Fill(d);
            int n=delta.GetN();
            delta.SetPoint(n,r.ref.b,a.b-b.b);
            delta.SetPointError(n,0,std::hypot(a.be,b.be));
            slope.SetPoint(n,r.id,a.s-b.s);
            slope.SetPointError(n,0,std::hypot(a.se,b.se));
        }
        improvement.Write();
        closure.Write();
        mask.Write();
        improvements.Write();
        closureDist.Write();
        delta.Write();
        slope.Write();
        before.Write();
        after.Write();
        TNamed("error_assumption","Before/after errors combined as independent; covariance between regenerated samples is unknown. Common-fit validity is supplied separately; invalid map cells are NaN.").Write();
        TH2D displayImprovement(improvement),displayClosure(closure);
        for(auto*h: {
            &displayImprovement,&displayClosure
        }) {
            double lo=std::numeric_limits<double>::infinity(),hi=-lo;
            for(int i=1;i<=h->GetNbinsX();++i)for(int j=1;j<=h->GetNbinsY();++j) {
                double z=h->GetBinContent(i,j);
                if(std::isfinite(z)) {
                    lo=std::min(lo,z);
                    hi=std::max(hi,z);
                }
            }
            if(!std::isfinite(lo)) {
                lo=0;
                hi=1;
            }
            if(!(hi>lo)) {
                lo-=.5;
                hi+=.5;
            }
            h->SetMinimum(lo);
            h->SetMaximum(hi);
            for(int i=0;i<h->GetNcells();++i)if(!std::isfinite(h->GetBinContent(i)))h->SetBinContent(i,lo-std::max(1.,hi-lo));
        }
        TCanvas cv("overview","HL intercept overview",1500,1000);
        cv.Divide(2,2);
        cv.cd(1);
        gPad->SetRightMargin(.18);
        displayImprovement.Draw("COLZ");
        cv.cd(2);
        gPad->SetRightMargin(.18);
        displayClosure.Draw("COLZ");
        cv.cd(3);
        delta.SetTitle("Input vs fitted shift;Input intercept [LG ADC];After - before intercept [LG ADC]");
        delta.SetMarkerStyle(20);
        delta.SetMarkerSize(.3);
        std::unique_ptr<TLine> diagonal;
        if(delta.GetN()) {
            delta.Draw("AP");
            double lo=1e10,hi=-1e10;
            for(int i=0;i<delta.GetN();++i) {
                lo=std::min(lo,delta.GetPointX(i));
                hi=std::max(hi,delta.GetPointX(i));
            }
            diagonal=std::make_unique<TLine>(lo,lo,hi,hi);
            diagonal->SetLineStyle(2);
            diagonal->Draw();
        }
        cv.cd(4);
        before.SetTitle("Data agreement;Physical CellID;MC - data mean LG [ADC]");
        before.SetMarkerStyle(20);
        after.SetMarkerStyle(24);
        before.SetMarkerColor(colors[1]);
        after.SetMarkerColor(colors[2]);
        before.SetMarkerSize(.3);
        after.SetMarkerSize(.3);
        if(before.GetN()) {
            double lo=0,hi=0;
            for(int i=0;i<before.GetN();++i) {
                lo=std::min( {
                    lo,before.GetPointY(i),after.GetPointY(i)
                });
                hi=std::max( {
                    hi,before.GetPointY(i),after.GetPointY(i)
                });
            }
            double margin=std::max(.5,.1*(hi-lo));
            before.SetMinimum(lo-margin);
            before.SetMaximum(hi+margin);
            before.Draw("AP");
            after.Draw("P SAME");
        }
        TLegend leg(.6,.75,.9,.9);
        leg.AddEntry(&before,"Before","p");
        leg.AddEntry(&after,"After","p");
        leg.Draw();
        cv.SaveAs((dir/"overview.png").c_str());
        cv.SaveAs((dir/"overview.pdf").c_str());
        cv.Write();
        out.Close();
    }
    void detail(const Row&r,std::array<std::unique_ptr<TFile>,3>&files,const fs::path&dir,const Options&o) {
        TCanvas c("detail","Channel comparison",1600,1100);
        c.Divide(2,2);
        std::array<TGraphErrors,3>g,dev;
        std::array<std::unique_ptr<TH1D>,3>rh;
        std::vector<std::unique_ptr<TLine>> lines;
        TLegend legend(.14,.72,.42,.88);
        legend.SetBorderSize(0);
        double ymin=1e10,ymax=-1e10;
        for(int k=0;k<3;++k) {
            auto*d=files[k]->GetDirectory(("channel_"+std::to_string(r.id)).c_str());
            if(!d)continue;
            auto p=hist<TProfile>(d,"hg_lg_zoom_profile");
            for(int i=1;i<=p->GetNbinsX();++i) {
                double x=p->GetBinCenter(i),y=p->GetBinContent(i);
                if(x<o.low||x>r.upper||p->GetBinEntries(i)<o.minBin)continue;
                int n=g[k].GetN();
                g[k].SetPoint(n,x,y);
                g[k].SetPointError(n,0,p->GetBinError(i));
                ymin=std::min(ymin,y);
                ymax=std::max(ymax,y);
                if(!r.fit[k].status) {
                    dev[k].SetPoint(n,x,y-r.fit[k].s*x-r.fit[k].b);
                    dev[k].SetPointError(n,0,p->GetBinError(i));
                }
            }
            for(auto*gr: {
                &g[k],&dev[k]
            }) {
                gr->SetMarkerStyle(20+k);
                gr->SetMarkerColor(colors[k]);
                gr->SetLineColor(colors[k]);
                gr->SetMarkerSize(.5);
            }
            legend.AddEntry(&g[k],labels[k].c_str(),"p");
            rh[k]=std::make_unique<TH1D>(("res_"+labels[k]).c_str(),"Binned residual, visible Y bins only;LG - fitted LG [ADC];Fraction",240,-30,30);
            rh[k]->SetDirectory(nullptr);
            for(auto[v,w]:r.res[k].points)rh[k]->Fill(v,w);
            if(r.res[k].n)rh[k]->Scale(1/r.res[k].n);
            rh[k]->SetLineColor(colors[k]);
        }
        c.cd(1);
        if(ymax>ymin) {
            gPad->DrawFrame(o.low,ymin-2,r.upper,ymax+2,("Cell "+std::to_string(r.id)+";HG - pedestal [ADC];Mean LG - pedestal [ADC]").c_str());
            for(int k=0;k<3;++k) {
                g[k].Draw("P SAME");
                if(!r.fit[k].status) {
                    auto&f=r.fit[k];
                    lines.emplace_back(new TLine(f.xlo,f.s*f.xlo+f.b,f.xhi,f.s*f.xhi+f.b));
                    lines.back()->SetLineColor(colors[k]);
                    lines.back()->Draw();
                }
            }
            legend.Draw();
        }
        c.cd(2);
        gPad->DrawFrame(o.low,-5,r.upper,5,"Profile residual (own fit);HG - pedestal [ADC];Mean residual [LG ADC]");
        for(auto&gr:dev)gr.Draw("P SAME");
        c.cd(3);
        gPad->SetLogy();
        bool drawn=false;
        double max=0;
        for(auto&h:rh)if(h)max=std::max(max,h->GetMaximum());
        for(auto&h:rh)if(h) {
            h->SetMaximum(std::max(.01,max*1.5));
            h->SetMinimum(1e-6);
            h->Draw(drawn?"HIST SAME":"HIST");
            drawn=true;
        }
        // Keep the lower-right comparison area, with intercept and slope side by side.
        auto* parameters = c.cd(4);
        parameters->Divide(2, 1, .01, .01);
        std::array<std::unique_ptr<TH1D>, 2> frames;
        std::array<std::array<TGraphErrors, 3>, 2> points;
        std::array<std::unique_ptr<TLine>, 2> references;
        std::array<std::unique_ptr<TLegend>, 2> referenceLegends;
        std::array<TLatex, 2> changes;
        for (int parameter = 0; parameter < 2; ++parameter) {
            parameters->cd(parameter + 1);
            gPad->SetLeftMargin(.22);
            gPad->SetRightMargin(.06);
            gPad->SetBottomMargin(.18);
            gPad->SetTopMargin(.23);
            const bool isSlope = parameter == 1;
            const double input = isSlope ? r.ref.slope : r.ref.b;
            double minimum = std::numeric_limits<double>::infinity();
            double maximum = -minimum;
            if (std::isfinite(input)) minimum = maximum = input;
            for (int k = 0; k < 3; ++k) if (!r.fit[k].status) {
                const auto& fit = r.fit[k];
                const double value = isSlope ? fit.s : fit.b;
                const double error = isSlope ? fit.se : fit.be;
                minimum = std::min(minimum, value - error);
                maximum = std::max(maximum, value + error);
                auto& graph = points[parameter][k];
                graph.SetPoint(0, k, value);
                graph.SetPointError(0, 0, error);
                graph.SetMarkerStyle(20 + k);
                graph.SetMarkerColor(colors[k]);
                graph.SetLineColor(colors[k]);
                graph.SetMarkerSize(.9);
            }
            if (!std::isfinite(minimum)) { minimum = 0; maximum = 1; }
            const double margin = std::max(isSlope ? 1e-5 : 1., .15 * (maximum - minimum));
            const std::string name = isSlope ? "slope_categories" : "intercept_categories";
            frames[parameter] = std::make_unique<TH1D>(name.c_str(),
                isSlope ? "Fit slope;;Slope [LG / HG]" : "Fit intercept;;Intercept [LG ADC]",
                3, -.5, 2.5);
            auto& frame = *frames[parameter];
            frame.SetDirectory(nullptr);
            frame.SetMinimum(minimum - margin);
            frame.SetMaximum(maximum + margin);
            frame.SetStats(false);
            for (int k = 0; k < 3; ++k) frame.GetXaxis()->SetBinLabel(k + 1, labels[k].c_str());
            frame.GetXaxis()->LabelsOption("h");
            frame.GetXaxis()->SetLabelFont(43);
            frame.GetXaxis()->SetLabelSize(19);
            frame.GetXaxis()->SetLabelOffset(.015);
            frame.GetYaxis()->SetLabelFont(43);
            frame.GetYaxis()->SetLabelSize(15);
            frame.GetYaxis()->SetTitleFont(43);
            frame.GetYaxis()->SetTitleSize(18);
            frame.GetYaxis()->SetTitleOffset(2.4);
            frame.GetYaxis()->SetNdivisions(505);
            frame.Draw("AXIS");
            for (auto& graph : points[parameter]) if (graph.GetN()) graph.Draw("P SAME");
            if (std::isfinite(input)) {
                references[parameter] = std::make_unique<TLine>(-.5, input, 2.5, input);
                references[parameter]->SetLineStyle(2);
                references[parameter]->Draw();
                referenceLegends[parameter] = std::make_unique<TLegend>(.23, .87, .92, .92);
                auto& legend = *referenceLegends[parameter];
                legend.SetBorderSize(0); legend.SetFillStyle(0);
                legend.SetTextFont(43); legend.SetTextSize(15);
                legend.AddEntry(references[parameter].get(), "Input coefficient", "l");
                legend.Draw();
            }
            auto& heading = changes[parameter];
            heading.SetNDC(); heading.SetTextFont(43); heading.SetTextSize(20);
            heading.DrawLatex(.23, .955, isSlope ? "Fit slope" : "Fit intercept");
            if (!r.fit[1].status && !r.fit[2].status) {
                const double delta = isSlope ? r.fit[2].s - r.fit[1].s : r.fit[2].b - r.fit[1].b;
                auto& text = changes[parameter];
                text.SetNDC(); text.SetTextFont(43); text.SetTextSize(15);
                text.DrawLatex(.23, .81, Form("after - before = %.3g", delta));
            }
        }
        c.SaveAs((dir/("cell_"+std::to_string(r.id)+".png")).c_str());
    }
    void run(Options o) {
        gROOT->SetBatch(true);
        gStyle->SetOptStat(0);
        TH1::AddDirectory(false);
        require(o.low<o.high&&o.minBins>=3&&o.minBin>=2&&o.minHits>0&&o.minSpan>0&&o.top>=0,"Invalid fit options");
        if(o.channels.empty())for(int l=0;l<30;++l)for(int c=0;c<6;++c)for(int ch=0;ch<(c==5?30:36);++ch)o.channels.push_back(l*100000+c*10000+ch);
        auto refs=readReference(o.reference);
        std::array<std::unique_ptr<TFile>,3>files;
        std::array<double,3>events;
        for(int k=0;k<3;++k) {
            files[k]=open(o.inputs[k]);
            events[k]=get<TParameter<Long64_t>>(files[k].get(),"processed_events")->GetVal();
            require(events[k]>0,"No processed events");
        }
        auto dataOld=open(o.dataBefore);
        require(named(*files[0],"selection")==named(*dataOld,"selection"),"Data selection changed");
        require(named(*files[0],"input_pairs")==named(*dataOld,"input_pairs"),"Data input manifest changed");
        require(events[0]==get<TParameter<Long64_t>>(dataOld.get(),"processed_events")->GetVal(),"Data event count changed");
        for(int k=1;k<3;++k)require(named(*files[0],"selection")==named(*files[k],"selection"),"Sample selection metadata differs");
        fs::path base=o.output;
        require(!fs::exists(base),"Output exists; choose a new staging directory");
        auto out=base/"comparison"/(o.energy+"GeV");
        fs::create_directories(out/"channels");
        std::array<std::unique_ptr<FitWriter>,3>writers;
        for(int k=0;k<3;++k)writers[k]=std::make_unique<FitWriter>(base/("InterCalib_"+labels[k]+"_"+o.energy+"GeV.root"),o,k);
        std::vector<Row>rows;
        std::ofstream meta(out/"provenance.txt");
        meta<<std::setprecision(17)<<"energy_GeV="<<o.energy<<"\nreference="<<o.reference<<"\nfit_HG_requested="<<o.low<<","<<o.high<<"\nmin_bin_entries="<<o.minBin<<"\nmin_bins="<<o.minBins<<"\nmin_span="<<o.minSpan<<"\nmin_hits="<<o.minHits<<"\npivot="<<o.pivot<<"\nswitch_margin="<<o.switchMargin<<"\ngain_switch_raw_ADC="<<o.gainSwitch<<"\ninner_fit=requested_low+200,requested_high-200 (same channel switch cap)\n";
        for(int k=0;k<3;++k)meta<<labels[k]<<"="<<o.inputs[k]<<"\n"<<labels[k]<<"_events="<<events[k]<<"\n";
        const std::array<std::string,6>plots= {
            "hg_lg","hg_lg_zoom","hg_energy","lg_energy","hg_selected_energy","lg_selected_energy"
        };
        for(int id:o.channels) {
            Row r;
            r.id=id;
            bool ref=refs.count(id);
            if(ref)r.ref=refs.at(id);
            std::array<TDirectory*,3>ds {
            };
            for(int k=0;k<3;++k)ds[k]=files[k]->GetDirectory(("channel_"+std::to_string(id)).c_str());
            auto*old=dataOld->GetDirectory(("channel_"+std::to_string(id)).c_str());
            require(bool(old)==bool(ds[0]),"Data channel presence changed: "+std::to_string(id));
            if(old) {
                for(auto&name:plots) {
                    auto a=hist<TH2D>(old,name),b=hist<TH2D>(ds[0],name);
                    sameData(*a,*b);
                    auto p=hist<TProfile>(old,name+"_profile"),q=hist<TProfile>(ds[0],name+"_profile");
                    sameData(*p,*q);
                }
                for(auto name: {
                    "pedestal_hg","pedestal_lg"
                })require(get<TParameter<double>>(old,name)->GetVal()==get<TParameter<double>>(ds[0],name)->GetVal(),"Data pedestal changed");
            }
            std::array<std::unique_ptr<TProfile>,3>ps,es;
            std::array<std::unique_ptr<TH2D>,3>hs,raw;
            for(int k=0;k<3;++k)if(ds[k]) {
                ps[k]=hist<TProfile>(ds[k],"hg_lg_zoom_profile");
                es[k]=hist<TProfile>(ds[k],"hg_selected_energy_profile");
                hs[k]=hist<TH2D>(ds[k],"hg_lg_zoom");
                raw[k]=hist<TH2D>(ds[k],"hg_lg");
                auto en=hist<TH2D>(ds[k],"hg_selected_energy");
                r.stat[k]=stats(*raw[k],*en,events[k]);
                r.sw[k]=switchCheck(*es[k],get<TParameter<double>>(ds[k],"pedestal_hg")->GetVal(),o);
            }
            if(raw[1]&&raw[2]) {
                r.hgTV=projectionTV(*raw[1],*raw[2],true);
                r.lgTV=projectionTV(*raw[1],*raw[2],false);
            }
            bool all=ds[0]&&ds[1]&&ds[2];
            if(all) {
                for(int k=1;k<3;++k) {
                    sameAxis(*ps[0]->GetXaxis(),*ps[k]->GetXaxis());
                    sameAxis(*es[0]->GetXaxis(),*es[k]->GetXaxis());
                }
                r.upper=o.high;
                // MCDigi compares RAW HG against XMax-600. Profiles have pedestal removed.
                if(ref&&std::isfinite(r.ref.xmax)) {
                    double bw=ps[0]->GetXaxis()->GetBinWidth(1);
                    for(int k=1;k<3;++k)r.upper=std::min(r.upper,r.ref.xmax-o.switchMargin-get<TParameter<double>>(ds[k],"pedestal_hg")->GetVal()-bw);
                }
                auto bins=commonBins(ps,o.low,r.upper,o);
                auto inner=commonBins(ps,o.low+200,std::min(o.high-200,r.upper),o);
                r.rawDiff=distances(ps,bins);
                auto ebins=commonBins(es,o.low,r.upper,o);
                r.energyDiff=distances(es,ebins);
                for(int k=0;k<3;++k) {
                    r.fit[k]=fitProfile(*ps[k],bins,o);
                    r.inner[k]=fitProfile(*ps[k],inner,o);
                    if(!ref||!std::isfinite(r.ref.xmax)) {
                        r.fit[k].status=6;
                        r.inner[k].status=6;
                    }
                    r.res[k]=residual(*hs[k],r.fit[k],bins);
                }
                double threshold=5*r.res[0].width;
                if(threshold>0)for(auto&res:r.res)if(res.n) {
                    double n=0;
                    for(auto[v,w]:res.points)if(std::abs(v-res.median)>threshold)n+=w;
                    res.tail=n/res.n;
                }
            }
            for(int k=0;k<3;++k)writers[k]->fill(r,k);
            rows.push_back(std::move(r));
            if(rows.size()%500==0)std::cout<<"Processed "<<rows.size()<<" / "<<o.channels.size()<<" channels"<<std::endl;
        }
        for(auto&w:writers)w->finish();
        writeTables(rows,out,o);
        overview(rows,out);
        std::vector<const Row*>rank;
        for(auto&r:rows)if(!r.fit[0].status&&!r.fit[1].status&&!r.fit[2].status)rank.push_back(&r);
        std::sort(rank.begin(),rank.end(),[](auto*a,auto*b) {
            return a->rawDiff[1].rms>b->rawDiff[1].rms;
        });
        std::set<int>draw;
        for(int i=0;i<std::min(o.top,int(rank.size()));++i)draw.insert(rank[i]->id);
        std::sort(rank.begin(),rank.end(),[](auto*a,auto*b) {
            return a->rawDiff[1].rms-a->rawDiff[0].rms>b->rawDiff[1].rms-b->rawDiff[0].rms;
        });
        for(int i=0;i<std::min(o.top,int(rank.size()));++i)draw.insert(rank[i]->id);
        for(auto&r:rows)if(draw.count(r.id))detail(r,files,out/"channels",o);
        std::ofstream sum(out/"summary.txt");
        sum<<"All requested channels: "<<rows.size()<<"\nData before/after: identical histogram contents, statistics, profiles and pedestal; manifest/selection/event counts equal\n";
        for(int k=0;k<3;++k) {
            std::map<int,int>counts;
            int poor=0;
            for(auto&r:rows) {
                counts[r.fit[k].status]++;
                if(!r.fit[k].status&&r.fit[k].chi/r.fit[k].ndf>5)poor++;
            }
            sum<<labels[k]<<":";
            for(auto [s,n]:counts)sum<<" "<<statusName(s)<<"="<<n;
            sum<<" chi2_ndf_gt5="<<poor<<"\n";
        }
        int better=0,worse=0,n=0;
        double sumBefore=0,sumAfter=0,sw=0,wb=0,wa=0;
        for(auto&r:rows)if(!r.fit[0].status&&!r.fit[1].status&&!r.fit[2].status) {
            double b=r.rawDiff[0].rms,a=r.rawDiff[1].rms,w=r.rawDiff[0].weight;
            n++;
            better+=a<b;
            worse+=a>b;
            sumBefore+=b;
            sumAfter+=a;
            sw+=w;
            wb+=w*b;
            wa+=w*a;
        }
        sum<<"Three-valid-fit channels: "<<n<<" improved="<<better<<" worsened="<<worse<<" (descriptive, not significance tests)\n";
        if(n)sum<<"Equal-channel mean RMS difference [ADC]: before="<<sumBefore/n<<" after="<<sumAfter/n<<"\nCommon-bin-hit-weighted mean RMS difference [ADC]: before="<<wb/sw<<" after="<<wa/sw<<"\n";
        sum<<"Detailed plots: "<<draw.size()<<"; union of largest remaining mismatch and largest worsening.\nResidual width/tails are binned approximations, exclude Y flow; profile means include Y flow.\nEnergy sums/means are display-inrange only; no event resolution or gain-selection fraction can be inferred.\nFit coefficients use profile bin centers and include tails. Check range stability and chi2 before using coefficients.\nXMax is copied from reference. Missing/failed channels are absent from InterCalib and remain in FitDiagnostics and TSV.\n";
        meta<<"data_before="<<o.dataBefore<<"\ncompleted_channels="<<rows.size()<<"\n";
        require(bool(meta)&&bool(sum),"Failed writing summary");
        std::cout<<"Completed "<<rows.size()<<" channels; all-three-valid="<<n<<" output="<<base<<std::endl;
    }
    Options parse(int argc,char**argv) {
        Options o;
        for(int i=1;i<argc;++i) {
            std::string a=argv[i];
            if(a=="--help") {
                std::cout<<"CompareHLIntercept --data FILE --before FILE --after FILE --data-before FILE --reference FILE --output NEW_DIR [--energy 100 --low 800 --high 2200 --min-bin 20 --min-bins 8 --min-span 400 --min-hits 200 --pivot 1500 --switch-margin 600 --gain-switch 2600 --top 20 --channels IDS]\n";
                std::exit(0);
            }
            require(i+1<argc,"Missing value for "+a);
            std::string v=argv[++i];
            if(a=="--data")o.inputs[0]=v;
            else if(a=="--before")o.inputs[1]=v;
            else if(a=="--after")o.inputs[2]=v;
            else if(a=="--data-before")o.dataBefore=v;
            else if(a=="--reference")o.reference=v;
            else if(a=="--output")o.output=v;
            else if(a=="--energy")o.energy=v;
            else if(a=="--low")o.low=std::stod(v);
            else if(a=="--high")o.high=std::stod(v);
            else if(a=="--min-bin")o.minBin=std::stod(v);
            else if(a=="--min-bins")o.minBins=std::stoi(v);
            else if(a=="--min-span")o.minSpan=std::stod(v);
            else if(a=="--min-hits")o.minHits=std::stod(v);
            else if(a=="--pivot")o.pivot=std::stod(v);
            else if(a=="--switch-margin")o.switchMargin=std::stod(v);
            else if(a=="--gain-switch")o.gainSwitch=std::stod(v);
            else if(a=="--top")o.top=std::stoi(v);
            else if(a=="--channels") {
                std::istringstream ss(v);
                std::string id;
                while(std::getline(ss,id,','))o.channels.push_back(std::stoi(id));
            }
            else throw std::runtime_error("Unknown option "+a);
        }
        for(auto&s:o.inputs)require(!s.empty(),"--data/--before/--after required");
        require(!o.dataBefore.empty()&&!o.reference.empty()&&!o.output.empty(),"--data-before/--reference/--output required");
        require(o.energy.find_first_not_of("0123456789.")==std::string::npos&&!o.energy.empty(),"Invalid energy label");
        std::set<int>ids;
        for(int id:o.channels) {
            int l=id/100000,c=id%100000/10000,ch=id%100;
            require(id>=0&&l<30&&c<6&&ch<(c==5?30:36)&&id==l*100000+c*10000+ch,"Invalid physical CellID");
            require(ids.insert(id).second,"Duplicate channel");
        }
        return o;
    }
}
// namespace hl
#ifndef HL_NO_MAIN
int main(int argc,char**argv) {
    try {
        hl::run(hl::parse(argc,argv));
        return 0;
    }
    catch(const std::exception&e) {
        std::cerr<<"ERROR: "<<e.what()<<std::endl;
        return 1;
    }
}
#endif
