#pragma once
#include "RootInput.hh"
#include <iomanip>

struct CalibrationParameter {
    double ph=0,pl=0,mip=1,gain=1,offset=0;
    double inputMIP=NAN,chi2ndf=NAN,inputSlope=NAN,inputIntercept=NAN;
    bool pedPresent=false,mipSeen=false,hlSeen=false,mipAccepted=false,hlAccepted=false;
    bool mipFallback=false,gainFallback=false,offsetFallback=false;
};
inline int sipmGroup(int id) { return id/100000>=4 && id/100000<=27; }
// Audit the current *_simulation extraction policy used for BOTH data and MC.
// These are candidate constants, not proof of saved-file production provenance.
inline std::map<int,CalibrationParameter> readCalibration(
        const std::string& pedestal,const std::string& mip,const std::string& hl) {
    std::map<int,CalibrationParameter> result;
    for(int l=0;l<30;++l) for(int c=0;c<6;++c) for(int ch=0;ch<(c==5?30:36);++ch)
        result[l*100000+c*10000+ch]=CalibrationParameter{};
    result.at(0).offset=1; // C++ array={1}; Init() does not reset intercepts.
    for(const auto& [id,p]:readPedestal(pedestal)) if(result.count(id)) {
        auto& a=result.at(id); a.ph=static_cast<float>(p.hg); a.pl=static_cast<float>(p.lg); a.pedPresent=true;
    }
    double sm[2]={},sg[2]={},sb[2]={}; long long nm[2]={},nh[2]={};
    {
        auto f=openFile(mip); auto*t=tree(*f,"MIP_Fit"); int id=0,ndf=0;double mpv=0,chi=0;
        bind(t,"CellID",&id);bind(t,"LandauMPV",&mpv);bind(t,"ChiSquare",&chi);bind(t,"NDF",&ndf);
        std::set<int> seen;
        for(Long64_t i=0;i<t->GetEntries();++i) {
            if(t->GetEntry(i)<=0) throw std::runtime_error("Cannot read MIP entry");
            int p=physicalID(id); if(p/100000>=30) continue;
            if(!seen.insert(p).second) throw std::runtime_error("Duplicate MIP channel");
            if(!std::isfinite(mpv)||!std::isfinite(chi)) throw std::runtime_error("Nonfinite MIP fit");
            int g=sipmGroup(p); double q=chi/static_cast<double>(ndf);
            bool accepted=!(q>(g?2.:1.7));
            if(accepted) { sm[g]+=mpv; ++nm[g]; } // Includes nonpositive MPVs, as production does.
            if(result.count(p)) {
                auto&a=result.at(p);a.mipSeen=true;a.inputMIP=mpv;a.chi2ndf=q;a.mipAccepted=accepted;
                if(accepted) a.mip=mpv;
            }
        }
        t->ResetBranchAddresses();
    }
    {
        auto f=openFile(hl);auto*t=tree(*f,"InterCalib");int id=0;double slope=0,intercept=0;
        bind(t,"CellID",&id);bind(t,"Slope",&slope);bind(t,"Intercept",&intercept);
        std::set<int> seen;
        for(Long64_t i=0;i<t->GetEntries();++i) {
            if(t->GetEntry(i)<=0) throw std::runtime_error("Cannot read HL entry");
            int p=physicalID(id);if(p/100000>=30)continue;
            if(!seen.insert(p).second)throw std::runtime_error("Duplicate HL channel");
            if(!std::isfinite(slope)||!std::isfinite(intercept))throw std::runtime_error("Nonfinite HL fit");
            bool accepted=slope>=.02&&slope<=.05;int g=sipmGroup(p);
            if(accepted){sg[g]+=1/slope;sb[g]+=-intercept/slope;++nh[g];}
            if(result.count(p)) {
                auto&a=result.at(p);a.hlSeen=true;a.inputSlope=slope;a.inputIntercept=intercept;a.hlAccepted=accepted;
                if(accepted){a.gain=1/slope;a.offset=-intercept/slope;}
            }
        }
        t->ResetBranchAddresses();
    }
    for(int g=0;g<2;++g)if(!nm[g]||!nh[g])throw std::runtime_error("No accepted constants in a SiPM group");
    for(auto& [id,a]:result) {
        int g=sipmGroup(id);a.mipFallback=a.mip==1;a.gainFallback=a.gain==1;a.offsetFallback=a.offset==1;
        if(a.mipFallback)a.mip=sm[g]/nm[g];
        if(a.gainFallback)a.gain=sg[g]/nh[g];
        if(a.offsetFallback)a.offset=sb[g]/nh[g];
    }
    return result;
}
inline void writeCalibration(const std::map<int,CalibrationParameter>& p,const std::string& path) {
    std::ofstream out(path);if(!out)throw std::runtime_error("Cannot write "+path);
    out<<std::setprecision(17)
       <<"cellid\tped_hg_ADC\tped_lg_ADC\tped_present\tped_policy\tmip_input_ADC\tchi2_ndf\tmip_adopted_ADC\tmip_status\tmip_fallback\tmip_nonpositive\thl_input_slope\thl_input_intercept_LG_ADC\thl_gain\thl_offset_HG_ADC\thl_status\thl_gain_fallback\thl_offset_fallback\thl_offset_unset_zero\n";
    for(const auto& [id,a]:p)out<<id<<'\t'<<a.ph<<'\t'<<a.pl<<'\t'<<a.pedPresent<<'\t'
        <<(a.pedPresent?"file_float_rounding":"missing_production_zero")<<'\t'<<a.inputMIP<<'\t'<<a.chi2ndf<<'\t'<<a.mip<<'\t'
        <<(!a.mipSeen?"missing":!a.mipAccepted?"quality_rejected":a.mipFallback?"sentinel_1":"accepted")<<'\t'
        <<a.mipFallback<<'\t'<<(a.mip<=0)<<'\t'<<a.inputSlope<<'\t'<<a.inputIntercept<<'\t'<<a.gain<<'\t'<<a.offset<<'\t'
        <<(!a.hlSeen?"missing":!a.hlAccepted?"slope_rejected":"accepted")<<'\t'<<a.gainFallback<<'\t'<<a.offsetFallback<<'\t'
        <<(!a.hlAccepted&&!a.offsetFallback)<<'\n';
    if(!out)throw std::runtime_error("Failed writing calibration audit");
}
