#include <TMath.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
 
#include <cassert>
 
using namespace ROOT::Math;
 
 
// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
   // a parametric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}
 
 
bool first = true;
 
// function Object to be minimized
struct SumDistance2 {
   // the TGraph is a data member of the object
   TGraph2DErrors *fGraph;
 
   SumDistance2(TGraph2DErrors *g) : fGraph(g) {}
 
   // calculate distance line-point
   double distance2(double x,double y,double z, const double *p) {
      // distance line point is D= | (xp-x0) cross  ux |
      // where ux is direction of line and x0 is a point in the line (like t = 0)
      XYZVector xp(x,y,z);
      XYZVector x0(p[0], p[2], 0. );
      XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
      XYZVector u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      return d2;
   }
 
   // implementation of the function to be minimized
   double operator() (const double *par) {
      assert(fGraph != 0);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i  = 0; i < npoints; ++i) {
         double d = distance2(x[i],y[i],z[i],par);
         sum += d;
      }
      if (first) {
         std::cout << "Total Initial distance square = " << sum << std::endl;
      }
      first = false;
      return sum;
   }
 
};
 
double* line3Dfit(vector<double> x_arr,vector<double> y_arr, vector<double> z_arr, TGraph2DErrors* gr)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit();
 
 
   //double e = 0.1;
   int nd = x_arr.size();
 
 
   // double xmin = 0; double ymin = 0;
   // double xmax = 10; double ymax = 10;
 
   // TGraph2D * gr = new TGraph2D();
 
   // Fill the 2D graph
   // double p0[4] = {10,20,1,2};
 
   // generate graph with the 3d points
   // set point in the graph
   for (int N=0; N<nd; N++) {
      double x,y,z = 0;
      // // Generate a random number
      // double t = gRandom->Uniform(0,10);
      // line(t,p0,x,y,z);
      // double err = 2.5;
      // // do a Gaussian smearing around the points in all coordinates
      x = x_arr[N];
      y = y_arr[N];
      z = z_arr[N];
      gr->SetPoint(N,x,y,z);
      gr->SetPointError(N,2.5,2.5,1);
   }
   // fit the graph now
 
   ROOT::Fit::Fitter  fitter;
 
   // make the functor objet
   SumDistance2 sdist(gr);
   ROOT::Math::Functor fcn(sdist,4);
   // set the function and the initial parameter values
   double pStart[4] = {1,1,1,1};
   fitter.SetFCN(fcn,pStart);
   // set step sizes different than default ones (0.3 times parameter values)
   for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
 
   bool ok = fitter.FitFCN(fcn,pStart,nd,true);
   // bool ok = fitter.FitFCN();
   // if (!ok) {
   //    Error("line3Dfit","Line3D Fit failed");
   //    return 1;
   // }

   const ROOT::Fit::FitResult & result = fitter.Result();
 
   // std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
   // result.Print(std::cout);
 
 
   // gr->Draw("p0");
 
   // get fit parameters
   const double * parFit = result.GetParams();
   const double * parFitError = result.GetErrors();
   const double chi2 = result.Chi2();
   const double ndf = result.Ndf();

   // static double par[] = {parFit[1], parFitError[1], parFit[0], parFitError[0], parFit[3],parFitError[3], parFit[2], parFitError[2]};
   // static double par[8] = {0};
   static double par[10] = {0};
   par[0] = parFit[1];
   par[1] = parFitError[1];
   par[2] = parFit[0];
   par[3] = parFitError[0];
   par[4] = parFit[3];
   par[5] = parFitError[3];
   par[6] = parFit[2];
   par[7] = parFitError[2];
   par[8] = chi2;
   par[9] = ndf;

   // cout << parFit[1]<< " " << parFitError[1]<< " " << parFit[0]<< " " << parFitError[0]<< " " << parFit[3]<< " " << parFitError[3]<< " " << parFit[2]<< " " << parFitError[2] << endl;

   return par;

   // // draw the fitted line
   // int n = 1000;
   // double t0 = 0;
   // double dt = 10;
   // TPolyLine3D *l = new TPolyLine3D(n);
   // for (int i = 0; i <n;++i) {
   //    double t = t0+ dt*i/n;
   //    double x,y,z;
   //    line(t,parFit,x,y,z);
   //    l->SetPoint(i,x,y,z);
   // }
   // l->SetLineColor(kRed);
   // l->Draw("same");
 
   // // draw original line
   // TPolyLine3D *l0 = new TPolyLine3D(n);
   // for (int i = 0; i <n;++i) {
   //    double t = t0+ dt*i/n;
   //    double x,y,z;
   //    line(t,p0,x,y,z);
   //    l0->SetPoint(i,x,y,z);
   // }
   // l0->SetLineColor(kBlue);
   // l0->Draw("same");
   // return 0;
}

double* line3Dfit(vector<double> x_arr,vector<double> y_arr, vector<double> z_arr)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit();
 
 
   //double e = 0.1;
   int nd = x_arr.size();
 
 
   // double xmin = 0; double ymin = 0;
   // double xmax = 10; double ymax = 10;
 
   TGraph2DErrors * gr = new TGraph2DErrors();
 
   // Fill the 2D graph
   // double p0[4] = {10,20,1,2};
 
   // generate graph with the 3d points
   // set point in the graph
   for (int N=0; N<nd; N++) {
      double x,y,z = 0;
      // // Generate a random number
      // double t = gRandom->Uniform(0,10);
      // line(t,p0,x,y,z);
      // double err = 2.5;
      // // do a Gaussian smearing around the points in all coordinates
      x = x_arr[N];
      y = y_arr[N];
      z = z_arr[N];
      gr->SetPoint(N,x,y,z);
      gr->SetPointError(N,2.5,2.5,1);
      //dt->SetPointError(N,0,0,err);
   }
   // fit the graph now
 
   ROOT::Fit::Fitter  fitter;
 
   // make the functor objet
   SumDistance2 sdist(gr);
   ROOT::Math::Functor fcn(sdist,4);
   // set the function and the initial parameter values
   double pStart[4] = {1,1,1,1};
   fitter.SetFCN(fcn,pStart);
   // set step sizes different than default ones (0.3 times parameter values)
   for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
 
   bool ok = fitter.FitFCN(fcn,pStart,nd,true);
   // bool ok = fitter.FitFCN();
   // if (!ok) {
   //    Error("line3Dfit","Line3D Fit failed");
   //    return 1;
   // }

   const ROOT::Fit::FitResult & result = fitter.Result();
 
   // std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
   // result.Print(std::cout);
 
 
   // gr->Draw("p0");
 
   // get fit parameters
   const double * parFit = result.GetParams();
   const double * parFitError = result.GetErrors();
   const double chi2 = result.Chi2();
   const double ndf = result.Ndf();

   // static double par[8] = {0};
   static double par[10] = {0};
   par[0] = parFit[1];
   par[1] = parFitError[1];
   par[2] = parFit[0];
   par[3] = parFitError[0];
   par[4] = parFit[3];
   par[5] = parFitError[3];
   par[6] = parFit[2];
   par[7] = parFitError[2];
   par[8] = chi2;
   par[9] = ndf;
 
   return par;

   // // draw the fitted line
   // int n = 1000;
   // double t0 = 0;
   // double dt = 10;
   // TPolyLine3D *l = new TPolyLine3D(n);
   // for (int i = 0; i <n;++i) {
   //    double t = t0+ dt*i/n;
   //    double x,y,z;
   //    line(t,parFit,x,y,z);
   //    l->SetPoint(i,x,y,z);
   // }
   // l->SetLineColor(kRed);
   // l->Draw("same");
 
   // // draw original line
   // TPolyLine3D *l0 = new TPolyLine3D(n);
   // for (int i = 0; i <n;++i) {
   //    double t = t0+ dt*i/n;
   //    double x,y,z;
   //    line(t,p0,x,y,z);
   //    l0->SetPoint(i,x,y,z);
   // }
   // l0->SetLineColor(kBlue);
   // l0->Draw("same");
   // return 0;
}
