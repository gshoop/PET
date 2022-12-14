#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TGraph.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Make a medium
  MediumSilicon si;
  const double rho = si.GetMassDensity();
  std::cout << "Density: " << rho << std::endl; 
  // Make a drift volume
  constexpr double length = 100.;
  GeometrySimple geo;
  SolidBox box(0, 0, length, length, length, length);
  geo.AddSolid(&box, &si);
  
  // Make a component with constant drift field
  ComponentConstant cmp;
  cmp.SetGeometry(&geo);
  constexpr double field = 10.;
  cmp.SetElectricField(0., 0., field);

  // Make a sensor
  Sensor sensor;
  sensor.AddComponent(&cmp);
  
  // Heed
  TrackHeed track;
  track.SetSensor(&sensor);
  
  // Histograms
  TH1::StatOverflows();
  constexpr int nBins = 1000;
  double zMax = 50.;
  TH1F hLong("hLong", ";#it{z} [#mum];", nBins, 0., zMax);
  hLong.SetLineColor(kBlue + 2);
  hLong.SetLineWidth(2);
  hLong.SetStats(false);
  TGraph gRange95;
  gRange95.SetMarkerStyle(20);
  gRange95.SetMarkerSize(1);
  gRange95.SetMarkerColor(kBlue + 2);
  TCanvas c1;
  
  const unsigned int nEvents = 1e6;
  double de = 1000.;
  constexpr double emin = 1000.;
  constexpr double emax = 5000.;
  double e0 = emin;
  while (e0 <= emax) {
    // Estimate the range [mg / cm2].
    double rp = 69.7e-3 * pow(e0 * 1.e-3, 1.6);
    rp /= rho;
    zMax = 2 * rp;
    if (e0 < 500.) zMax += rp;
    if (e0 < 100.) zMax += 3 * rp;
    hLong.SetBins(nBins, 0., zMax);
    hLong.Reset("ICE");
    hLong.ResetStats();
    std::cout << "Primary energy: " << e0 << " eV" << std::endl;
    double nEntries = 0.;
    for (unsigned int i = 0; i < nEvents; ++i) {
      int np = 0;
      track.TransportDeltaElectron(0, 0, 0, 0, e0, 0, 0, 1, np);
      double x1, y1, z1, t1, e1, dx1, dy1, dz1;
      if (np <= 1) continue;
      for (int j = np - 1; j--;) {
        track.GetElectron(j, x1, y1, z1, t1, e1, dx1, dy1, dz1);
        if (fabs(z1) < 1.e-8) continue;
        hLong.Fill(z1 * 1.e4);
        nEntries += 1.; 
      }
    }

    constexpr double fraction = 0.95;

    double sum = 0.;
    int iUp = 0;
    for (int j = 1; j <= nBins; ++j) {
      sum += hLong.GetBinContent(j);
      if (sum >= fraction * nEntries) {
        iUp = j; 
        break;
      }
    }
    if (iUp <= 1) {
      std::cerr << "    Unable to determine the range.\n";
      // continue;
    }
    int iLow = nBins;
    sum = 0.;
    for (int j = nBins; j--;) {
      sum += hLong.GetBinContent(j);
      if (sum >= (1. - fraction) * nEntries) {
        iLow = j;
        break;
      }
    }
    const double r95 = 0.5 * (hLong.GetBinCenter(iUp) +
                              hLong.GetBinCenter(iLow));
    const double r95rho = r95 * rho * 1.e2;
    std::cout << "    Range: " << r95rho << " ug/cm2\n";
    const double y95 = 0.5 * (hLong.GetBinContent(iUp) +
                              hLong.GetBinContent(iLow));
    c1.cd();
    c1.Clear();
    hLong.Draw();
    gRange95.DrawGraph(1, &r95, &y95, "psame");
    c1.Update();
    gSystem->ProcessEvents();

    std::ofstream outfile;  
    outfile.open("r95_Heed_Si.txt", std::ios::out | std::ios::app);
    outfile << e0 << "  " << r95 << "  " << r95rho << "\n";
    outfile.close();

    e0 += de;
  }

  app.Run(true);

}
