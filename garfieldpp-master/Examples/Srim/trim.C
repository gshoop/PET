#include <iostream>

#include <TApplication.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackTrim.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  // Application
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Define the medium.
  MediumSilicon si;

  // Define the geometry.
  ComponentConstant cmp;
  cmp.SetElectricField(0., 1000., 0.);
  // Thickness of the silicon layer [cm]
  constexpr double d = 100.e-4;
  cmp.SetArea(-d, 0., -d, d, d, d);
  cmp.SetMedium(&si); 
  // cmp.SetMagneticField(0., 0., 1.);

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  
  // Create a track class.
  TrackTrim tr;
  // Connect the track to a sensor.
  tr.SetSensor(&sensor);
  // Read the TRIM output file.
  const std::string filename = "EXYZ.txt";
  // Import 100 ions, skip the first 200 in the list.
  const unsigned int nIons = 100;
  const unsigned int nSkip = 200;
  if (!tr.ReadFile(filename, nIons, nSkip)) {
    std::cerr << "Reading TRIM EXYZ file failed.\n";
    return 1;
  }
  tr.Print();

  // Plot the tracks.
  ViewDrift driftView;
  tr.EnablePlotting(&driftView);
  
  // Generate tracks.
  for (unsigned int i = 0; i < nIons; ++i) {
    if (!tr.NewTrack(0., 0., 0., 0., 0., 1., 0.)) {
      std::cerr << "Generating clusters failed; skipping this track.\n";
      continue;
    }
    unsigned int netot = 0; 
    // Retrieve the clusters.
    double xc, yc, zc, tc, ec, ekin;
    int ne = 0;
    while (tr.GetCluster(xc, yc, zc, tc, ne, ec, ekin)) {
      // Count the total number of electrons.
      netot += ne;
    }
  }
  driftView.SetArea(-10.e-4, 0., 10.e-4, 50.e-4);
  driftView.Plot(true);

  app.Run();
  return 0;
}

