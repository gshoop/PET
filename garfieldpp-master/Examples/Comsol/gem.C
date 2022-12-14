#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "Garfield/ComponentComsol.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Load the field map.
  ComponentComsol fm;
  fm.Initialise("mesh.mphtxt", "dielectrics.dat", "field.txt", "mm");
  fm.EnableMirrorPeriodicityX();
  fm.EnableMirrorPeriodicityY();
  fm.PrintRange();

  // Dimensions of the GEM [cm]
  constexpr double pitch = 0.014;

  ViewField fieldView;
  constexpr bool plotField = true;
  if (plotField) {
    fieldView.SetComponent(&fm);
    // Set the normal vector of the viewing plane (xz plane).
    fieldView.SetPlane(0, -1, 0, 0, 0, 0);
    // Set the plot limits in the current viewing plane.
    fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    fieldView.SetVoltageRange(-160., 160.);
    TCanvas* cf = new TCanvas("cf", "", 600, 600);
    cf->SetLeftMargin(0.16);
    fieldView.SetCanvas(cf);
    fieldView.PlotContour();
  }

  // Setup the gas.
  MediumMagboltz gas;
  gas.SetComposition("ar", 80., "co2", 20.);
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  gas.Initialise(true);  
  // Set the Penning transfer efficiency.
  constexpr double rPenning = 0.51;
  constexpr double lambdaPenning = 0.;
  gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // Load the ion mobilities.
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
  // Associate the gas with the corresponding field map material.
  fm.SetGas(&gas); 
  fm.PrintMaterials();
 
  // Create the sensor.
  Sensor sensor;
  sensor.AddComponent(&fm);
  sensor.SetArea(-5 * pitch, -5 * pitch, -0.01,
                  5 * pitch,  5 * pitch,  0.025);

  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);

  AvalancheMC drift;
  drift.SetSensor(&sensor);
  drift.SetDistanceSteps(2.e-4);

  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    aval.EnablePlotting(&driftView);
    drift.EnablePlotting(&driftView);
  }

  constexpr unsigned int nEvents = 10;
  for (unsigned int i = 0; i < nEvents; ++i) { 
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position. 
    const double x0 = -0.5 * pitch + RndmUniform() * pitch;
    const double y0 = -0.5 * pitch + RndmUniform() * pitch;
    const double z0 = 0.02; 
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    const unsigned int np = aval.GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (unsigned int j = 0; j < np; ++j) {
      aval.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                  xe2, ye2, ze2, te2, e2, status);
      drift.DriftIon(xe1, ye1, ze1, te1);
      drift.GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
                              xi2, yi2, zi2, ti2, status);
    }
  }
  if (plotDrift) {
    TCanvas* cd = new TCanvas();
    constexpr bool plotMesh = true;
    if (plotMesh) {
      ViewFEMesh* meshView = new ViewFEMesh();
      meshView->SetArea(-2 * pitch, -2 * pitch, -0.02, 
                         2 * pitch,  2 * pitch, 0.02);
      meshView->SetCanvas(cd);
      meshView->SetComponent(&fm);
      // x-z projection.
      meshView->SetPlane(0, -1, 0, 0, 0, 0);
      meshView->SetFillMesh(true);
      // Set the color of the kapton and the metal.
      meshView->SetColor(1, kYellow + 3);
      meshView->SetColor(2, kGray);
      meshView->EnableAxes();
      meshView->SetViewDrift(&driftView);
      meshView->Plot();
    } else {
      driftView.SetPlane(0, -1, 0, 0, 0, 0);
      driftView.SetArea(-2 * pitch, -0.02, 2 * pitch, 0.02);
      driftView.SetCanvas(cd);
      constexpr bool twod = true;
      driftView.Plot(twod);
    }
  }

  app.Run(kTRUE);
}
