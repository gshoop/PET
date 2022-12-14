#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/Random.hh"

TH1F* hZ = nullptr;

void userHandle(double x, double y, double z, double /*t*/,
                int type, int level, Garfield::Medium* /*m*/) {

  // Skip inelastic collisions that are not excitations.
  if (type != 4) return;

  // Do something with x, y, z.
  // For instance, fill a histogram of the z coordinates.
  hZ->Fill(z);
}

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  // Make a gas medium.
  MediumMagboltz gas("xe");
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  if (!gas.Initialise()) return 1;

  // Make a component with constant drift field.
  constexpr double gap = 0.7;
  ComponentConstant cmp;
  cmp.SetArea(-2, -2, -gap, 2, 2, gap);
  cmp.SetMedium(&gas);
  cmp.SetElectricField(0, 0, -8000);

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);

  // Microscopic tracking.
  AvalancheMicroscopic aval;
  aval.SetSensor(&sensor);
  aval.SetUserHandleInelastic(userHandle);
  
  // Histogram of the z coordinates of the excitations. 
  hZ = new TH1F("hZ", "", 50, 0., gap);
 
  const unsigned int nEvents = 10;
  for (unsigned int i = 0; i < nEvents; ++i) {
    constexpr double e0 = 1.;
    std::cout << i << "/" << nEvents << "\n";
    aval.AvalancheElectron(0, 0, 0, 0, e0, 0, 0, 0);
  }
  TCanvas cZ("cZ", "", 600, 600);
  cZ.cd();
  hZ->Draw();
  hZ->GetXaxis()->SetTitle("z [cm]");
  cZ.Update();

  app.Run(true);
}
