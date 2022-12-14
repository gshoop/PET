#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  MediumMagboltz gas("ne", 90., "co2", 10.);

  // Define the cell layout.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  cmp.AddPlaneX(-0.5,    0., "p1");
  cmp.AddPlaneX( 0.5, 1000., "p2");
  // Add a readout strip along z.
  cmp.AddStripOnPlaneX('z', -0.5, -0.1, 0.1, "strip"); 
  cmp.AddReadout("strip");
  cmp.PrintCell();

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "strip");

  // Plot the weighting potential.
  const double xmin = -0.6;
  const double xmax =  0.6;
  const double ymin = -0.6;
  const double ymax =  0.6;
  TCanvas canvas("c", "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetSensor(&sensor);
  fieldView.SetArea(xmin, ymin, xmax, ymax);
  fieldView.PlotContourWeightingField("strip", "v");

  ViewCell cellView;
  cellView.SetCanvas(&canvas);
  cellView.SetComponent(&cmp);
  cellView.SetArea(xmin, ymin, xmax, ymax);
  cellView.Plot2d();
  gSystem->ProcessEvents();
  app.Run(true);
}
