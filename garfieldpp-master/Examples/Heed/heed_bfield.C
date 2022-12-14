#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  MediumMagboltz gas("ar", 90., "co2", 10.);

  ComponentConstant cmp;
  cmp.SetArea(-10., -10., -10., 10., 10., 10.);
  cmp.SetMedium(&gas);
  cmp.SetElectricField(100., 0., 0.);
  cmp.SetMagneticField(0., 1., 0.);

  Sensor sensor;
  sensor.AddComponent(&cmp);

  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle("muon");
  track.SetMomentum(1.e9);
  ViewDrift view;
  track.EnablePlotting(&view);

  ///*
  double maxrange = 0., rforstraight = 0., stepstraight = 0., stepcurved = 0.;
  track.GetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);
  stepcurved = 0.01;
  track.SetSteppingLimits(maxrange, rforstraight, stepstraight, stepcurved);
  //*/

  constexpr double z0 = -10.;
  track.NewTrack(0., 0., z0, 0., 0., 0., 1.);
  double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
  int nc = 0;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {

  }
  view.Plot();
  app.Run();
  return 0;
}
