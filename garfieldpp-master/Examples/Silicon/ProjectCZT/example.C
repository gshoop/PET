#include <iostream>

#include <TApplication.h>

#include "Garfield/SolidBox.hh" 
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/Sensor.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

    TApplication app("app", &argc, argv);

    // Create Silicon Medium and Set temperature
    MediumSilicon si;
    si.SetTemperature(293.15);

    ViewMedium mediumView;
    mediumView.SetMedium(&si);
    mediumView.PlotElectronVelocity('e');
    const bool same = true;
    mediumView.PlotHoleVelocity('e',same);

    // sensor thickness [cm]
    constexpr double d = 100.e-4;
    // Bias voltage [V]
    constexpr double vbias = -50.;

    ComponentConstant uniformField;
    uniformField.SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
    uniformField.SetMedium(&si);
    uniformField.SetElectricField(0, vbias / d, 0);
    uniformField.SetWeightingField(0, -1. / d, 0, "pad");

    auto eLinear = [](const double /*x*/, const double y, const double /*y*/, double& ex, double& ey, double& ez) {
        // Depletion voltage [V]
        constexpr double vdep = -20;
        ex = ez = 0.;
        ey = (vbias - vdep) / d + 2 * y * vdep / (d * d);   // Not sure where this expression comes from
    };

    ComponentUser linearField;
    linearField.SetElectricField(eLinear);
    linearField.SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
    linearField.SetMedium(&si);

    Sensor sensor;
    sensor.AddComponent(&linearField);

    constexpr double pitch = 55.e-4;
    constexpr double halfpitch = 0.5 * pitch;
    ComponentAnalyticField wField;
    wField.SetMedium(&si);
    wField.AddPlaneY(0, vbias, "back");
    wField.AddPlaneY(d, 0, "front");
    wField.AddStripOnPlaneY('z', d, -halfpitch, halfpitch, "strip");
    wField.AddReadout("strip");
    wField.AddPixelOnPlaneY(d, -halfpitch, halfpitch, -halfpitch, halfpitch, "pixel");
    wField.AddReadout("pixel");

    ViewField* wfieldView = new ViewField();
    wfieldView->SetComponent(&wField);
    wfieldView->SetArea(-0.5 * d, 0, 0.5 * d, d);
    wfieldView->PlotContourWeightingField("strip", "v");

    ViewField* wfieldViewA = new ViewField();
    wfieldViewA->SetComponent(&wField);
    wfieldViewA->SetArea(-0.5 * d, 0, 0.5 * d, d);
    wfieldViewA->PlotContourWeightingField("strip", "ex");

    ViewField* wfieldViewB = new ViewField();
    wfieldViewB->SetSensor(&sensor);
    wfieldViewB->SetArea(-0.5 * d, 0, 0.5 * d, d);
    wfieldViewB->PlotContour("ey");

    // sensor.AddElectrode(&wField, "strip");

    // const unsigned int nTimeBins = 1000;
    // const double tmin = 0.;
    // const double tmax = 10.;
    // const double tstep = (tmax - tmin) / nTimeBins;
    // sensor.SetTimeWindow(tmin, tstep, nTimeBins);

    app.Run();
}