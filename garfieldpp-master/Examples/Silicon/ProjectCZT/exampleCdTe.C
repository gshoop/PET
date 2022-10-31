#include <iostream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1D.h>

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
#include "Garfield/MediumCdTe.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

    TApplication app("app", &argc, argv);

    MediumCdTe cdte;

    ViewMedium mediumView;
    mediumView.SetMedium(&cdte);
    constexpr bool plotVelocity = true;
    if (plotVelocity) {
        mediumView.PlotElectronVelocity('e');
        mediumView.PlotHoleVelocity('e',true);
    }

    constexpr double d = 100e-4;
    // Bias voltage [V]
    constexpr double vbias = -50.;

    ComponentConstant uniformField;
    uniformField.SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
    uniformField.SetMedium(&cdte);
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
    linearField.SetMedium(&cdte);


    constexpr double pitch = 50.e-4;
    constexpr double halfpitch = 0.5 * pitch;
    ComponentAnalyticField wField;
    wField.SetMedium(&cdte);
    wField.AddPlaneY(0, vbias, "back");
    wField.AddPlaneY(d, 0, "front");
    wField.AddStripOnPlaneY('z', d, -halfpitch, halfpitch, "strip");
    wField.AddReadout("strip");
    wField.AddPixelOnPlaneY(d, -halfpitch, halfpitch, -halfpitch, halfpitch, "pixel");
    wField.AddReadout("pixel");

    Sensor sensor;
    const std::string label = "strip";
    sensor.AddComponent(&linearField);
    sensor.AddElectrode(&wField, "strip");
    
    if (plotVelocity) {
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
    }

    const unsigned int nTimeBins = 1000;
    const double tmin = 0.;
    const double tmax = 10.;
    const double tstep = (tmax - tmin) / nTimeBins;
    sensor.SetTimeWindow(tmin, tstep, nTimeBins);

    TrackHeed track;
    track.SetSensor(&sensor);
    // Set the particle type and momentum [ev/c]
    track.SetParticle("pion");
    track.SetMomentum(180.e9);

    AvalancheMC drift;
    drift.SetSensor(&sensor);
    // Set the step size [cm]
    drift.SetDistanceSteps(1.e-4);
    drift.EnableSignalCalculation();

    constexpr bool plotSignal = true;
    ViewSignal* signalView = nullptr;
    TCanvas* cSignal = nullptr;
    if (plotSignal) { 
        cSignal = new TCanvas("cSignal", "", 600, 600);
        signalView = new ViewSignal();
        signalView->SetCanvas(cSignal);
        signalView->SetSensor(&sensor);
    }

    constexpr bool plotDrift = true;
    ViewDrift* driftView = nullptr;
    TCanvas* cDrift = nullptr;
    if (plotDrift) {
        cDrift = new TCanvas("cDrift", "", 600, 600);
        driftView = new ViewDrift();
        driftView->SetArea(-0.5 * d, 0, -0.5 * d, 0.5 * d, d, 0.5 * d);
        driftView->SetCanvas(cDrift);
        track.EnablePlotting(driftView);
    }

    constexpr bool smearx = true;
    double x0 = 0., y0 = 0., z0 = 0., t0 = 0.;
    double dx = 0., dy = 1., dz = 0.;
    constexpr unsigned int nEvents = 10;

    double xt = 0.;
    if (smearx) xt = -0.5 * pitch + RndmUniform() * pitch;

    track.NewTrack(xt, y0, z0, t0, dx, dy, dz);
    double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
    int ne = 0;
    //Retrieve the clusters along the track.
    while (track.GetCluster(xc, yc, zc, tc, ne, ec, extra)) {
        // Loop over the electrons in the cluster.
        for (int j = 0; j < ne; ++j) {
            double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
            double dxe = 0., dye = 0., dze = 0.;
            track.GetElectron(j, xe, ye, ze, te, ee, dxe, dye, dze);
            // Simulate the electron and hole drift lines.
            drift.DisablePlotting();
           // if (RndmUniform() < 0.01) drift.EnablePlotting(driftView);
            drift.EnablePlotting(driftView);
            drift.DriftElectron(xe, ye, ze, te);
            drift.DriftHole(xe, ye, ze, te);
        }
    }
    constexpr bool twod = true;
    driftView->Plot(twod);
    cDrift->Update();
    gSystem->ProcessEvents();

    signalView->PlotSignal(label);
    cSignal->Update();
    gSystem->ProcessEvents();

    app.Run();
}