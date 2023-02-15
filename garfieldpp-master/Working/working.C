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
    std::cout << "GET TEMPERATURE RETURNS: " << cdte.GetTemperature() << "\n";

    // const double ex = -1, ey = 0, ez = 0, bx = 0, by = 0, bz = 0;
    // double vx = 0, vy = 0, vz = 0;

    // cdte.ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
    
    // // Velocitites are in [cm/s]
    // std::cout << "ELECTRON VELOCITY: \n"
    //           << "\t vx:\t" << vx << "\n"
    //           << "\t vy:\t" << vy << "\n"
    //           << "\t vz:\t" << vz << "\n";

    constexpr bool plotVelocity = false;
    if (plotVelocity) {
        ViewMedium* mediumView = new ViewMedium();
        mediumView->SetRangeE(100.,1.e8,false);
        mediumView->SetMedium(&cdte);
        mediumView->PlotElectronVelocity('e');
        mediumView->PlotHoleVelocity('e',true);
    }

    constexpr double d = 0.45;
    constexpr double x = 0.1;
    constexpr double vbias = -500;
    ComponentConstant uniformField;
    uniformField.SetArea(-x,0.,x,x,d,x);
    uniformField.SetMedium(&cdte);
    uniformField.SetElectricField(0, vbias / d, 0 );
    uniformField.SetWeightingField(0, -1. / d, 0,"pad");

    auto eLinear = [](const double /*x*/, const double y, const double /*z*/,
                      double& ex, double& ey, double& ez) {
        // Depletion voltage [V]
        constexpr double vdep = -400.;
        ex = ez = 0.;
        ey = (vbias - vdep) / d + 2 * y * vdep / (d * d);
    };
    ComponentUser linearField;
    linearField.SetArea(-d, 0., -d, d, d, d);
    linearField.SetMedium(&cdte);
    linearField.SetElectricField(eLinear);

    // Making a component with analytic weighting field for a strip/pixel
    constexpr double pitch = 0.01;
    constexpr double halfpitch = 0.5 * pitch;
    ComponentAnalyticField wField;
    wField.SetMedium(&cdte);
    wField.AddPlaneY(0, vbias, "back");
    wField.AddPlaneY(d,0, "front");
    wField.AddStripOnPlaneY('z',d,-halfpitch,halfpitch,"strip");
    wField.AddPixelOnPlaneY(d,-halfpitch,halfpitch,
                            -halfpitch,halfpitch, "pixel");
    wField.AddReadout("strip");
    wField.AddReadout("pixel");
    wField.AddReadout("front");

    Sensor sensor;
    sensor.AddComponent(&linearField);
    const std::string label = "strip";
    sensor.AddElectrode(&wField, label);

    constexpr bool plotField =  true;
    if(plotField) {
      ViewField* fieldView = new ViewField();
      fieldView->SetSensor(&sensor);
      fieldView->SetArea(-x,0,x,d);
      fieldView->PlotContour("ey");
    }
    
    // Photoabsorption cross-section
    const double e = 511000.;
    double sigma = 12;
    cdte.GetPhotoAbsorptionCrossSection(e,sigma);

    std::cout << "PHOTOABSORPTION CROSS-SECTION: " << sigma << "\n";

    if (plotVelocity || plotField) {
    app.Run();
  }
}
