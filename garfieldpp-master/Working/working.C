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

    // TApplication app("app", &argc, argv);

    MediumCdTe cdte;
    std::cout << "GET TEMPERATURE RETURNS: " << cdte.GetTemperature() << "\n";

    const double ex = -1, ey = 0, ez = 0, bx = 0, by = 0, bz = 0;
    double vx = 0, vy = 0, vz = 0;

    cdte.ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
    
    // Velocitites are in [cm/s]
    std::cout << "ELECTRON VELOCITY: \n"
              << "\t vx:\t" << vx << "\n"
              << "\t vy:\t" << vy << "\n"
              << "\t vz:\t" << vz << "\n";

    // Photoabsorption cross-section
    const double e = 511000.;
    double sigma = 12;
    cdte.GetPhotoAbsorptionCrossSection(e,sigma);

    std::cout << "PHOTOABSORPTION CROSS-SECTION: " << sigma << "\n";

    // app.Run();
}
