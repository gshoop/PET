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

    app.Run();
}