#include "GeoMagLib.h"

#pragma comment(lib, "GeoMagLib.lib")

CGeoMagnetic* GeoMagnetic = nullptr;

double sdate = 2015.15;

char projection = 'M';

double alt = 0.0;

double longitude = -81.440772;
double latitude = 39.342776;

char unit = 'F';

int main()
{
    GeoMagnetic = new CGeoMagnetic(sdate, projection);

    for (int i = 0; i < 20 * 500; i+=500)
    {
        alt += i;

        GeoMagnetic->CalculateFieldElements(latitude, longitude, alt, unit);

        printf("\nAlt:%9.2f\tDecl:%f", alt, GeoMagnetic->GeoMagneticElements->Decl);
    }

    delete GeoMagnetic;

    return 0;
}
