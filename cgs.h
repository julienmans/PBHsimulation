#ifndef INCLUDE_CONSTANTS_H
#define INCLUDE_CONSTANTS_H

#include <string>

// CGS units
static const double centimeter = 1;
static const double cm = centimeter;
static const double gram = 1;
static const double second = 1;
static const double s = second;
static const double dyne = 1;
static const double erg = 1;
static const double esu = 1;
static const double gauss = 1;
static const double kelvin = 1;
static const double K = kelvin;

// derived units
static const double meter = 1e2 * centimeter;
static const double kilometer = 1e3 * meter;
static const double km = kilometer;
static const double cm2 = cm * cm;
static const double cm3 = cm * cm * cm;
static const double kilogram = 1e3 * gram;
static const double newton = 1e5 * dyne;
static const double joule = 1e7 * erg;
static const double coulomb = 2997924580. * esu;
static const double tesla = 1e4 * gauss;
static const double microgauss = 1e-6 * gauss;
static const double nanogauss = 1e-9 * gauss;
static const double muG = microgauss;
static const double nG = nanogauss;

// physical constants
static const double c_light = 2.99792458e10 * centimeter / second;
static const double c_squared = c_light * c_light;
static const double mass_proton = 1.67262158e-24 * gram;
static const double mass_neutron = 1.67492735e-24 * gram;
static const double mass_electron = 9.10938291e-28 * gram;
static const double mass_sun = 1.989e33 * gram;
static const double h_planck = 6.62606957e-27 * erg * second;
static const double k_boltzmann = 1.3806488e-16 * erg / kelvin;
static const double G_Newton = 6.693 * 1.e-11 *pow(meter,3.)*(1./kilogram)*(1./(second*second));

// electron volt
static const double electronvolt = 1.60217657e-12 * erg;
static const double kiloelectronvolt = 1e3 * electronvolt;
static const double megaelectronvolt = 1e6 * electronvolt;
static const double gigaelectronvolt = 1e9 * electronvolt;
static const double teraelectronvolt = 1e12 * electronvolt;
static const double petaelectronvolt = 1e15 * electronvolt;
static const double exaelectronvolt = 1e18 * electronvolt;
static const double eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double MeV = megaelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;
static const double EeV = exaelectronvolt;

// time
static const double year = 3.15569e7 * second;
static const double kiloyear = 1e3 * year;
static const double megayear = 1e6 * year;
static const double gigayear = 1e9 * year;
static const double kyr = kiloyear;
static const double Myr = megayear;
static const double Gyr = gigayear;

// parsec
static const double parsec = 3.0856775807e18 * centimeter;
static const double kiloparsec = 1e3 * parsec;
static const double megaparsec = 1e6 * parsec;
static const double gigaparsec = 1e9 * parsec;
static const double pc = parsec;
static const double pc3 = pc*pc*pc;
static const double kpc = kiloparsec;
static const double Mpc = megaparsec;
static const double Gpc = gigaparsec;

/*const double cYear = 365.25*24.*3600.;
const double cKpc = 3.08568e21;
const double cKpc2pc = 1e3;
const double cKm = 1e5;
const double cLight = 2.99792458e10; // [cm/s]
const double hPlanck = 1; // TBD
const double kBoltzmann = 1; // TBD
const double cProtonMassInGeV = 0.938;
const double cElectronMassInMeV = 0.5109990615;
const double cElectronMassInGeV = 1e-3*cElectronMassInMeV;
const double cMinNucleiLifetime = 1e3; // [seconds] Lifetime below which a nucleus is not considered
const double cPi = 3.14;
*/

const double cXSun = 8.5 * kpc;
const double cYSun = 0.0 * kpc;
const double cZSun = 0.0 * kpc;

const  double degreeToRad = M_PI / 180.0;
const  double rad = 1.;
const  double srad = 1.;
const  double deg = 0.0174533;
const  double arcmin = 0.000290888;
const double squareDeg = 3.0462e-4*srad;

const std::string isotopeListFilename = "data/particleList.txt";
//const std::string sourceParamsDefaultFilename = "data/default.source";

#endif
