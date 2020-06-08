
#ifndef I3CORSIKAREADER_UTILS_H_INCLUDED
#define I3CORSIKAREADER_UTILS_H_INCLUDED

#include <cstdint>
#include <icetray/I3Units.h>
#include <dataclasses/I3Constants.h>

class I3Direction;
class I3Position;

namespace I3CORSIKAReaderUtils {

int32_t CorsikaToPDG(int corsika_id);
int32_t PDGToCorsika(int32_t pdg_id);

const double EarthRadius = 637131500*I3Units::cm; // as in CORSIKA

double GetSlantDepth(const I3Direction &dir, const I3Position &pos, double altitude=I3Constants::SurfaceElev);

double
LocalZenith(double reference_zenith, double reference_elevation, double target_elevation);

};

#endif // I3CORSIKAREADER_UTILS_H_INCLUDED
