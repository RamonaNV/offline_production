
#include "corsika-reader/I3CORSIKAReaderUtils.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/I3Direction.h"
#include <phys-services/surfaces/Sphere.h>

#include <boost/bimap.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace I3CORSIKAReaderUtils {

typedef std::map<int, int> CORSIKAtoPDGMap_t;

static const CORSIKAtoPDGMap_t CORSIKAtoPDGMap = {
	{1,	22},	// Gamma
	{2,	-11},	// e+
	{3,	11},	// e-
	{5,	-13},	// mu-
	{6,	13},	// mu+
	{7,	111},	// pi_0
	{8,	211},	// pi+
	{9,	-211},	// pi-
	{10,	130},	// K_L
	{11,	321},	// K+
	{12,	-321},	// K-
	{13,	2112},	// neutron
	{14,	2212},	// proton
	{15,	-2212},	// anti-proton
	{16,	310},	// K_S
	{17,	221},	// eta
	{18,	3122},	// lambda
	{19,	3222},	// sigma+
	{20,	3212},	// sigma 0
	{21,	3112},	// sigma-
	{22,	3322},	// xi 0
	{23,	3312},	// xi-
	{24,	3334},	// omega-
	{25,	-2112},	// anti-neutron
	{26,	-3122},	// lambda-bar
	{27,	-3112},	// sigma- bar
	{28,	-3212},	// sigma 0 bar
	{29,	-3222},	// sigma+ bar
	{30,	-3322},	// xi 0 bar
	{31,	-3312},	// xi+
	{32,	-3334},	// omega+
	{50,	223},	// small omega
	{51,	113},	// rho 0
	{52,	213},	// rho-plus
	{53,	-213},	// rho-minus
	{54,	2224},	// delta ++
	{55,	2214},	// delta +
	{56,	2114},	// delta 0
	{57,	1114},	// delta -
	{58,	-2224},	// delta ++ bar
	{59,	-2214},	// delta + bar
	{60,	-2114},	// delta 0 bar
	{61,	-1114},	// delta - bar
	{62,	313},	// k*0 {892},
	{63,	323},	// k*+ {892},
	{64,	-323},	// k*- {892},
	{65,	-313},	// k*0 bar {892},
	{66,	12},	// nu_e
	{67,	-12},	// anti-nu_e
	{68,	14},	// nu_mu
	{69,	-14},	// anti-nu_mu
	// 71-76 are special CORSIKA interaction tracking
     // Had to add these extra entries to map to accommodate EHISTORY showers.
	{71,	221},	// eta
	{72,	221},	// eta
	{73,	221},	// eta
	{74,	221},	// eta
	{116,	421},	// D0
	{117,	411},	// D+
	{118,	-411},	// D-
	{119,	-421},	// D0 bar
	{120,	431},	// D_s +
	{121,	-431},	// D_s -
	{122,	441},	// eta_c
	{123,	423},	// D*0 {2007},
	{124,	413},	// D*+ {2010},
	{125,	-413},	// D*- {2010},
	{126,	-423},	// D*0 bar {2007},
	{127,	433},	// D*s+
	{128,	-433},	// D*s-
	{130,	443},	// J/psi
	{131,	-15},	// tau+
	{132,	15},	// tau-
	{133,	16},	// nu_tau
	{134,	-16},	// anti-nu_tau
	{137,	4122},	// lambda_c+
	{138,	4232},	// xi_c+
	{139,	4132},	// xi_c0
	{140,	4222},	// sigma_c++
	{141,	4212},	// sigma_c+
	{142,	4112},	// sigma_c0
	{143,	4322},	// xi_c prime +
	{144,	4312},	// xi_c prime 0
	{145,	4332},	// omega_c 0
	{149,	-4122},	// lambda_c-
	{150,	-4232},	// xi_c-
	{151,	-4132},	// xi_c0 bar
	{152,	-4222},	// sigma_c--
	{153,	-4212},	// sigma_c-
	{154,	-4112},	// sigma_c0 bar
	{155,	-4322},	// xi_c prime -
	{156,	-4312},	// xi_c prime 0 bar
	{157,	-4332},	// omega_c 0 bar
	{161,	4224},	// sigma_c * ++
	{162,	4214},	// sigma_c * +
	{163,	4114},	// sigma_c * 0
	{171,	-4224},	// sigma_c * --
	{172,	-4214},	// sigma_c * -
	{173,	-4114},	// sigma_c * 0 bar
	{176,	511},	//B 0
	{179,	-511},	//B 0 bar
	{177,	521},	//B +
	{178,	-521},	//B -
	{180,	531},	//B_s 0
	{181,	-531},	//B_s 0 bar
	{182,	541},	//B_c +
	{183,	-541},	//B_c -
	{184,	5122},	//lambda_b 0
	{185,	5112},	//sigma_b -
	{186,	5222},	//sigma_b +
	{187,	5232},	//xi_b 0
	{188,	5132},	//xi_b -
	{189,	5332},	//omega_b -
	{190,	-5122},	//lambda_b 0 bar
	{191,	-5112},	//sigma_b + bar
	{192,	-5222},	//sigma_b - bar
	{193,	-5232},	//xi_b 0 bar
	{194,	-5132},	//xi_b +
	{195,	-5332}	//omega_b +
};

template <typename T1, typename T2>
std::pair<T2,T1> reverse(const std::pair<T1,T2> &p)
{
	return std::make_pair(p.second, p.first);
}

static CORSIKAtoPDGMap_t PDGtoCORSIKAMap(
    boost::make_transform_iterator(CORSIKAtoPDGMap.begin(), reverse<int,int>),
    boost::make_transform_iterator(CORSIKAtoPDGMap.end(), reverse<int,int>)
);

int32_t
CorsikaToPDG(int corsika_id)
{
	// Some event generators (e.g. DPMJET) emit history entries
	// with PDG codes that CORSIKA's numbering scheme can't
	// accomodate (e.g. phi mesons).
	if (corsika_id == 9999)
		return I3Particle::unknown;
	if (corsika_id >= 200)
		return 1000000000 + /* PDG Nucleus code */
		    10000 * (corsika_id % 100) + /* A */
		    10 * (corsika_id / 100); /* Z */
	
	auto i = CORSIKAtoPDGMap.find(corsika_id);
	if (i == CORSIKAtoPDGMap.cend())
		log_fatal("Unknown CORSIKA ID %d", corsika_id);

	return i->second;
}

int32_t
PDGToCorsika(int pdg_id)
{
	if (pdg_id >= 1000000000)
		return (pdg_id % 10000)*10 /* A */
		    + (pdg_id/10000)%1000; /* Z */
	auto i = PDGtoCORSIKAMap.find(pdg_id);
	if (i == PDGtoCORSIKAMap.cend())
		log_fatal("No CORSIKA equivalent for PDG ID %d", pdg_id);

	return i->second;
}

// Distance to the surface
double
GetSlantDepth(const I3Direction &dir, const I3Position &pos, double altitude)
{
	return -I3Surfaces::Sphere(altitude-I3Constants::OriginElev,
	    EarthRadius+altitude)
	    .GetIntersection(pos, dir).first;
}

double
LocalZenith(double reference_zenith, double reference_elevation, double target_elevation)
{
	return std::asin(std::sin(reference_zenith)*(EarthRadius+reference_elevation)/(EarthRadius+target_elevation));
}

}
