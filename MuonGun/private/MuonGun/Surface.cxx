/** $Id: Surface.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 12:24:47 -0600 (Mon, 31 Aug 2015) $
 */

#include <MuonGun/I3MuonGun.h>
#include <MuonGun/Surface.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <phys-services/I3RandomService.h>
#include <boost/bind.hpp>

namespace {

inline void
sort(std::pair<double, double> &pair)
{
	if (pair.first > pair.second) {
		double aux = pair.first;
		pair.first = pair.second;
		pair.second = aux;
	}
}

}

namespace simclasses {

Surface::~Surface() {}

template <typename Archive>
void
Surface::serialize(Archive &ar __attribute__ ((unused)), unsigned version __attribute__ ((unused)))
{}

SamplingSurface::~SamplingSurface() {}

double
SamplingSurface::SampleImpactRay(I3Position &impact, I3Direction &dir, I3RandomService &rng,
    double cosMin, double cosMax) const
{
	dir = SampleDirection(rng, cosMin, cosMax);
	impact = SampleImpactPosition(dir, rng);

	// Calculate projected area
	return GetArea(dir);
}

I3Direction
SamplingSurface::SampleDirection(I3RandomService &rng,
    double cosMin, double cosMax) const
{
	// Sample a direction proportional to the projected area 
	// of the surface.
	double maxarea = GetMaximumArea();
	I3Direction sampled_dir;
	do {
		sampled_dir = I3Direction(acos(rng.Uniform(cosMin, cosMax)),
		    rng.Uniform(0, 2*M_PI));
	} while (rng.Uniform(0, maxarea) > GetArea(sampled_dir));
	
	return sampled_dir;
}

template <typename Archive>
void
SamplingSurface::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Surface", base_object<Surface>(*this));
}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<Base>(*this));
}

Cylinder::~Cylinder() {}

AxialCylinder::AxialCylinder(double length, double radius, I3Position center)
    : length_(length/2.,length/2.), radius_(radius), center_(center)
{}

AxialCylinder::AxialCylinder(double lengthBefore, double lengthAfter, double radius, I3Position center)
    : length_(lengthBefore,lengthAfter), radius_(radius), center_(center)
{}

std::pair<double, double>
AxialCylinder::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	// Distance to point of closest approach to the center
	double to_center = (center_ - p)*dir;
	// Check distance of closest approach against cylinder radius
	if ((p + to_center*dir - center_).Magnitude() > radius_)
		return no_intersection();
	else
		return std::make_pair(to_center-length_.first, to_center+length_.second);
}

double
AxialCylinder::GetArea(const I3Direction &dir __attribute__((unused))) const
{
	return M_PI*radius_*radius_;
}

double
AxialCylinder::GetAcceptance(double cosMin, double cosMax) const
{
	return M_PI*radius_*radius_*(cosMax-cosMin);
}

double
AxialCylinder::GetMaximumArea() const
{
	return M_PI*radius_*radius_;
}

I3Direction
AxialCylinder::SampleDirection(I3RandomService &rng, double cosMin, double cosMax) const
{
	return I3Direction(std::acos(rng.Uniform(cosMin, cosMax)), rng.Uniform(0, 2*M_PI));
}

I3Position
AxialCylinder::SampleImpactPosition(const I3Direction &dir, I3RandomService &rng) const
{
	// Choose a position in a circle in axis-centered coordinates
	I3Position impact(std::sqrt(rng.Uniform(0, radius_*radius_)), 0, 0);
	impact.RotateZ(rng.Uniform(0, 2*M_PI));
	
	// Rotate into the transverse plane
	impact.RotateY(dir.GetZenith());
	impact.RotateZ(dir.GetAzimuth());
	// Shift from cylinder-centered to real coordinates
	impact += center_;
	// Shift back to the entry point
	impact -= length_.first*dir;
	
	return impact;
}

template <typename Archive>
void
AxialCylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("SamplingSurface", base_object<SamplingSurface>(*this));
	ar & make_nvp("Length", length_);
	ar & make_nvp("Radius", radius_);
	ar & make_nvp("Center", center_);
}

AxialCylinder::~AxialCylinder() {}

std::pair<double, double>
Sphere::GetIntersection(const I3Position &p, const I3Direction &dir) const
{
	std::pair<double, double> h(no_intersection());
	
	double x = p.GetX();
	double y = p.GetY();
	double z = p.GetZ() - originDepth_;
	
	double sinph = sin(dir.GetAzimuth());
	double cosph = cos(dir.GetAzimuth());
	double sinth = sin(dir.GetZenith());
	double costh = cos(dir.GetZenith());
	
	double b = (x*cosph + y*sinph)*sinth + (z + radius_)*costh;
	double d = b*b - (x*x + y*y + z*(z + 2*radius_));
	
	if (d > 0) {
		d = sqrt(d);
		h.first = b - d;
		h.second = b + d;
	}
	
	return h;
}

template <typename Archive>
void
Sphere::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	ar & make_nvp("Surface", base_object<Surface>(*this));
	ar & make_nvp("OriginDepth", originDepth_);
	ar & make_nvp("Radius", radius_);
}

Sphere::~Sphere() {}

}

namespace I3MuonGun {

template <typename Archive>
void
SamplingSurface::serialize(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");

	ar & make_nvp("Base", base_object<I3Surfaces::SamplingSurface>(*this));
}

SamplingSurface::~SamplingSurface() {}

std::pair<double, double>
Cylinder::GetZRange() const
{
	return std::make_pair(GetCenter().GetZ() - GetLength()/2., GetCenter().GetZ() + GetLength()/2.);
}

double
Cylinder::GetTopArea() const
{
	return M_PI*GetRadius()*GetRadius();
}

double
Cylinder::GetSideArea() const
{
	return 2*GetRadius()*GetLength();
}

bool
Cylinder::operator==(const SamplingSurface &s) const
{
	const Cylinder *other = dynamic_cast<const Cylinder*>(&s);
	if (!other)
		return false;
	else 
		return (GetRadius() == other->GetRadius() &&
		    GetLength() == other->GetLength() &&
		    GetCenter() == other->GetCenter());
}

template <typename Archive>
void
Cylinder::serialize(Archive &ar, unsigned version)
{
	if (version > 1)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	if (version == 0) {
		// Class hierarchy changed from v0 to v1, so we have to
		// deserialize by hand
		double radius, length;
		I3Position center;
		ar & make_nvp("SamplingSurface", base_object<I3Surfaces::SamplingSurface>(*this));
		ar & make_nvp("Length", length);
		ar & make_nvp("Radius", radius);
		ar & make_nvp("Center", center);
		SetLength(length);
		SetRadius(radius);
		SetCenter(center);
	} else {
		ar & make_nvp("Base", base_object<Base>(*this));
	}
}

}

// explicitly instantiate the base classes used
template class I3Surfaces::detail::CylinderBase<I3Surfaces::SamplingSurface>;
template class I3Surfaces::detail::CylinderBase<I3MuonGun::SamplingSurface>;
template class I3MuonGun::detail::UprightSurface<I3Surfaces::detail::CylinderBase<I3MuonGun::SamplingSurface> >;

I3_SERIALIZABLE(I3Surfaces::Surface);
I3_SERIALIZABLE(I3Surfaces::SamplingSurface);
I3_SERIALIZABLE(I3Surfaces::Cylinder);
I3_SERIALIZABLE(I3Surfaces::Sphere);
I3_SERIALIZABLE(I3Surfaces::AxialCylinder);

I3_SERIALIZABLE(I3MuonGun::SamplingSurface);
I3_SERIALIZABLE(I3MuonGun::Cylinder);


