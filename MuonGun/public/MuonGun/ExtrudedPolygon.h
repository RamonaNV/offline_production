/** $Id: ExtrudedPolygon.h 146654 2016-06-01 18:31:33Z cweaver $
 * @file
 * @author Jakob van Santen <jakob.van.santen@desy.de>
 *
 * $Revision: 146654 $
 * $Date: 2016-06-01 12:31:33 -0600 (Wed, 01 Jun 2016) $
 */

#ifndef I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED
#define I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED

#include <MuonGun/SamplingSurface.h>
#include <phys-services/surfaces/detail/ExtrudedPolygonBase.h>
#include <MuonGun/UprightSurface.h>

namespace I3MuonGun {

class ExtrudedPolygon : public detail::UprightSurface<I3Surfaces::ExtrudedPolygonBase<SamplingSurface > > {
private:
	typedef I3Surfaces::ExtrudedPolygonBase<SamplingSurface > ExtrudedPolygonBase;
	typedef detail::UprightSurface<ExtrudedPolygonBase > Base;
public:
	virtual ~ExtrudedPolygon();
	ExtrudedPolygon(const std::vector<I3Position> &points, double padding=0.) : Base(points, padding) {};
	
	virtual bool operator==(const SamplingSurface&) const
	{
		return false;
	}

protected:
	// UprightSurface interface
	double GetTopArea() const { return ExtrudedPolygonBase::GetCapArea(); };
	double GetSideArea() const { return ExtrudedPolygonBase::GetAverageSideArea(); };
	double GetLength() const { return ExtrudedPolygonBase::GetLength(); };
	std::pair<double, double> GetZRange() const { return ExtrudedPolygonBase::GetZRange(); };

private:
	ExtrudedPolygon() {}
	friend class icecube::serialization::access;
	template <typename Archive>
	void serialize(Archive &, unsigned);
};

I3_POINTER_TYPEDEFS(ExtrudedPolygon);

}

I3_CLASS_VERSION(I3MuonGun::ExtrudedPolygon, 0);

#endif // I3MUONGUN_EXTRUDEDPOLYGON_H_INCLUDED
