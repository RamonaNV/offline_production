
#ifndef MUONGUN_SPLINETABLE_H_INCLUDED
#define MUONGUN_SPLINETABLE_H_INCLUDED

#include <string>

extern "C" {
	#include <photospline/splinetable.h>
}

#include "icetray/I3FrameObject.h"
#include "icetray/serialization.h"

namespace I3MuonGun {

/**
 * @brief An encapsulated, serializable interface to splinetable
 */
class SplineTable : public I3FrameObject {
public:
	/**
	 * @brief Read a spline table from a FITS file on disk
	 *
	 * @param[in] path The filesystem path to the FITS file
	 * @throws std::runtime_error if the file does not 
	 *         exist or is corrupt.
	 */
	SplineTable(const std::string &path);
	virtual ~SplineTable();
	
	/** @brief Default constructor, to be used only in serialization */
	SplineTable();

	/**
	 * @brief Evaulate the spline surface at the given coordinates
	 *
	 * @param[in]  x      Coordinates at which to evaluate
	 * @param[out] result Where to store the result
	 * @returns 0 on success.
	 */
	int Eval(double *x, double *result) const;

	/** @brief Return the number of dimensions of the spline surface */
	unsigned GetNDim() const { return unsigned(table_.ndim); };
	
	/**
	 * @brief Return the smallest and largest coordinates along the
	 * given axis where the spline surface has full support.
	 * 
	 * @param[in] dim the dimension to query
	 * @throws std::out_of_range if dim is < 0 or
	 *         >= the value returned by GetNDim()
	 */
	std::pair<double, double> GetExtents(int dim) const;
	
	/** @brief Deep comparison */
	bool operator==(const SplineTable &) const;
private:
	struct splinetable table_;
	double bias_;
	
	friend class icecube::serialization::access;
	template <typename Archive>
	void save(Archive &, unsigned) const;
	template <typename Archive>
	void load(Archive &, unsigned);
	I3_SERIALIZATION_SPLIT_MEMBER();
};

}

I3_CLASS_VERSION(I3MuonGun::SplineTable, 0);

#endif
