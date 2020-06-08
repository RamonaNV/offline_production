#include <cerrno>
#include <stdexcept>
#include <MuonGun/SplineTable.h>
#include <icetray/I3Logging.h>
#include <serialization/binary_object.hpp>

extern "C" {
	#include <photospline/bspline.h>
}

namespace I3MuonGun {

SplineTable::SplineTable() : bias_(0)
{
  memset(&table_, 0, sizeof(struct splinetable));
}

SplineTable::SplineTable(const std::string &path)
{
	if (readsplinefitstable(path.c_str(), &table_) != 0)
		throw std::runtime_error("Couldn't read spline table " + path);
	if (splinetable_read_key(&table_, SPLINETABLE_DOUBLE, "BIAS", &bias_))
		bias_ = 0;
}

SplineTable::~SplineTable()
{
	splinetable_free(&table_);
}

bool
SplineTable::operator==(const SplineTable &other) const
{
	if (bias_ != other.bias_)
		return false;
	// Same dimensions
	if (table_.ndim != other.table_.ndim)
		return false;
	// Same spline order
	if (!std::equal(table_.order, table_.order + table_.ndim, other.table_.order))
		return false;
	// Same knot grid
	if (!std::equal(table_.nknots, table_.nknots + table_.ndim, other.table_.nknots))
		return false;
	// Same region of support
	for (int i=0; i < table_.ndim; i++)
		if (table_.extents[i][0] != other.table_.extents[i][0] || table_.extents[i][1] != other.table_.extents[i][1])
			return false;
	// Same size of coefficient grid
	if (!std::equal(table_.naxes, table_.naxes + table_.ndim, other.table_.naxes))
		return false;
	size_t size = 1;
	for (int i=0; i < table_.ndim; i++)
		size *= size_t(table_.naxes[i]);
	// Same coefficient grid
	if (!std::equal(table_.coefficients, table_.coefficients + size, other.table_.coefficients))
		return false;
	
	return true;
}

int
SplineTable::Eval(double *coordinates, double *result) const
{
	std::vector<int> centers(unsigned(table_.ndim));
	
	if (tablesearchcenters(&table_, coordinates, &centers[0]) == 0)
		*result = ndsplineeval(&table_, coordinates, &centers[0], 0);
	else
		return EINVAL;
	
	*result -= bias_;
	
	return 0;
}

std::pair<double, double>
SplineTable::GetExtents(int dim) const
{
	if (dim < 0 || dim >= table_.ndim)
		throw std::out_of_range("Dimension index out of range");
	return std::make_pair(table_.extents[dim][0], table_.extents[dim][1]);
}

template <typename Archive>
void
SplineTable::save(Archive &ar, unsigned version) const
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	splinetable_buffer buf;
	buf.mem_alloc = &malloc;
	buf.mem_realloc = &realloc;
	writesplinefitstable_mem(&buf, &table_);
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("NBytes", buf.size);
	ar & make_nvp("FITSFile", icecube::serialization::make_binary_object(buf.data, buf.size));
	free(buf.data);
}

template <typename Archive>
void
SplineTable::load(Archive &ar, unsigned version)
{
	if (version > 0)
		log_fatal_stream("Version "<<version<<" is from the future");
	
	splinetable_buffer buf;
	buf.mem_alloc = &malloc;
	buf.mem_realloc = &realloc;
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("NBytes", buf.size);
	buf.data = buf.mem_alloc(buf.size);
	ar & make_nvp("FITSFile", icecube::serialization::make_binary_object(buf.data, buf.size));
	readsplinefitstable_mem(&buf, &table_);
	free(buf.data);
}

}

I3_SPLIT_SERIALIZABLE(I3MuonGun::SplineTable);

