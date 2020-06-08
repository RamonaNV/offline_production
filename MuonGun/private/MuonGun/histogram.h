/** $Id: histogram.h 128654 2015-02-04 18:34:51Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 128654 $
 * $Date: 2015-02-04 11:34:51 -0700 (Wed, 04 Feb 2015) $
 */

#ifndef I3MUONGUN_HISTOGRAM_H_INCLUDED
#define I3MUONGUN_HISTOGRAM_H_INCLUDED

#include <boost/function.hpp>
#include <boost/make_shared.hpp>
#include <boost/variant.hpp>

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>
#include <boost/array.hpp>

// An n-dimensional histogram class, loosely inspired by dashi.histogram

namespace I3MuonGun {
	
namespace histogram {

namespace binning {

/**
 * @brief Base class for binning schemes
 */ 
class scheme {
public:
	/** Get the bin index of the given value */
	virtual long index(double value) const = 0;
	
	virtual ~scheme();
	
	/** Return the edges of the bins */
	const std::vector<double>& edges() const
	{ return edges_; }
protected:
	std::vector<double> edges_;
};

/**
 * @brief A non-equispaced binning scheme
 *
 * In the general case, the proper bin can be found in logarithmic time
 * by binary search
 */
class general : public scheme {
public:
	/**
	 * Construct a binning scheme from the given ordered list
	 * of bin edges, inserting under- and overflow bins as necessary
	 */
	general(const std::vector<double> &edges)
	{
		if (edges.front() > -std::numeric_limits<double>::infinity())
			edges_.push_back(-std::numeric_limits<double>::infinity());
		std::copy(edges.begin(), edges.end(), std::back_inserter(edges_));
		if (edges.back() < std::numeric_limits<double>::infinity())
			edges_.push_back(std::numeric_limits<double>::infinity());
	}
	virtual ~general();
	
	static boost::shared_ptr<general> create(const std::vector<double> &edges)
	{
		return boost::make_shared<general>(edges);
	}
	
	long index(double value) const
	{
		long j = /*static_cast<long>*/(std::distance(edges_.begin(),
		    std::upper_bound(edges_.begin(),
		    edges_.end(), value)));
		assert(j > 0);
		return j-1;
	}
};

/**
 * @brief Trivial linear mapping
 */
struct identity {
	static inline double map(double v) { return v; }
	static inline double imap(double v) { return v; }
};

/**
 * @brief Bin edges linear in @f$ \log_{10}{x} @f$
 */
struct log10 {
	static inline double map(double v) { return std::pow(10, v); }
	static inline double imap(double v) { return std::log10(v); }
};

/**
 * @brief Bin edges linear in @f$ \cos{\theta} @f$
 */
struct cosine {
	static inline double map(double v) { return std::acos(v); }
	static inline double imap(double v) { return std::cos(v); }
};

/**
 * @brief Bin edges linear in @f$ x^N @f$
 *
 * @tparam N the exponent of the power law
 */
template <int N>
struct power {
	static inline double map(double v) { return std::pow(v, N); }
	static inline double imap(double v) { return std::pow(v, 1./N); }
};

/** @cond */
template <>
struct power<2> {
	static inline double map(double v) { return v*v; }
	static inline double imap(double v) { return std::sqrt(v); }
};
/** @endcond */

/**
 * @brief An equispaced binning scheme
 *
 * In this optimal case the bin edges are uniform under some
 * transformation between set limits and the bin index can be
 * found in constant time
 *
 * @tparam Transformation the transformation that makes the bin edges
 *                        equispaced
 */
template <typename Transformation = identity >
class uniform : public scheme {
public:
	uniform(double low, double high, size_t nsteps)
	    : offset_(Transformation::imap(low)),
	    range_(Transformation::imap(high)-Transformation::imap(low)),
	    min_(map(0)), max_(map(1)), nsteps_(nsteps)
	{
		edges_.reserve(nsteps+2);
		edges_.push_back(-std::numeric_limits<double>::infinity());
		for (size_t i = 0; i < nsteps; i++)
			edges_.push_back(map(i/double(nsteps_-1)));
		edges_.push_back(std::numeric_limits<double>::infinity());
	}
	virtual ~uniform() {};
	
	static boost::shared_ptr<uniform<Transformation> > create(double low,
	    double high, size_t nsteps)
	{
		return boost::make_shared<uniform<Transformation> >(low, high, nsteps);
	}
	
	long index(double value) const
	{
		if (value < min_)
			return 0;
		else if (value >= max_)
			return static_cast<long>(edges_.size()-2);
		else {
			return static_cast<long>(floor((nsteps_-1)*imap(value)))+1;
		}
	}
private:
	inline double map(double value) const
	{
		return Transformation::map(range_*value + offset_);
	}
	inline double imap(double value) const
	{
		return (Transformation::imap(value)-offset_)/range_;
	}
	
	double offset_, range_, min_, max_;
	size_t nsteps_;
};

};

/**
 * @brief A base class for histograms of arbitrary type
 *        and dimensionality
 */
class histogram_base {
public:
	typedef boost::multi_array_types::index index;
	typedef boost::multi_array_types::size_type size_type;
	
	virtual ~histogram_base();
	
	/** Get a pointer to the memory backing the bin contents */
	virtual const char* raw_bincontent() const = 0;
	/** Get a pointer to the memory backing the sum of squared weights */
	virtual const char* raw_squaredweights() const = 0;
	
	/** Get the bin edges */
	virtual std::vector<std::vector<double> > binedges() const = 0;
	
	/** Get the dimensionality of the histogram */
	virtual size_type ndim() const = 0;
	/** Get the number of elements along each dimension */
	virtual std::vector<size_type> shape() const = 0;
	/** Get the stride required to move sequentially along each dimension */
	virtual std::vector<index> strides() const = 0;

};

/**
 * @brief An N-dimensional histogram
 *
 * @tparam N number of dimensions
 * @tparam T type weight stored in the histogram
 */
template <size_t N, typename T = double>
class histogram : public histogram_base {
public:
	typedef boost::array<boost::variant< std::vector<double>, boost::shared_ptr<binning::scheme> >, N> bin_specification;
public:
	/** Construct with non-uniform bins in all dimensions */
	histogram(const boost::array<std::vector<double>, N> &edges)
	{
		for (size_t i=0; i < N; i++)
			binners_[i] = boost::make_shared<binning::general>(edges[i]);
	
		make_datacube();
	}
	
	/** Construct with uniform bins in all dimensions */
	histogram(const boost::array<boost::shared_ptr<binning::scheme>, N> &schemes)
	    : binners_(schemes)
	{
		make_datacube();
	}
	
	/** Construct with a mix of uniform and non-uniform bins in different dimensions */
	histogram(const boost::array<boost::variant< std::vector<double>, boost::shared_ptr<binning::scheme> >, N> &schemes)
	{
		for (size_t i=0; i < N; i++)
			binners_[i] = boost::apply_visitor(bin_visitor(), schemes[i]);
		make_datacube();
	}
	
	/**
	 * Add a sample to the histogram
	 *
	 * @param[in] values value of the sample in each dimension
	 * @param[in] weight weight to assign to the sample
	 */
	void fill(const boost::array<double, N> &values, T weight=1)
	{
		boost::array<typename datacube_type::index, N> idx;
		if (!std::isfinite(weight))
			return;
		for (size_t i=0; i < N; i++) {
			if (std::isnan(values[i]))
				return;
			idx[i] = binners_[i]->index(values[i]);
		}
		
		bincontent_(idx) += weight;
		squaredweights_(idx) += weight*weight;
	}
	
	const boost::array<std::vector<double>, N> & edges() const
	{
		return edges_;
	}

	const T* bincontent() const
	{
		return bincontent_.data();
	}
	
	const T* squaredweights() const
	{
		return squaredweights_.data();
	}
	
public:
	// template-agnostic part of the interface
	const char* raw_bincontent() const
	{
		return reinterpret_cast<const char*>(bincontent_.data());
	}
	
	const char* raw_squaredweights() const
	{
		return reinterpret_cast<const char*>(squaredweights_.data());
	}
	
	std::vector<std::vector<double> > binedges() const
	{
		std::vector<std::vector<double> > edges;
		edges.assign(edges_.begin(), edges_.end());
		return edges;
	}
	
	size_type ndim() const
	{
		return N;
	}
	
	std::vector<size_type> shape() const
	{
		std::vector<size_type> shape;
		shape.assign(bincontent_.shape(), bincontent_.shape()+N);
		
		return shape;
	}
	
	std::vector<index> strides() const
	{
		std::vector<index> strides;
		strides.assign(bincontent_.strides(), bincontent_.strides()+N);
		
		return strides;
	}
	
private:
	void make_datacube()
	{
		boost::array<size_t, N> dims;
		for (size_t i=0; i < N; i++) {
			edges_[i] = binners_[i]->edges();
			dims[i] = edges_[i].size()-1u;			
		}
		
		bincontent_.resize(dims);
		squaredweights_.resize(dims);
		
		std::fill(bincontent_.data(), bincontent_.data()+bincontent_.size(), 0);
		std::fill(squaredweights_.data(), squaredweights_.data()+squaredweights_.size(), 0);
	}
	
	/** @cond */
	typedef boost::shared_ptr<binning::scheme > scheme_ptr;
	struct bin_visitor : public boost::static_visitor<scheme_ptr> {
		scheme_ptr operator()(const scheme_ptr & v) const
		{ return v; }
		
		scheme_ptr operator()(const std::vector<double> & v) const
		{ return binning::general::create(v); }
	};
	/** @endcond */

private:
	boost::array<scheme_ptr, N> binners_;
	boost::array<std::vector<double>, N> edges_;
	
	typedef boost::multi_array<T, N> datacube_type;
	datacube_type bincontent_;
	datacube_type squaredweights_;

};

} // namespace histogram

} // namespace I3MuonGun

#endif
