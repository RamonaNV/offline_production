
extern "C" {
	#include "photospline/splinetable.h"
	#include "photospline/bspline.h"
}

#include <I3Test.h>
#include <boost/filesystem.hpp>
#include <sys/time.h>
#include <limits>

namespace fs = boost::filesystem;

struct TableSet {
	fs::path abs, prob;
};

static void
splinetable_destructor(struct splinetable *table) {
	if (!table) return;
	
	splinetable_free(table);
	delete table;
	
}

static TableSet
get_splinetables()
{
	ENSURE(getenv("I3_TESTDATA") != NULL,
	    "I3_TESTDATA must be defined in the parent shell.");

	const std::string I3_TESTDATA(getenv("I3_TESTDATA"));
	
	fs::path abs_table(I3_TESTDATA + "/photospline/ems_z0_a0.pt.abs.fits");
	fs::path prob_table(I3_TESTDATA + "/photospline/ems_z0_a0.pt.prob.fits");
		
	ENSURE(fs::exists(abs_table), "Amplitude table does not exist.");
	ENSURE(fs::exists(prob_table), "Quantile table does not exists.");
	
	TableSet tabset;
	tabset.abs = abs_table;
	tabset.prob = prob_table;
	
	return tabset;
}

static boost::shared_ptr<struct splinetable>
load_splinetable(fs::path &fname)
{
	boost::shared_ptr<struct splinetable> table(new struct splinetable, splinetable_destructor);
	ENSURE(readsplinefitstable(fname.string().c_str(), table.get()) == 0, "Table can be read.");
	
	return table;
}

static void
compare_tables(struct splinetable *oldtable, struct splinetable *newtable)
{
	ENSURE_EQUAL(oldtable->ndim, newtable->ndim, "Number of dimensions match");
	size_t arraysize = 1;
	for (int i=0; i < oldtable->ndim; i++) {
		ENSURE_EQUAL(oldtable->naxes[i], newtable->naxes[i], "Coefficient array sizes match");
		arraysize *= oldtable->naxes[i];
		ENSURE_EQUAL(oldtable->order[i], newtable->order[i], "Spline orders match");
		ENSURE_EQUAL(oldtable->periods[i], newtable->periods[i], "Spline periods match");
		ENSURE_EQUAL(oldtable->nknots[i], newtable->nknots[i], "Size of knot vectors match");
		for (int j=0; j < oldtable->nknots[i]; j++)
			ENSURE_EQUAL(oldtable->knots[i][j], newtable->knots[i][j], "Knot positions match");
		ENSURE_EQUAL(oldtable->extents[i][0], newtable->extents[i][0], "Lower extents match");
		ENSURE_EQUAL(oldtable->extents[i][1], newtable->extents[i][1], "Upper extents match");
	}
	
	for (size_t i=0; i < arraysize; i++)
		ENSURE_EQUAL(oldtable->coefficients[i], newtable->coefficients[i], "Coefficients match");
	
	ENSURE_EQUAL(oldtable->naux, newtable->naux, "Number of auxiliary keywords match");
	for (int i=0; i < oldtable->naux; i++) {
		ENSURE_EQUAL(std::string(oldtable->aux[i][0]), std::string(newtable->aux[i][0]), "Keys match");
		ENSURE_EQUAL(std::string(oldtable->aux[i][1]), std::string(newtable->aux[i][1]), "Values match");
	}
}


TEST_GROUP(BoundaryIssues);

TEST(DiskFile)
{
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> oldtable = load_splinetable(tables.abs);
	
	fs::path tmp("photospline-write-test.fits");
	if (fs::exists(tmp))
		fs::remove(tmp);
	ENSURE_EQUAL(writesplinefitstable(tmp.string().c_str(), oldtable.get()), 0, "Table can be written");
	
	boost::shared_ptr<struct splinetable> newtable = load_splinetable(tmp);
	
	compare_tables(oldtable.get(), newtable.get());
	
	fs::remove(tmp);
}

TEST(MemoryFile)
{
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> oldtable = load_splinetable(tables.abs);
	
	struct splinetable_buffer buf;
	buf.mem_alloc = &malloc;
	buf.mem_realloc = &realloc;
	ENSURE_EQUAL(writesplinefitstable_mem(&buf, oldtable.get()), 0, "Table can be written");
	
	boost::shared_ptr<struct splinetable> newtable(new struct splinetable, splinetable_destructor);
	ENSURE_EQUAL(readsplinefitstable_mem(&buf, newtable.get()), 0, "Table can be read.");
	free(buf.data);
	
	compare_tables(oldtable.get(), newtable.get());
}

/*
 * Check that analytic convolution works and preserves the monotonicity
 * of the arrival-time CDF.
 */
TEST(Convolution)
{
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table = load_splinetable(tables.prob);
	const unsigned time_dim = 3;
	double sigma = 10.0;
	int n_knots = 3;
	double knots[3] = {-2*sigma, 0, 2*sigma};
	unsigned i;
	unsigned ndim = table->ndim;
	double tablecoords[ndim], base, nudge, q0_raw, q0_conv, q1_raw, q1_conv;
	const double eps = std::numeric_limits<double>::epsilon();
	const double tmax = table->extents[time_dim][1];
	int centers[ndim];
	
	for (i = 0; i < ndim; i++) {
		double low, high;
		low = table->extents[i][0];
		high = table->extents[i][1];
		tablecoords[i] = (low+high)/2.0;
	}
	
	tablecoords[time_dim] = table->knots[time_dim][0]+eps;
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q0_raw = ndsplineeval(table.get(), tablecoords, centers, 0);
	tablecoords[time_dim] = table->extents[time_dim][1];
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q1_raw = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE_DISTANCE(q1_raw-q0_raw, 1, 1e-3,
	    "Arrival-time CDF is normalized to within 0.1%%.");
	
	base = q0_raw;
	for (i = 1; i+table->order[time_dim]+1 < table->nknots[time_dim]; i++) {
		tablecoords[time_dim] = table->knots[time_dim][i];
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge - base >= 0,
		    "Arrival-time CDF is monotonic.");
		base = nudge;
	}
	
	int err = splinetable_convolve(table.get(), time_dim, knots, n_knots);
	ENSURE_EQUAL(err, 0, "Convolution succeeds.");
	
	tablecoords[time_dim] = table->knots[time_dim][0]*(1-eps);
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q0_conv = ndsplineeval(table.get(), tablecoords, centers, 0);

	base = q0_conv;
	for (i = 1; i < table->nknots[time_dim]; i++) {
		tablecoords[time_dim] = table->knots[time_dim][i];
		if (tablecoords[time_dim] > table->extents[time_dim][1])
			break;
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge - base >= 0,
		    "Arrival-time CDF remains monotonic.");
		base = nudge;
	}
	
	tablecoords[time_dim] = table->extents[time_dim][1];
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0);
	q1_conv = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(table->extents[time_dim][1] < tmax,
	    "t_max is slightly smaller after convolution.");
	ENSURE(q1_conv - q1_raw > -10*std::numeric_limits<float>::epsilon(),
	    "Maximum quantile (at the end of support) is not diminished.");
	
	ENSURE_DISTANCE(q1_raw-q0_raw, q1_conv-q0_conv, 10*std::numeric_limits<float>::epsilon(),
	    "Arrival-time CDF remains normalized after convolution.");
}

/*
 * The quantiles in the time dimension are continous. This can only be true
 * if the spline evaluation code is capable of handling points near the
 * edges of the knot fields that are supported by fewer than (order+1)
 * splines.
 */
TEST(QuantileContinuity)
{
	const unsigned time_dim = 3;
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table = load_splinetable(tables.prob);
	
	unsigned i;
	unsigned ndim = table->ndim;
	double tablecoords[ndim], base, nudge;
	const double eps = std::numeric_limits<double>::epsilon();
	int centers[ndim];
	
	for (i = 0; i < ndim; i++) {
		double low, high;
		low = table->extents[i][0];
		high = table->extents[i][1];
		tablecoords[i] = (low+high)/2.0;
	}
	
	/* Check the transition into the knot field from the left. */
	tablecoords[time_dim] = table->knots[time_dim][0];
	ENSURE(tablesearchcenters(table.get(), tablecoords, centers) != 0,
	    "tablesearchcenters() fails right at the left edge of the knot field.");
	
	tablecoords[time_dim] += eps;
	
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
	    "tablesearchcenters() succeeds just inside the left edge of the knot field.");
	ENSURE_EQUAL(centers[time_dim], table->order[time_dim],
	    "centers[time_dim] holds the first fully-supported knot index.");
	
	base = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(base >= 0, "Base quantile is positive.");
	ENSURE_DISTANCE(0, base, std::numeric_limits<float>::epsilon(),
	    "Time quantile is continuous "
	    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
		
	/*
	 * Now, step over the intertior knots, checking for continuity
	 * at every knot crossing.
	 */
	for (i = 1; i+1 < table->nknots[time_dim]; i++) {
		tablecoords[time_dim] = table->knots[time_dim][i]*(1-eps);
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
		    "tablesearchcenters() succeeds inside the knot field.");
		
		base = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(base >= 0);
		
		tablecoords[time_dim] = table->knots[time_dim][i];
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
		    "tablesearchcenters() succeeds inside the knot field.");
		
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge >= 0);
		
		ENSURE_DISTANCE(base, nudge, std::numeric_limits<float>::epsilon(),
		    "Time quantile is continuous "
		    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
		
		base = nudge;
		
		tablecoords[time_dim] = table->knots[time_dim][i]*(1+eps);
		ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
		    "tablesearchcenters() succeeds inside the knot field.");
		
		nudge = ndsplineeval(table.get(), tablecoords, centers, 0);
		ENSURE(nudge >= 0);
		
		ENSURE_DISTANCE(base, nudge, std::numeric_limits<float>::epsilon(),
		    "Time quantile is continuous "
		    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
	}

	/* Check the transition into the knot field from the right */
	tablecoords[time_dim] = table->knots[time_dim][table->nknots[time_dim]-1]*(1+eps);
	ENSURE(tablecoords[time_dim] > table->knots[time_dim][table->nknots[time_dim]-1]);
	ENSURE(tablesearchcenters(table.get(), tablecoords, centers) != 0,
	    "tablesearchcenters() fails right at the right edge of the knot field.");
	
	tablecoords[time_dim] = table->knots[time_dim][table->nknots[time_dim]-1];
	
	ENSURE_EQUAL(tablesearchcenters(table.get(), tablecoords, centers), 0,
	    "tablesearchcenters() succeeds just inside the right edge of the knot field.");
	ENSURE_EQUAL(centers[time_dim],
	    table->nknots[time_dim]-table->order[time_dim]-2,
	    "centers[time_dim] holds the first fully-supported knot index.");
	
	base = ndsplineeval(table.get(), tablecoords, centers, 0);
	
	ENSURE(base >= 0, "Base quantile is positive.");
	ENSURE_DISTANCE(0, base, std::numeric_limits<float>::epsilon(),
	    "Time quantile is continuous "
	    "(slope is less than FLOAT_EPSILON/DBL_EPSILON)");
}

TEST(ndssplineeval_vs_ndssplineeval_gradient)
{
	srand(42);
	
	TableSet tables = get_splinetables();
	boost::shared_ptr<struct splinetable> table = load_splinetable(tables.prob);
	
	const int ndim = table->ndim;
	ENSURE(ndim < 6);
	for (int i=0; i < 1000; i++) {
		/* Generate coordinates uniformly in the support of the spline */
		double x[ndim];
		int centers[ndim];
		for (int j=0; j < ndim; j++) {
			double p = double(rand())/double(RAND_MAX);
			x[j] = table->extents[j][0] + p*(table->extents[j][1]-table->extents[j][0]);
		}
		ENSURE_EQUAL(tablesearchcenters(table.get(), x, centers), 0);
		double evaluate;
		double evaluate_with_gradient[ndim+1];
		evaluate = ndsplineeval(table.get(), x, centers, 0);
		ndsplineeval_gradient(table.get(), x, centers, evaluate_with_gradient);
		ENSURE_EQUAL(evaluate, evaluate_with_gradient[0],
		    "ndsplineeval() and ndssplineeval_gradient() yield identical evaluates");
		for (int j=0; j < ndim; j++) {
			evaluate = ndsplineeval(table.get(), x, centers, 1 << j);
			ENSURE_EQUAL(evaluate, evaluate_with_gradient[1+j],
			    "ndsplineeval() and ndssplineeval_gradient() yield identical derivatives");
		}
	}
}

/*
 * bsplvb_simple() can be made to return sensical values anywhere
 * in the knot field.
 */

TEST(bsplvb_simple_vs_bspline)
{
	unsigned i;
	const unsigned n_knots = 10;
	const int order = 2;
	double x, *knots;
	int center, offset;
	float localbasis_bsplvb[order+1], localbasis_bspline[order+1];
	struct timeval thetime;
	std::vector<double> knotvec;
	
	/* Seed our crappy little PSRNG. */
	gettimeofday(&thetime, NULL);
	srand(thetime.tv_sec);
	/* Generate a random knot field. */
	for (i = 0; i < n_knots; i++)
		knotvec.push_back(double(rand())/double(RAND_MAX));
	std::sort(knotvec.begin(), knotvec.end());
	knots = &knotvec.front();
	
	/* Before the first fully-supported knot. */
	for (i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = order; /* First fully-supported spline. */
	
		ENSURE(int(i) <= center);
	
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order+1 /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
			std::numeric_limits<double>::epsilon());
		}
	}
	
	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order+1 /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
			std::numeric_limits<double>::epsilon());
		}
	}
	
	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order+1 /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
			std::numeric_limits<double>::epsilon());
		}
	}
}

/*
 * bspline_deriv_nonzero() can be made to return sensical values anywhere
 * in the knot field.
 */

TEST(bspline_deriv_nonzero_vs_bspline_deriv)
{
	unsigned i;
	const unsigned n_knots = 10;
	const int order = 2;
	double x, *knots;
	int center, offset;
	float localbasis_bsplvb[order+1], localbasis_bspline[order+1];
	struct timeval thetime;
	std::vector<double> knotvec;
	/* This calculation is less stable. */
	double tol = 10*std::numeric_limits<float>::epsilon();
	
	/* Seed our crappy little PSRNG. */
	gettimeofday(&thetime, NULL);
	srand(thetime.tv_sec);
	/* Generate a random knot field. */
	for (i = 0; i < n_knots; i++)
		knotvec.push_back(double(rand())/double(RAND_MAX));
	std::sort(knotvec.begin(), knotvec.end());
	knots = &knotvec.front();
	
	/* Before the first fully-supported knot. */
	for (i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = order; /* First fully-supported spline. */
	
		ENSURE(int(i) <= center);
	
		bspline_deriv_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline_deriv(knots, x, center + offset, order, 1);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
		}
	}
	
	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bspline_deriv_nonzero(knots, n_knots, x, center /* left */,
		    order, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline_deriv(knots, x, center + offset, order, 1);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
		}
	}

	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bspline_deriv_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb /* biatx */);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline_deriv(knots, x, center + offset, order, 1);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
		}
	}
}

/*
 * bspline_nonzero() can be made to return sensical values anywhere
 * in the knot field.
 */

TEST(bspline_nonzero_vs_bspline)
{
	unsigned i;
	const unsigned n_knots = 10;
	const int order = 2;
	double x, *knots;
	int center, offset;
	float localbasis_bsplvb[order+1], localbasis_bspline[order+1],
	    localbasis_bsplvb_deriv[order+1], localbasis_bspline_deriv[order+1];
	struct timeval thetime;
	std::vector<double> knotvec;
	/* This calculation is less stable. */
	double tol = 10*std::numeric_limits<float>::epsilon();
	
	/* Seed our crappy little PSRNG. */
	gettimeofday(&thetime, NULL);
	srand(thetime.tv_sec);
	/* Generate a random knot field. */
	for (i = 0; i < n_knots; i++)
		knotvec.push_back(double(rand())/double(RAND_MAX));
	std::sort(knotvec.begin(), knotvec.end());
	knots = &knotvec.front();

	/* Before the first fully-supported knot. */
	for (i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = order; /* First fully-supported spline. */
	
		ENSURE(int(i) <= center);
	
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb, localbasis_bsplvb_deriv);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			    bspline_deriv(knots, x, center + offset, order, 1);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
			
		}
	}
	
	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb, localbasis_bsplvb_deriv);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			    bspline_deriv(knots, x, center + offset, order, 1);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
			
		}
	}

	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bsplvb, localbasis_bsplvb_deriv);
	
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			    bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			    bspline_deriv(knots, x, center + offset, order, 1);
		}
	
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
			
		}
	}
}

/*
 * bspline_nonzero() gives the same result as bsplvb_simple(), ergo the
 * value bases used in ndsplineeval() and ndsplineeval_gradient() are
 * identical. Roundoff error in the derivative basis calculation
 * scales with the spline order, so we use a large number here.
 */

TEST(single_basis_vs_multi)
{
	unsigned i;
	const unsigned n_knots = 100;
	const int order = 5;
	double x, *knots;
	int center, offset;
	float localbasis_bspline_nonzero[order+1], localbasis_bspline_nonzero_deriv[order+1],
	    localbasis_bsplvb_simple[order+1], localbasis_bspline_deriv_nonzero[order+1];
	// bsplvb() may access up to *order* elements off either end
	// Pad accordingly.
	std::vector<double> knotvec(n_knots + 2*order, 0.);
	
	srand(0);
	/* Generate a random knot field. */
	for (i = 0; i < n_knots; i++)
		knotvec[i+order] = (double(rand())/double(RAND_MAX));
	std::sort(knotvec.begin()+order, knotvec.end()-order);
	knots = &knotvec[order];

	/* Before the first fully-supported knot. */
	for (i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = order; /* First fully-supported spline. */
	
		ENSURE(int(i) <= center);
	
		/* As used in ndssplineeval_gradient() */
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bspline_nonzero, localbasis_bspline_nonzero_deriv);
		/* As used in ndsplineeval() */
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order + 1 /* jhigh */, localbasis_bsplvb_simple);
		bspline_deriv_nonzero(knots, n_knots, x, center,
		    order, localbasis_bspline_deriv_nonzero);
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_EQUAL(localbasis_bspline_nonzero[offset], localbasis_bsplvb_simple[offset]);
			ENSURE_EQUAL(localbasis_bspline_nonzero_deriv[offset], localbasis_bspline_deriv_nonzero[offset]);
		
		}

	}

	/* Within the support. */
	for (i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bspline_nonzero, localbasis_bspline_nonzero_deriv);
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order + 1 /* jhigh */, localbasis_bsplvb_simple);
		bspline_deriv_nonzero(knots, n_knots, x, center,
		    order, localbasis_bspline_deriv_nonzero);		
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_EQUAL(localbasis_bspline_nonzero[offset], localbasis_bsplvb_simple[offset]);
			ENSURE_EQUAL(localbasis_bspline_nonzero_deriv[offset], localbasis_bspline_deriv_nonzero[offset]);
		}
	}

	/* After the last first fully-supported knot. */
	for (i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; /* Less than fully-supported */
		center = n_knots-order-2; /* Last fully-supported spline. */
	
		ENSURE(int(i) >= center);
	
		bspline_nonzero(knots, n_knots, x, center /* left */,
		    order /* jhigh */, localbasis_bspline_nonzero, localbasis_bspline_nonzero_deriv);
		bsplvb_simple(knots, n_knots, x, center /* left */,
		    order + 1 /* jhigh */, localbasis_bsplvb_simple);
		bspline_deriv_nonzero(knots, n_knots, x, center,
		    order, localbasis_bspline_deriv_nonzero);		
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_EQUAL(localbasis_bspline_nonzero[offset], localbasis_bsplvb_simple[offset]);
			ENSURE_EQUAL(localbasis_bspline_nonzero_deriv[offset], localbasis_bspline_deriv_nonzero[offset]);
		}
	}
}

