#include <iostream>
#include <cmath>
#include <string>
#include <sim-services/I3SumGenerator.h>

using namespace std;
/**
 * Constructor with arguments, just calls Initialise method
 */
I3SumGenerator::I3SumGenerator(I3RandomServicePtr r,boost::function<double(double)> fun,
			       double xlo, double xhi, int nbins, 
			       int switchgauss, double PLow, int nBinsLow, 
			       double PHigh, int nBinsHigh)
{
  Initialise(r,fun,xlo,xhi,nbins,switchgauss,PLow,nBinsLow,PHigh,nBinsHigh);
}
/**
 * Initialise I3SumGenerator
 */
void I3SumGenerator::Initialise(I3RandomServicePtr r, boost::function<double(double)> fun,
				double xlo, double xhi, int nbins, 
                                int switchgauss, double PLow, int nBinsLow, 
                                double PHigh, int nBinsHigh)
{
  int bin,bin1;
  vector<double>total;
  vector< vector<double> > P;

  /**
   * Initialise class variables
   */
  random_ = r;
  xLo_ = xlo;
  xHi_ = xhi;
  switchGauss_ = switchgauss;
  nBins_ = nbins;
  if(nBins_<2) nBins_=2;
  PLow_ = PLow;
  nBinsLow_ = nBinsLow;
  PHigh_ = PHigh;
  nBinsHigh_ = nBinsHigh;

  double dx( (xhi - xlo) / nBins_ );

  /** 
   * Sum (integral) for different number of terms
   */
  total.resize(switchgauss);
  /**
   * Vectors holding probabilities for different values of sum 
   */
  P.resize(switchgauss);
  /**
   * Unnormalised integral (sum over sampling points, really) for single term
   */
  double x( xlo+0.5*dx );
  double sumvalue(0),fvalue;
  for(bin=1; bin <= nBins_ ; bin++) {
    fvalue = fun(x);
    sumvalue += fvalue;
    x+=dx;
  }

  /**
   * Normalise probabilities, and calculate mean and sigma for single term
   */
  expectVal_=0;
  double sumsq(0);
  x = xlo+0.5*dx;
  P[1].resize(nBins_+1);
  P[1][0]=0;
  for(bin=1; bin <= nBins_ ; bin++) {
    fvalue = fun(x) / sumvalue;
    expectVal_ += x*fvalue;
    sumsq += x*x*fvalue;
    P[1][bin] = fvalue;
    x+=dx;
  }
  stdDev_ = sqrt(sumsq - expectVal_*expectVal_);
  // We normalised to unity:
  total[1] = 1;

  /**
   * Loop over number of terms from 2 up, recursively filling the vectors of 
   * sampled probability densities. The number of points is the same for different
   * number of terms, so the spacing is equal to the number of terms times dx.
   */
  for(int terms=2; terms<switchgauss; terms++) 
    {
      P[terms].resize(nBins_+1);
      total[terms] = 0;
      P[terms][0]=0;
      int bin,binold;
      // Loop over bins (samples) for current number of terms
      for(bin=0; bin <= nBins_ ; bin++) {
        P[terms][bin] = 0;
	// For each bin, loop over bins (samples) in terms-1 vector, stopping when it
	// becomes too big for the wanted value of bin
	for(binold=0; binold <= nBins_ && (terms-1)*binold <= terms*bin ;binold++) {
	  // Work out the bin (sample) number in the single term vector corresponding
	  // to the added term, and check that a solution exists
	  bin1 = terms*bin - (terms-1)*binold;
          if(bin1 <= nBins_) 
	    {
	      // Add the probabilities contributing to bin
	      P[terms][bin] += P[terms-1][binold]*P[1][bin1];
	    }
	}
	// Sum the probabilities for terms
	total[terms] += P[terms][bin];
      }
    }

  /**
   * Now calculate the normalised cumulants and invert to get a quick mapping from
   * cumulative probability to value of sum (X_)
   */
  X_.resize(switchgauss);  
  XLow_.resize(switchgauss);  
  XHigh_.resize(switchgauss);
  binStepLow_ = pow(PLow_,1./3.)/(double)nBinsLow_;
  binStepHigh_ = pow((1-PHigh_),1./3.)/(double)nBinsHigh_;
  for(int terms=1; terms<switchgauss; terms++) {
    // Convert to cumulant probabilities
    for(bin=1; bin <= nBins_ ; bin++) P[terms][bin] = P[terms][bin]/total[terms] + P[terms][bin-1] ;
    // Make sure to avoid round-off errors in probability sum
    P[terms][nBins_] = 1;
    //Fill the X_-vector for this number of terms
    X_[terms].resize(nBins_+2);
    //Start with the first x bin (x is the sum of terms)
    int binxBefore(0),binxAfter(1);
    double prob,dprob(1/(double)nBins_);
    /**
    * Loop over probability bins (sample points) and move the x bin along until 
    * it matches current probability
    */
    for(int binP=0;binP<nBins_;binP++) {
      prob=dprob*binP;
      while (P[terms][binxAfter] < prob) binxAfter++;
      binxBefore=binxAfter-1;
      // Find x value for current probability value by interpolation
      X_[terms][binP] = terms*(xlo + dx*( binxBefore + (prob-P[terms][binxBefore])/
				   (P[terms][binxAfter] - P[terms][binxBefore]) ));
    }
    X_[terms][nBins_] = terms*xhi;
    /**
     * Add an extra element just in case uniform random number sometimes equals one
     * in Generate()
     */
    X_[terms][nBins_+1] = X_[terms][nBins_];
    /**
     * Fill XLow_ array with values corresponding to low probabilities (lower tail of 
     * distribution to be generated
     */
    XLow_[terms].resize(nBinsLow_+2);
    binxAfter=1;
    for(int binP=0;binP<=nBinsLow_;binP++) {
      // Pick probabilities with non-uniform (cubically increasing) spacing
      prob=binStepLow_*binP;
      prob=prob*prob*prob;
      while (P[terms][binxAfter] <= prob) binxAfter++;
      binxBefore=binxAfter-1;
      XLow_[terms][binP] = terms*(xlo + dx*(binxBefore + (prob-P[terms][binxBefore])/
				     (P[terms][binxAfter] - P[terms][binxBefore])));
    }
    // Add extra element in case of rounding problems at generation:
    XLow_[terms][nBinsLow_+1] = XLow_[terms][nBinsLow_];
    /**
     * Upper tail. Sample cumulative probabilities at a non-uniform (cubically increasing) 
     * distance from probability=1.
     */
    XHigh_[terms].resize(nBinsHigh_+2);
    binxAfter = nBins_;
    binxBefore = binxAfter-1;
    for(int binP=0;binP<=nBinsHigh_;binP++) {
      prob=binStepHigh_*binP;
      prob=1-prob*prob*prob;
      while (P[terms][binxBefore] > prob) binxBefore--;
      binxAfter=binxBefore+1;
      XHigh_[terms][nBinsHigh_-binP] = terms*(xlo + dx*( binxBefore + (prob-P[terms][binxBefore])/
						  (P[terms][binxAfter] - P[terms][binxBefore])));
    }
    // For rounding errors:
    XHigh_[terms][nBinsHigh_+1] = XHigh_[terms][nBinsHigh_];
  }
}
 
double I3SumGenerator::Generate(int terms)
{
  double retval;
  if(terms < switchGauss_){
    /**
     * Low number of terms, read off sum from lookup table for a random
     * value of cumulative probability
     */
    double xi = random_->Uniform(1.);
    /**
     * If probability in lower or upper tail, use finer, cubically spaced tables
     * and interpolate linearly
     */
    if(xi < PLow_) {
      int bin = static_cast<int> (pow(xi,1./3.)/binStepLow_);
      double binfrac = (xi - pow(bin*binStepLow_,3))/
        (pow((bin+1)*binStepLow_,3) - pow((bin)*binStepLow_,3) );
      retval = XLow_[terms][bin] + binfrac*(XLow_[terms][bin+1] - XLow_[terms][bin]); 
    } else if (xi > PHigh_) {
      int bin = nBinsHigh_ - 1 - static_cast<int> (pow(1-xi,1./3.)/binStepHigh_);
      double Pbin = 1 - pow((nBinsHigh_-bin)*binStepHigh_,3);
      double PbinNext = 1 - pow((nBinsHigh_-(bin+1))*binStepHigh_,3);
      double binfrac = (xi - Pbin)/(PbinNext-Pbin);
      retval = XHigh_[terms][bin] + binfrac*(XHigh_[terms][bin+1] - XHigh_[terms][bin]); 
    } else {
      // Use full-range, linearly spaced, table
      xi*=nBins_;
      int bin = static_cast<int>(xi);
      double binfrac = xi - bin;
      retval = X_[terms][bin] + binfrac*(X_[terms][bin+1] - X_[terms][bin]);
    }
  } else {
    /**
     * Large number of terms, rely on central limit theorem
     */
    double sigma = stdDev_*sqrt((double)terms);
    retval = random_->Gaus(expectVal_*terms ,sigma);
    while(retval<terms*xLo_ || retval>terms*xHi_)
    retval = random_->Gaus(expectVal_*terms ,sigma);
  }
  return retval;
}
