/**
 * copyright  (C) 2004
 * the IceCube Collaboration
 * $Id: I3GaussianPMTPulse.h 20757 2006-06-13 03:56:16Z olivas $
 *
 * @file I3GaussianPMTPulse.h
 * @version $Revision: 1.5 $
 * @date $Date: 2006-06-12 23:56:16 -0400 (Mon, 12 Jun 2006) $
 * @author deyoung
 *
 */
#ifndef I3GAUSSIANPMTPULSE_H_INCLUDED
#define I3GAUSSIANPMTPULSE_H_INCLUDED

#include "icetray/I3FrameObject.h"
#include "icetray/I3Tray.h"
#include "icetray/OMKey.h"
#include "dataclasses/I3Map.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Constants.h"
#include <vector>

using namespace I3Constants;
using namespace I3Units;

/**
 * @brief Implementation of PMTPulse for Gaussian pulse shape
 * 
 * This class records the true (simulated) voltage, as a function of
 * time, produced by a single photoelectron (hit), using a Gaussian
 * model. 
 */
class I3GaussianPMTPulse : public I3FrameObject {
    
public:
  /**
   * constructor
   */
  I3GaussianPMTPulse() : normalization_(0.), 
			 sigma_(0.),
			 timeZero_(0.),
			 pedestal_(0.) {};
  
  /** 
   * Returns the voltage at the given time.  Voltage is negative.
   */  
   double GetVoltage(const double time) {
    return (- normalization_ / (sqrt(2 * pi) * sigma_)
	    * exp( - (time - timeZero_) * (time - timeZero_) 
		   / (2 * sigma_ * sigma_)) + pedestal_);
  };

  /** 
   * Returns the peak voltage for this pulse.  The pulse is negative,
   * so the peak value is a negative number.
   */ 
   double GetPeakVoltage() {
    return ( -normalization_ / (sqrt(2 * pi) * sigma_) + pedestal_);
  };

  /** 
   * Returns the time at which the pulse reaches its peak voltage.
   */ 
   double GetPeakTime() {
    return timeZero_;
  };

  /** 
   * Set the zero time for the Gaussian describing the pulse.  This is
   * the time at which the pulse reaches its peak voltage.
   */  
   void SetPeakTime(const double time) {
    timeZero_ = time;
  };

  /**
   * Get the time at which the pulse first crosses the given voltage
   * threshold.  Remember that the pulse is negative.
   */
   double GetStartTime(const double threshold) {
    float thres = (threshold >= GetPedestal() ? 
		     GetPedestal() - 1. * microvolt : 
		     threshold);
    return timeZero_ - sigma_ * 
      sqrt(-2 * log((pedestal_-thres) / normalization_ * sigma_ * sqrt(2*pi)));
  };

  /** 
   * Return the pedestal on which the Gaussian peak lies.  This is a
   * linear offset to the pulse voltage curve.
   */
   double GetPedestal() {
    return pedestal_;
  };

  /**
   * Set the pedestal for this channel.
   */
   void SetPedestal(const double ped) {
    pedestal_ = ped;
  };

  /**
   * Return the Gaussian width of the pulse.
   */
   double GetSigma() {
    return sigma_;
  };

  /**
   * Set the Gaussian width of the pulse.
   */
   void SetSigma(const double sigma) {
    sigma_ = sigma;
  };

  /**
   * Get the normalization parameter of the Gaussian describing this
   * pulse.  This is *not* the peak voltage.
   */
   double GetNormalization() {
    return normalization_;
  };

  /** 
   * Set the normalization parameter for this pulse.  This is *not*
   * the peak voltage.
   */
   void SetNormalization(const double norm) {
    normalization_ = norm;
  };

private:

  double normalization_;
  double sigma_;
  double timeZero_;
  double pedestal_;

  friend class icecube::serialization::access;

  template <class Archive> void serialize(Archive & ar, unsigned version);
};

I3_POINTER_TYPEDEFS(I3GaussianPMTPulse);

typedef std::vector<I3GaussianPMTPulse> I3GaussianPMTPulseList;
I3_POINTER_TYPEDEFS(I3GaussianPMTPulseList);

typedef I3Map<OMKey, std::vector<I3GaussianPMTPulse> > I3GaussianPMTPulseListMap;
I3_POINTER_TYPEDEFS(I3GaussianPMTPulseListMap);

#endif

