/**
 *
 * @file I3MCNKGInfo.h
 *
 * @brief Store the data points for the NKG samples provided by CORSIKA
 * @author Peter Niessen Fri Feb 23 10:09:26 EST 2007
 * @date $Date: 2017-10-19 17:26:27 -0600 (Thu, 19 Oct 2017) $
 * @version $Revision:$
 *
 * copyright (C) 2007
 * the IceCube collaboration
 *
 * $Id: I3MCNKGInfo.h 158928 2017-10-19 23:26:27Z cweaver $
 */

// multiple include protection
#ifndef I3MNKGINFO_H_INCLUDED
#define I3MNKGINFO_H_INCLUDED

#include <vector>
#include <map>

#include <dataclasses/Utility.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/physics/I3Particle.h>

struct I3MCNKGPoint {

  /**
   * where is this point?
   */
  I3Position Position;

  /**
   * see CORSIKA manual for details on this
   */
  double LD1; /// lateral distribution density observation level 1, in m^-2
  double LD2; /// lateral distribution density observation level 2, in m^-2

  /**
   * a constructor;
   * sets contents to NAN
   */
  I3MCNKGPoint ();

  /**
   * a virtual destructor
   * weiss nix? Macht nix!
   */
  virtual ~I3MCNKGPoint ();

  /**
   * Write a human readable representation
   */
  std::ostream& Print(std::ostream&) const;
  
  /**
   * the streamer necessary for boosting it to a file
   */
  template <class Archive> void serialize (Archive &ar, unsigned version);

};

/**
 * aggregate many points into a map
 */
typedef I3Vector<I3MCNKGPoint> I3MCNKGInfoList;

/**
 * A class to do the interpolation of values
 */
class I3MCNKGInterpolation {

 public:

  // constructor
  I3MCNKGInterpolation (const I3MCNKGInfoList &nkg_values,
			const I3Particle &primary);

  // destructor
  ~I3MCNKGInterpolation ();
  
  /**
   * Write a human readable representation
   */
  std::ostream& Print(std::ostream&) const;
  
  /**
   * Performs an interpolation
   */
  const double Interpolate (const I3Position &tank_position,
			    const I3Particle &primary) const;

  /**
   * Test the interpolation by sampling in between the grid
   */
  void TestInterpolate ();

 private:

  // need some typedefs
  typedef std::vector<int> f_table_index_t;
  typedef std::map<f_table_index_t, double> f_table_t;

  // semi - constant stuff
  double R_0;
  double R_MAX;
  double PHI_0;
  double PHI_MAX;

  // the grid spacing
  double H_LR;
  double H_PHI;

  // the function value of the nkg density as well as the derivatives
  // necessary for interpolation.
  f_table_t f_nkg;
  f_table_t df_dlr;
  f_table_t df_dphi;
  f_table_t d2f_dlr_dlr;
  f_table_t d2f_dlr_dphi;
  f_table_t d2f_dphi_dphi;

  // fill the grid to enable the interpolations
  void fillGrid (const I3MCNKGInfoList &nkg_values,
		 const I3Particle &primary);

  
  // this function will translate the coordinate r, phi into the right
  // index. 
  const f_table_index_t polar2index (const double r,
				     const double phi,
				     const double roundoff = 0.) const;

  // translate an index into carthesian coordinates. Basically used
  // for testing.
  void index2cartesian (const f_table_index_t &index,
			float &x, float &y) const;


  // the dimensions of the interpolating grid
  static const int N_R;
  static const int N_PHI;

  // the order to which to interpolate
  static const int ORDER;

  static int N_TEST;

  // looging
  SET_LOGGER ("I3MCNKGInterpolation");

};

std::ostream& operator<<(std::ostream&, const I3MCNKGPoint&);

I3_POINTER_TYPEDEFS (I3MCNKGPoint);
I3_POINTER_TYPEDEFS (I3MCNKGInfoList);

#endif
