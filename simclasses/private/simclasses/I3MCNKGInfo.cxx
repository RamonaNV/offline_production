/*******************************************************************/
/* I3MCNKGInfo.cxx                                                 */
/* Implementation of the NKG info                                  */
/* PN Fri Feb 23 10:21:09 EST 2007                                 */
/* $Id: I3MCNKGInfo.cxx 158928 2017-10-19 23:26:27Z cweaver $ */
/*******************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

#include <icetray/I3TrayHeaders.h>
#include <simclasses/I3MCNKGInfo.h>



/*******************************************************************/

// a constructor, set the values to NAN
I3MCNKGPoint::I3MCNKGPoint ()
  : Position (NAN, NAN, NAN),
    LD1 (NAN),
    LD2 (NAN) {
}

/*******************************************************************/

// a virtual destructor
I3MCNKGPoint::~I3MCNKGPoint () {
  // wer nix weiss und wer nix kann, f"angt zu simulieren an.
}

/*******************************************************************/

// the streamer necessary for boosting it to a file
template <class Archive> void I3MCNKGPoint::serialize (Archive &ar,
						      unsigned version) {

  ar &make_nvp ("Position", Position);
  ar &make_nvp ("LD1", LD1);
  ar &make_nvp ("LD2", LD2);

}

std::ostream& I3MCNKGPoint::Print(std::ostream& os) const{
  os << "I3MCNKGPoint: Position=" << Position << " LD1=" << LD1 << " LD2=" << LD2;
  return os;
}

std::ostream& operator<<(std::ostream& os, const I3MCNKGPoint& p){
  return(p.Print(os));
}

/*******************************************************************/


I3_SERIALIZABLE (I3MCNKGPoint);
I3_SERIALIZABLE (I3MCNKGInfoList);

/*******************************************************************/
/* here follows the interpolation implementation                   */
/*******************************************************************/

// some ugly hardcoding: In corsika, we have 10 r bins and 8 phi
// bins.
const int I3MCNKGInterpolation::N_R (10);
const int I3MCNKGInterpolation::N_PHI (8);

// the order of the interpolation
const int I3MCNKGInterpolation::ORDER (1);

int I3MCNKGInterpolation::N_TEST (0);

// constructor
I3MCNKGInterpolation::I3MCNKGInterpolation (const I3MCNKGInfoList &nkg_values,
					    const I3Particle &primary) {

  fillGrid (nkg_values, primary);

}

/*******************************************************************/

// destructor
I3MCNKGInterpolation::~I3MCNKGInterpolation () {
}
  
/*******************************************************************/

const double
I3MCNKGInterpolation::Interpolate (const I3Position &tank_position,
				   const I3Particle &primary) const {

  
  double x; // cartesian coordinates
  double y;
  double r; // polar coordinates
  double phi;

  
  x = tank_position.GetX () - primary.GetPos ().GetX ();
  y = tank_position.GetY () - primary.GetPos ().GetY ();

  r = sqrt (x * x + y * y);
  phi = atan2 (y, x);

  // make sure that r is contained
  if (r < R_0)
    r = R_0;

  double lr = log (r / R_0);

  // need to calculate the h_lr and h_phi
  // get the lowest value below the desired coordinate
  f_table_index_t grid_pos (polar2index (r, phi));
 
  // grid_pos will be maximally the last bin. so from above these
  // values, it will be extrapolation
  double h_lr = lr - H_LR * grid_pos[0];
  double h_phi = phi - H_PHI * grid_pos[1];

  log_debug ("i, r, h_lr/H_LR=%d, %f, %f",
	     grid_pos[0], r, h_lr / H_LR);
  //fprintf (stderr, "i, r, h_lr/H_LR=%d, %f, %f\n",
//	   grid_pos[0], r, h_lr / H_LR);

  // look up the derivatives. (Do it this way so that the method can
  // be const.)
  f_table_t::const_iterator f = f_nkg.find (grid_pos);
  f_table_t::const_iterator dlr = df_dlr.find (grid_pos);
  f_table_t::const_iterator dphi = df_dphi.find (grid_pos);
  f_table_t::const_iterator d2lrlr = d2f_dlr_dlr.find (grid_pos);
  f_table_t::const_iterator d2phiphi = d2f_dphi_dphi.find (grid_pos);
  f_table_t::const_iterator d2lrphi = d2f_dlr_dphi.find (grid_pos);

  // check if the necessary values are available
  if (f_nkg.end () == f)
    return -1.;
  if (df_dlr.end () == dlr)
    return -1.;
  if (df_dphi.end () == dphi)
    return -1.;
  if (d2f_dlr_dlr.end () == d2lrlr)
    return -1.;
  if (d2f_dphi_dphi.end () == d2phiphi)
    return -1.;
  if (d2f_dlr_dphi.end () == d2lrphi)
    return -1.;
  

  // f_nkg etc. store logarithm, thus, exponential is required
  return exp (f->second
	      + (ORDER > 0 ? h_lr * dlr->second : 0.)
	      + (ORDER > 0 ? h_phi * dphi->second : 0.)
	      + (ORDER > 1 ? 0.5 * (h_lr * h_lr * d2lrlr->second
				    + h_phi * h_phi * d2phiphi->second
				    + 2 * h_lr * h_phi * d2lrphi->second)
		 : 0.)
	      );

}

/*******************************************************************/

void I3MCNKGInterpolation::TestInterpolate () {

  // do this primitive: write a textfile

  std::ostringstream test_file_name;
  test_file_name << "interpolate_test." << N_TEST << ".txt";

  std::ofstream test_file (test_file_name.str ().c_str ());
  test_file << "# test the interpolation" << std::endl;
  test_file << "# X Y LD1" << std::endl;

  float x;
  float y;

  // write the values known from CORSIKA
  for (f_table_t::const_iterator i_nkg = f_nkg.begin ();
       i_nkg != f_nkg.end ();
       i_nkg++) {

    index2cartesian (i_nkg->first, x, y);
    test_file << (i_nkg->first)[0]
	      << " " << (i_nkg->first)[1]
	      << " " << x
	      << " " << y
	      << " " << exp (i_nkg->second) << std::endl;
    
  
  }

  I3Particle zero;
  zero.SetPos (I3Position (0., 0., 0.));

  // interpolate along the the positive x-axis
  for (float x = R_0; x <= 5. * R_MAX; x += 0.5) {

    I3Position pos (x, 0.,  0.);

    test_file << -1
	      << " " << -1
	      << " " << pos.GetX ()
	      << " " << pos.GetY ()
	      << " " << Interpolate (pos, zero) << std::endl;
    

  }

  // interpolate along the the positive y-axis
  for (float y = R_0; y <= 5. * R_MAX; y += 0.5) {

    I3Position pos (0., y,  0.);

    test_file << -1
	      << " " << -1
	      << " " << pos.GetX ()
	      << " " << pos.GetY ()
	      << " " << Interpolate (pos, zero) << std::endl;
    

  }

  // interpolate along the the negative x-axis
  for (float x = R_0; x <= 5. * R_MAX; x += 0.5) {

    I3Position pos (-x, 0.,  0.);

    test_file << -1
	      << " " << -1
	      << " " << pos.GetX ()
	      << " " << pos.GetY ()
	      << " " << Interpolate (pos, zero) << std::endl;
    

  }

  // interpolate along the the negative y-axis
  for (float y = R_0; y <= 5. * R_MAX; y += 0.5) {

    I3Position pos (0., -y,  0.);

    test_file << -1
	      << " " << -1
	      << " " << pos.GetX ()
	      << " " << pos.GetY ()
	      << " " << Interpolate (pos, zero) << std::endl;
    

  }



  // interpolate in between the angles
  for (float phi = 0.; phi <= 2. * M_PI; phi += M_PI / 20) {
    for (float r = R_0; r <= 5. * R_MAX; r += 0.5) {

      x = r * cos (phi);
      y = r * sin (phi);
      
      I3Position pos (x, y,  0.);
      
      test_file << -2
		<< " " << -2
		<< " " << pos.GetX ()
		<< " " << pos.GetY ()
		<< " " << Interpolate (pos, zero) << std::endl;
      
      
    }
    
  }


  

  N_TEST++;
}

/*******************************************************************/


// fill the grid to enable the interpolations. Calculate the first and
// second derivatives as well
void I3MCNKGInterpolation::fillGrid (const I3MCNKGInfoList &nkg_values,
				     const I3Particle &primary) {

  // we store the logarithm (natural, to avoid the non-existing 10^)
  // of the function value as a function of log (r/r_0) and phi.

  // Don't forget to subtract the core position!

  // r_0 is the smallest value, usually the negative of the first nkg
  // values. The largest value is at the tenth position

  for (I3MCNKGInfoList::const_iterator i_nkg = nkg_values.begin ();
       i_nkg != nkg_values.end ();
       i_nkg++) {

    //    static int i (0);
//    fprintf (stderr, "\n\n\nX/Y[%d]=%6.2f/%6.2f\n", 
//	     i++,
//	     i_nkg->Position.GetX () - primary.GetPos ().GetX (),
//	     i_nkg->Position.GetY () - primary.GetPos ().GetY ());
  }

  R_0 = -nkg_values[0].Position.GetX ()
    + primary.GetPos ().GetX ();
  R_MAX = -nkg_values[N_R-1].Position.GetX ()
    + primary.GetPos ().GetX ();

  log_debug ("R_0/R_MAX=%6.2f/%6.2f", R_0, R_MAX);
//  fprintf (stderr, "\n\n\nnkg.X[0]/m/nkg.X[N_R-1]/m=%6.2f/%6.2f\n",
//	   nkg_values[0].Position.GetX () / I3Units::m,
//	   nkg_values[N_R-1].Position.GetX () / I3Units::m);
  //fprintf (stderr, "\n\n\nR_0/R_MAX=%6.2f/%6.2f\n", R_0, R_MAX);

  // need to make r grid equidistant, thus use logarithm of r
  double lr_0 = log (R_0 / R_0);
  double lr_max = log (R_MAX / R_0);
  H_LR = (lr_max - lr_0) / (N_R - 1);

  // here again, we need to hardcode by the values of atan2 (y, x)
  // which is +pi maximally and -3/4 pi minimally
  PHI_0 = atan2 (-1., -1.); // -3/4 pi
  PHI_MAX = atan2 (0., -1.); // pi
  H_PHI = (PHI_MAX - PHI_0) / (N_PHI - 1);

  log_debug ("PHI_0/PHI_MAX=%6.2f/%6.2f", PHI_0, PHI_MAX);
  //fprintf (stderr, "PHI_0/PHI_MAX=%6.2f/%6.2f\n", PHI_0, PHI_MAX);

  log_debug ("H_LR/H_PHI=%6.2f/%6.2f", H_LR, H_PHI);
  //fprintf (stderr, "H_LR/H_PHI=%6.2f/%6.2f\n", H_LR, H_PHI);

  double x; // cartesian coordinates
  double y;
  double r; // polar coordinates
  double phi;
  double l; // logarithm of the density

  for (I3MCNKGInfoList::const_iterator i_nkg = nkg_values.begin ();
       i_nkg != nkg_values.end ();
       i_nkg++) {

    l = log (i_nkg->LD1);
    x = i_nkg->Position.GetX () - primary.GetPos ().GetX ();
    y = i_nkg->Position.GetY () - primary.GetPos ().GetY ();

    r = sqrt (x * x + y * y);
    phi = atan2 (y, x);

    f_nkg[polar2index (r, phi, 0.5)] = l;

    log_debug ("f_nkg=%6.2f", f_nkg[polar2index (r, phi, 0.5)]);
    //fprintf (stderr,
	//     "x/y/r/phi/f_nkg=%6.2f/%6.2f/%6.2f/%6.2f/%6.2f/%6.2f\n",
	//     x, y, r, phi,
	//     f_nkg[polar2index (r, phi, 0.5)],
	//     exp (f_nkg[polar2index (r, phi, 0.5)])
	//     );

  }

  // now we have the function values ordered by rho = log (r/r_0) and
  // phi. Let's start calculating the first derivatives. Take care of
  // the cyclic nature of the phi coordinate and the edges of the r
  // coordinate.
  // if not at the border:
  // df/dr = 0.5 * (f (x+h, phi) - f (x-h, phi)) / h
  // df/dphi = 0.5 * (f (x, phi+h) - f (x, phi-h)) / h
  // if at the upper border: (i.e. no f(.+h) available)
  // df/dr = (f(x, phi) - f(x-h, phi)) / h
  // df/dphi = (f(r, phi) - f(r, phi-h)) / h
  
  for (f_table_t::const_iterator i_nkg = f_nkg.begin ();
       i_nkg != f_nkg.end ();
       i_nkg++) {

    f_table_index_t lr1 (i_nkg->first);
    f_table_index_t lr2 (i_nkg->first);

    f_table_index_t phi1 (i_nkg->first);
    f_table_index_t phi2 (i_nkg->first);

    //fprintf (stderr, "i_nkg->first: index for df_dlr: %02d/%02d\n",
//	     (i_nkg->first)[0], (i_nkg->first)[1]);

//    fprintf (stderr, "lr1: index for df_dlr: %02d/%02d\n",
//	     lr1[0], lr1[1]);


    // are we at the lower/upper r border?
    switch ((i_nkg->first)[0]) {
    case 0:
      lr1[0]++;
      df_dlr[i_nkg->first] = (f_nkg[lr1] - i_nkg->second) / H_LR;
      break;
    case N_R - 1:
      lr1[0]--;
      df_dlr[i_nkg->first] = (i_nkg->second - f_nkg[lr1]) / H_LR;
      break;
    default:
      lr1[0]--;
      lr2[0]++;
      df_dlr[i_nkg->first] = 0.5 / H_LR *  (f_nkg[lr2] - f_nkg[lr1]);
      break;
    } // switch

    // for the phi interpolation we can make use of the cyclic nature
    phi1[1] = (phi1[1] + N_PHI - 1) % N_PHI;
    phi2[1] = (phi2[1] + 1) % N_PHI;
    df_dphi[i_nkg->first] = 0.5 / H_PHI * (f_nkg[phi2] - f_nkg[phi1]);

    log_debug ("df_dlr=%6.2f", df_dlr[i_nkg->first]);
    log_debug ("df_dphi=%6.2f", df_dphi[i_nkg->first]);
    //fprintf (stderr, "df_dlr=%6.2f\n", df_dlr[i_nkg->first]);
    //fprintf (stderr, "df_dphi=%6.2f\n", df_dphi[i_nkg->first]);

  } // loop to fill the first order derivatives

  // now the second order derivatives
  // no values exist for the border values of lr, so set them to 0 in
  // order to ignore them when interpolating.
  for (f_table_t::const_iterator i_nkg = f_nkg.begin ();
       i_nkg != f_nkg.end ();
       i_nkg++) {

    f_table_index_t lr1 (i_nkg->first);
    f_table_index_t lr2 (i_nkg->first);

    f_table_index_t phi1 (i_nkg->first);
    f_table_index_t phi2 (i_nkg->first);

    // are we at the lower/upper r border?
    switch ((i_nkg->first)[0]) {
    case 0:
    case N_R - 1:
      d2f_dlr_dlr[i_nkg->first] = 0.;
      break;
    default:
      lr1[0]--;
      lr2[0]++;
      d2f_dlr_dlr[i_nkg->first] =
	(f_nkg[lr2] - 2. * i_nkg->second + f_nkg[lr1]) / (H_LR * H_LR);
      break;
    } // switch

    // for the phi interpolation we can make use of the cyclic nature
    phi1[1] = (phi1[1] + N_PHI - 1) % N_PHI;
    phi2[1] = (phi2[1] + 1) % N_PHI;
    d2f_dphi_dphi[i_nkg->first] =
      (f_nkg[phi2] - 2. * i_nkg->second + f_nkg[phi1]) / (H_PHI * H_PHI);

    log_debug ("d2f_dlr_dlr=%6.2f", d2f_dlr_dlr[i_nkg->first]);
    log_debug ("d2f_dphi_dphi=%6.2f", d2f_dphi_dphi[i_nkg->first]);
    //fprintf (stderr, "d2f_dlr_dlr=%6.2f\n", d2f_dlr_dlr[i_nkg->first]);
    //fprintf (stderr, "d2f_dphi_dphi=%6.2f\n", d2f_dphi_dphi[i_nkg->first]);


  } // loop to fill the first order derivatives

  // finally, the mixed derivative
  // no values exist for the border values of lr, so set them to 0 in
  // order to ignore them when interpolating.
  for (f_table_t::const_iterator i_nkg = f_nkg.begin ();
       i_nkg != f_nkg.end ();
       i_nkg++) {

    f_table_index_t philr1 (i_nkg->first);
    f_table_index_t philr2 (i_nkg->first);
    f_table_index_t philr3 (i_nkg->first);
    f_table_index_t philr4 (i_nkg->first);


    // are we at the lower/upper r border?
    switch ((i_nkg->first)[0]) {
    case 0:
    case N_R - 1:
      d2f_dlr_dphi[i_nkg->first] = 0.;
      break;
    default:
      philr1[0]++;
      philr1[1] = ( philr1[1] + 1 ) % N_PHI;
      philr2[0]--;
      philr2[1] = ( philr2[1] + 1 ) % N_PHI;
      philr3[0]++;
      philr3[1] = ( philr3[1] + N_PHI - 1 ) % N_PHI;
      philr4[0]--;
      philr4[1] = ( philr4[1] + N_PHI - 1 ) % N_PHI;
      d2f_dlr_dphi[i_nkg->first] =
	(f_nkg[philr1] - f_nkg[philr2] - f_nkg[philr3] + f_nkg[philr4])
	/ (H_LR * H_PHI);
      break;
    } // switch

    log_debug ("d2f_dlr_dphi=%6.2f", d2f_dlr_dphi[i_nkg->first]);
    //fprintf (stderr, "d2f_dlr_dphi=%6.2f\n", d2f_dlr_dphi[i_nkg->first]);

  } // loop to fill the first order derivatives



}

  
/*******************************************************************/

const I3MCNKGInterpolation::f_table_index_t
I3MCNKGInterpolation::polar2index (const double r,
				   const double phi,
				   const double roundoff) const {

  f_table_index_t coords (2);

  //fprintf (stderr, "normalised log(r)=%6.2f, phi=%6.2f\n",
	//   (log (r / R_0) - log (R_0 / R_0))
	//   / (log (R_MAX / R_0) - log (R_0 / R_0))
	//   );
  
  // this makes the equidistant and eases interpolation
  coords[0] = int ((N_R - 1) * (log (r / R_0) - log (R_0 / R_0))
		   / (log (R_MAX / R_0) - log (R_0 / R_0)) + roundoff);
  coords[1] = int ((N_PHI - 1) * (phi - PHI_0) / (PHI_MAX - PHI_0) + roundoff);


  // make sure that the coordinates don't exceed the grid:
  if (coords[0] >= N_R)
    coords[0] = N_R - 1;

  if (coords[1] >= N_PHI)
    coords[1] = N_PHI - 1;

  //fprintf (stderr, "normal phi %f (%d)\n",
	//   (N_PHI - 1.) * (phi - PHI_0) / (PHI_MAX - PHI_0),
	//   int ((N_PHI - 1.) * (phi - PHI_0) / (PHI_MAX - PHI_0) + roundoff)
	//   );
  //fprintf (stderr, "r/phi -> i/j: %6.2f/%6.2f -> %02d/%02d\n",
	//   r, phi, coords[0], coords[1]);

  return coords;

}

/*******************************************************************/

void I3MCNKGInterpolation::index2cartesian (const f_table_index_t &index,
					    float &x, float &y) const {

  // calculate the log of r first, then 
  float r  = float (index[0]) / float (N_R - 1)
    * (log (R_MAX / R_0) - log (R_0 / R_0)) + log (R_0 / R_0);
  r = R_0 * exp (r);

  float phi = float (index[1]) /  float (N_PHI - 1)
    * (PHI_MAX - PHI_0) + PHI_0;

  x = r * cos (phi);
  y = r * sin (phi);

}
