/**
 *  $Id: I3SumGenerator.cxx 108207 2013-07-13 15:19:02Z nwhitehorn $
 *  
 *  Copyright (C) 2007, 2008
 *  The IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */

#include <sim-services/I3SumGenerator.h>

using namespace boost::python;

struct thunk_to_python
{
  object callable;
  thunk_to_python(object callable_) : callable(callable_) { }

  double operator()(double d)
  {
    object r_obj = callable(d);
    extract<double> extractor(r_obj);
    if (!extractor.check())
      log_fatal("Function passed to I3SumGenerator returned a value not convertible to 'double'");
  
    double r_dub = extractor();
    return r_dub;
  }
};


struct I3SumGeneratorProxy
{
  I3SumGeneratorPtr impl;

  I3SumGeneratorProxy(I3RandomServicePtr r,
		      object fun,
		      const double &xlo, const double &xhi, const int &nbins, 
		      const int &switchgauss, const double &PLow, const int &nBinsLow, 
		      const double &PHigh, const int &nBinsHigh)
    : impl(new I3SumGenerator)
  {
    impl->Initialise(r, thunk_to_python(fun),
		     xlo,xhi,nbins,switchgauss,
		     PLow,nBinsLow,PHigh,nBinsHigh);
  }

  double Generate(int terms) { return impl->Generate(terms); }
};

void register_I3SumGenerator()
{
  class_<I3SumGeneratorProxy>("I3SumGenerator",
			      init<I3RandomServicePtr, object, 
			      double, double, int, int, double, int, double, int>())
    .def("Generate", &I3SumGeneratorProxy::Generate)
    ;

}

