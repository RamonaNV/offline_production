//
//   Copyright (c) 2015   Juan Carlos Diaz Velez and the IceCube Collaboration 
//   
//   This is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 3 of the License, or
//   (at your option) any later version.
//
//   This sofware is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <icetray/python/context_suite.hpp>
#include <sim-services/I3GeneratorService.h>

using namespace boost::python;
namespace bp = boost::python;

struct I3GeneratorServiceWrapper : I3GeneratorService, wrapper<I3GeneratorService>
{
  public:
    double GetRate() { return this->get_override("GetRate")(); }
    I3MCTreePtr GetNextEvent() { return this->get_override("GetNextEvent")(); }
    I3FramePtr GetNextFrame() { return this->get_override("GetNextFrame")(); }
};




template <typename T, typename Init>
class_<T, boost::shared_ptr<T>, boost::python::bases<I3GeneratorService>, boost::noncopyable>
register_coincidenteventservice(const char* name, const char* doc, const Init& init)
{
  return class_<T, boost::shared_ptr<T>, boost::python::bases<I3GeneratorService>, boost::noncopyable>(name,
							     doc,
							     init)
    ;
}


void register_I3GeneratorService()
{
	bp::class_<I3GeneratorServiceWrapper, boost::shared_ptr<I3GeneratorServiceWrapper>, boost::noncopyable>(
	    "I3GeneratorService", "Base class for coincident event services") 
	        .def("get_next", pure_virtual(&I3GeneratorService::GetNextEvent))
	        .def("get_next_frame", pure_virtual(&I3GeneratorService::GetNextFrame))
	        .def("get_rate", pure_virtual(&I3GeneratorService::GetRate))
	;
}
