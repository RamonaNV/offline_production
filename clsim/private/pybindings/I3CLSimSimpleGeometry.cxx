/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3CLSimSimpleGeometry.cxx 178169 2019-12-19 21:40:42Z jvansanten $
 *
 * @file I3CLSimSimpleGeometry.cxx
 * @version $Revision: 178169 $
 * @date $Date: 2019-12-19 14:40:42 -0700 (Thu, 19 Dec 2019) $
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/I3CLSimSimpleGeometry.h>
#include <clsim/I3CLSimSimpleGeometryTextFile.h>
#include <clsim/I3CLSimSimpleGeometryFromI3Geometry.h>

#include <boost/preprocessor/seq.hpp>

#include "python_gil_holder.h"

using namespace boost::python;
namespace bp = boost::python;

struct I3CLSimSimpleGeometryWrapper : I3CLSimSimpleGeometry, bp::wrapper<I3CLSimSimpleGeometry>
{
    // pure virtual
    virtual std::size_t size() const {utils::python_gil_holder gil; return this->get_override("size")();}
    virtual double GetOMRadius() const {utils::python_gil_holder gil; return this->get_override("GetOMRadius")();}
    
    virtual const std::vector<int32_t> &GetStringIDVector() const {utils::python_gil_holder gil; return this->get_override("GetStringIDVector")();}
    virtual const std::vector<uint32_t> &GetDomIDVector() const {utils::python_gil_holder gil; return this->get_override("GetDomIDVector")();}
    virtual const std::vector<double> &GetPosXVector() const {utils::python_gil_holder gil; return this->get_override("GetPosXVector")();}
    virtual const std::vector<double> &GetPosYVector() const {utils::python_gil_holder gil; return this->get_override("GetPosYVector")();}
    virtual const std::vector<double> &GetPosZVector() const {utils::python_gil_holder gil; return this->get_override("GetPosZVector")();}
    virtual const std::vector<std::string> &GetSubdetectorVector() const  {utils::python_gil_holder gil; return this->get_override("GetSubdetectorVector")();}
    
    virtual int32_t GetStringID(std::size_t pos) const {utils::python_gil_holder gil; return this->get_override("GetStringID")();}
    virtual uint32_t GetDomID(std::size_t pos) const {utils::python_gil_holder gil; return this->get_override("GetDomID")();}
    virtual double GetPosX(std::size_t pos) const {utils::python_gil_holder gil; return this->get_override("GetPosX")();}
    virtual double GetPosY(std::size_t pos) const {utils::python_gil_holder gil; return this->get_override("GetPosY")();}
    virtual double GetPosZ(std::size_t pos) const {utils::python_gil_holder gil; return this->get_override("GetPosZ")();}
    virtual std::string GetSubdetector(std::size_t pos) const {utils::python_gil_holder gil; return this->get_override("GetSubdetector")();}

};

static boost::shared_ptr<I3CLSimSimpleGeometry>
MakeSimpleGeometrySimply(double OMRadius, double oversizeFactor,
                         const I3FramePtr &frame,
                         const std::vector<int> &ignoreStrings,
                         const std::vector<unsigned int> &ignoreDomIDs,
                         const std::vector<std::string> &ignoreSubdetectors,
                         int32_t ignoreStringIDsSmallerThan,
                         int32_t ignoreStringIDsLargerThan,
                         uint32_t ignoreDomIDsSmallerThan,
                         uint32_t ignoreDomIDsLargerThan,
                         bool splitIntoPartsAccordingToPosition,
                         bool useHardcodedDeepCoreSubdetector)
{
	std::set<int> ignoreStringss(ignoreStrings.begin(), ignoreStrings.end());
	std::set<unsigned int> ignoreDomIDss(ignoreDomIDs.begin(), ignoreDomIDs.end());
	std::set<std::string> ignoreSubdetectorss(ignoreSubdetectors.begin(), ignoreSubdetectors.end());
	return boost::make_shared<I3CLSimSimpleGeometry>(I3CLSimSimpleGeometryFromI3Geometry(
		OMRadius, oversizeFactor,
		frame,
		ignoreStringss,
		ignoreDomIDss,
		ignoreSubdetectorss,
		ignoreStringIDsSmallerThan,
		ignoreStringIDsLargerThan,
		ignoreDomIDsSmallerThan,
		ignoreDomIDsLargerThan,
		splitIntoPartsAccordingToPosition,
		useHardcodedDeepCoreSubdetector));
}

void register_I3CLSimSimpleGeometry()
{
    {
        bp::scope I3CLSimSimpleGeometry_scope = 
        bp::class_<I3CLSimSimpleGeometry, boost::shared_ptr<I3CLSimSimpleGeometry>, boost::noncopyable>("I3CLSimSimpleGeometry", bp::init<double,double>((bp::arg("DOMRadius"), "CableRadius")))
        .def("size", bp::pure_virtual(&I3CLSimSimpleGeometry::size))
        .def("__len__", bp::pure_virtual(&I3CLSimSimpleGeometry::size))
        
        .def("GetOMRadius", bp::pure_virtual(&I3CLSimSimpleGeometry::GetOMRadius))

        .def("GetStringIDVector", &I3CLSimSimpleGeometry::GetStringIDVector, bp::return_value_policy<bp::copy_const_reference>())
        .def("GetDomIDVector", &I3CLSimSimpleGeometry::GetDomIDVector, bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosXVector", &I3CLSimSimpleGeometry::GetPosXVector, bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosYVector", &I3CLSimSimpleGeometry::GetPosYVector, bp::return_value_policy<bp::copy_const_reference>())
        .def("GetPosZVector", &I3CLSimSimpleGeometry::GetPosZVector, bp::return_value_policy<bp::copy_const_reference>())
        .def("GetSubdetectorVector", &I3CLSimSimpleGeometry::GetSubdetectorVector, bp::return_value_policy<bp::copy_const_reference>())

        .add_property("stringIDs", bp::make_function(&I3CLSimSimpleGeometry::GetStringIDVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("domIDs", bp::make_function(&I3CLSimSimpleGeometry::GetDomIDVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("posX", bp::make_function(&I3CLSimSimpleGeometry::GetPosXVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("posY", bp::make_function(&I3CLSimSimpleGeometry::GetPosYVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("posZ", bp::make_function(&I3CLSimSimpleGeometry::GetPosZVector, bp::return_value_policy<bp::copy_const_reference>()))
        .add_property("subdetectors", bp::make_function(&I3CLSimSimpleGeometry::GetSubdetectorVector, bp::return_value_policy<bp::copy_const_reference>()))
        
        .def("GetStringID", &I3CLSimSimpleGeometry::GetStringID)
        .def("GetDomID", &I3CLSimSimpleGeometry::GetDomID)
        .def("GetPosX", &I3CLSimSimpleGeometry::GetPosX)
        .def("GetPosY", &I3CLSimSimpleGeometry::GetPosY)
        .def("GetPosZ", &I3CLSimSimpleGeometry::GetPosZ)
        .def("GetSubdetector", &I3CLSimSimpleGeometry::GetSubdetector)

        .def("AddModule", &I3CLSimSimpleGeometry::AddModule, (bp::arg("string"), "dom", "x", "y", "z", "subdetector", bp::arg("cableAngle")=NAN))
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3CLSimSimpleGeometry>, boost::shared_ptr<const I3CLSimSimpleGeometry> >();
    
    // I3CLSimSimpleGeometryTextFile
    {
        bp::def("I3CLSimSimpleGeometryTextFile",
            &I3CLSimSimpleGeometryTextFile,
            (
            bp::arg("OMRadius"),
            bp::arg("filename"),
            bp::arg("ignoreStringIDsSmallerThan")=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreStringIDsSmallerThan,
            bp::arg("ignoreStringIDsLargerThan")=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreStringIDsLargerThan,
            bp::arg("ignoreDomIDsSmallerThan")=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreDomIDsSmallerThan,
            bp::arg("ignoreDomIDsLargerThan")=clsim::defaults::I3CLSimSimpleGeometryTextFile::ignoreDomIDsLargerThan
            )
        );
    }

    // I3CLSimSimpleGeometryFromI3Geometry
    {
        bp::def("I3CLSimSimpleGeometryFromI3Geometry",
            &MakeSimpleGeometrySimply,
            (
            bp::arg("OMRadius"),
            bp::arg("oversizeFactor"),
            bp::arg("frame"),
#define default(name, type) std::vector<type>(clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::name.begin(), clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::name.end()) 
            bp::arg("ignoreStrings")=default(ignoreStrings, int),
            bp::arg("ignoreDomIDs")=default(ignoreDomIDs, unsigned),
            bp::arg("ignoreSubdetectors")=default(ignoreSubdetectors, std::string),
#undef default
            bp::arg("ignoreStringIDsSmallerThan")=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreStringIDsSmallerThan,
            bp::arg("ignoreStringIDsLargerThan")=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreStringIDsLargerThan,
            bp::arg("ignoreDomIDsSmallerThan")=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDsSmallerThan,
            bp::arg("ignoreDomIDsLargerThan")=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::ignoreDomIDsLargerThan,
            bp::arg("splitIntoPartsAccordingToPosition")=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::splitIntoPartsAccordingToPosition,
            bp::arg("useHardcodedDeepCoreSubdetector")=clsim::defaults::I3CLSimSimpleGeometryFromI3Geometry::useHardcodedDeepCoreSubdetector
            )
        );
    }

}
