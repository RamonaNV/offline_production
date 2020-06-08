/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3MetaService.cxx 169291 2012-11-04 16:17:04Z nwhitehorn $
 *
 * @version $Revision: 169291 $
 * @date $LastChangedDate: 2012-11-04 09:17:04 -0700 (Sun, 04 Nov 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: nwhitehorn $
 */

#include "interfaces/I3MetaService.h"
#include "icetray/I3Frame.h"

#include "install.h"

namespace bp = boost::python;

struct I3MetaServiceWrapper : I3MetaService, bp::wrapper<I3MetaService>
{
  I3FramePtr PopMeta() {  return this->get_override("PopMeta")(); }
};

void
register_I3MetaService()
{	
	bp::class_<I3MetaServiceWrapper, boost::shared_ptr<I3MetaServiceWrapper>,
	    boost::noncopyable>("I3MetaService", bp::init<>())
		.def("pop_meta", &I3MetaServiceWrapper::PopMeta)
		.def("install", &I3InstallService<I3MetaService>().func)
	;
}
