/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3CSVTableService.cxx 94948 2012-11-04 16:21:52Z nwhitehorn $
 *
 * @version $$
 * @date $LastChangedDate: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: nwhitehorn $
 */

#include "tableio/textwriter/I3CSVTableService.h"

namespace bp = boost::python;

void register_I3CSVTableService() {
    bp::class_<I3CSVTableService, 
               boost::shared_ptr<I3CSVTableService>,
               bp::bases<I3TableService> >
               ("I3CSVTableService",
                bp::init<const std::string>(bp::args("folder_name")))
               ;
}
