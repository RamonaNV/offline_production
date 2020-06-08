/**
 * tableio pybindings
 *
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: module.cxx 87951 2012-05-04 17:19:19Z olivas $
 *
 * @version $Revision: 87951 $
 * @date $LastChangedDate: 2012-05-04 11:19:19 -0600 (Fri, 04 May 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: olivas $
 */

#include <icetray/load_project.h>

namespace bp = boost::python;

#define REGISTER_THESE_THINGS \
   (I3TableRowDescription)(I3TableRow)(I3Converter)(I3TableService)     \
   (I3TableWriter)(I3TableTranscriber)(I3ConverterBundle)(I3Datatype)   \
   (dataclasses_converters)(I3Table)(I3BroadcastTableService)           \
   (I3CSVTableService)
#define I3_REGISTRATION_FN_DECL(r, data, t) void BOOST_PP_CAT(register_,t)();
#define I3_REGISTER(r, data, t) BOOST_PP_CAT(register_,t)();
BOOST_PP_SEQ_FOR_EACH(I3_REGISTRATION_FN_DECL, ~, REGISTER_THESE_THINGS)

I3_PYTHON_MODULE(tableio)
{

  load_project("tableio", false); 

  BOOST_PP_SEQ_FOR_EACH(I3_REGISTER, ~, REGISTER_THESE_THINGS);
}
