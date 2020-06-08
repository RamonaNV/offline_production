//
//   Copyright (c) 2009
//   the icecube collaboration
//   $Id: module.cxx 171049 2019-01-30 12:45:22Z kjmeagher $
//                                                                                                       
//   @version $Revision: 171049 $                                                                                                    
//   @author Henrike Wissing 
//

#include <icetray/I3FrameObject.h>
#include <icetray/load_project.h>
#include <sim-services/I3CombineMCPE.h>

using namespace boost::python;
namespace bp = boost::python;

void register_I3SumGenerator();
void register_I3CosmicEventGenerator();
void register_I3GeneratorService();
void register_I3PropagatorService();
void register_I3SimConstants();
void register_ShowerParameters();
void register_I3GeneratorService();

BOOST_PYTHON_MODULE(sim_services)
{
  load_project("sim-services", false);
    
  def("MergeMCPEs", MergeMCPEs);

  register_I3SumGenerator();
  register_I3CosmicEventGenerator();
  register_I3PropagatorService();
  register_I3GeneratorService();
  register_I3SimConstants();
  register_ShowerParameters();
}

