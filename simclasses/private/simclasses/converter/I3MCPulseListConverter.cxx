/*
 * copyright  (C) 2014
 * The Icecube Collaboration
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Joao Pedro Athayde Marcondes de Andre <jpa14@psu.edu> Last changed by: $LastChangedBy$
 *
 * Based on tpx/private/tpx/converter/convert_I3IceTopBaseline.cxx
 */

#include "I3MCPulseListConverter.h"
#include <icetray/I3Units.h>
#include <tableio/converter/I3MapConverter.h>

namespace convert {
	void I3MCPulseList::AddFields(I3TableRowDescriptionPtr desc, const booked_type&){
		MAKE_ENUM_VECTOR(source, ::I3MCPulse, PulseSource, I3MCPULSE_H_I3MCPulse_PulseSource);
		desc->AddEnumField< ::I3MCPulse::PulseSource >("source", source, "", "Process that originated the pulse");
		desc->AddField<float>("time", "ns", "time");
		desc->AddField<float>("charge", "pe", "charge");
	}

	void I3MCPulseList::FillSingleRow(const booked_type& bl, I3TableRowPtr row){
		row->Set< ::I3MCPulse::PulseSource >("source", bl.source);
		row->Set<float>("time", bl.time/I3Units::ns);
		row->Set<float>("charge", bl.charge);
	}

}
