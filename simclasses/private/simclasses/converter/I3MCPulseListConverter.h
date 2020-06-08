#ifndef I3MCPULSELISTCONVERTER_H_INCLUDED
#define I3MCPULSELISTCONVERTER_H_INCLUDED

/*
 * copyright  (C) 2014
 * The Icecube Collaboration
 *
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Joao Pedro Athayde Marcondes de Andre <jpa14@psu.edu> Last changed by: $LastChangedBy$
 *
 * Based on tpx/private/tpx/converter/convert_I3IceTopBaseline.h
 */

#include "simclasses/I3MCPulse.h"
#include <tableio/I3TableRowDescription.h>
#include <tableio/I3TableRow.h>

namespace convert {

	struct I3MCPulseList {
		typedef ::I3MCPulse booked_type;

		void AddFields(I3TableRowDescriptionPtr desc, const booked_type& = booked_type());
		void FillSingleRow(const booked_type& dl, I3TableRowPtr row);
	};

}

#endif  // I3MCPULSELISTCONVERTER_H_INCLUDED
