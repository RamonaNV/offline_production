/*
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3Converter.cxx 85585 2012-02-25 17:29:16Z claudio.kopper $
 *
 * @version $Revision: 85585 $
 * @date $LastChangedDate: 2012-02-25 10:29:16 -0700 (Sat, 25 Feb 2012) $
 * @author Eike Middell <eike.middell@desy.de> Last changed by: $LastChangedBy: claudio.kopper $
 */

#include <icetray/I3FrameObject.h>

#include "tableio/I3Converter.h"

/******************************************************************************/

I3Converter::~I3Converter() {;}

I3TableRowDescriptionConstPtr I3Converter::GetDescription() {
    if (description_)
        return description_;
    else
        log_fatal("description has not been created yet!");
}

/******************************************************************************/
