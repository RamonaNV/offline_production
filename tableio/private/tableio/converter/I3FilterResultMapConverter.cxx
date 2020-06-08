/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3FilterResultMapConverter.cxx 94948 2012-11-04 16:21:52Z nwhitehorn $
 *
 * @version $Revision: 94948 $
 * @date $LastChangedDate: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: nwhitehorn $
 */

#include "tableio/converter/I3FilterResultMapConverter.h"

I3TableRowDescriptionPtr I3FilterResultMapConverter::CreateDescription(const I3FilterResultMap& frmap) {
    I3TableRowDescriptionPtr desc = 
        I3TableRowDescriptionPtr(new I3TableRowDescription() );
    I3FilterResultMap::const_iterator it;
    for (it = frmap.begin(); it != frmap.end(); it++) {
        desc->AddField<bool>(it->first,"bool","Field 0: condition passed, Field 1: prescale passed",2);
    }
    return desc;
}
        
size_t I3FilterResultMapConverter::FillRows(const I3FilterResultMap& frmap, I3TableRowPtr rows) {
    I3FilterResultMap::const_iterator it;
    bool* filter_result;
    for (it = frmap.begin(); it != frmap.end(); it++) {
        filter_result = rows->GetPointer<bool>(it->first);
        filter_result[0] = it->second.conditionPassed;
        filter_result[1] = it->second.prescalePassed;
    }
    return 1;
}
