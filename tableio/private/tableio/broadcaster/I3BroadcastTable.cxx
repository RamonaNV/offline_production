/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3BroadcastTable.cxx 94948 2012-11-04 16:21:52Z nwhitehorn $
 *
 * @version $Revision: 94948 $
 * @date $LastChangedDate: 2012-11-04 09:21:52 -0700 (Sun, 04 Nov 2012) $
 * @author Jakob van Santen <vansanten@wisc.edu> Last changed by: $LastChangedBy: nwhitehorn $
 */

#include <tableio/broadcaster/I3BroadcastTable.h>

I3BroadcastTable::I3BroadcastTable(I3TableService& service, std::string name,
        I3TableRowDescriptionConstPtr description, std::vector<I3TablePtr>& clients)
  : I3Table(service,name,description), clients_(clients)
{};

void I3BroadcastTable::WriteRows(I3TableRowConstPtr rows) {
    log_fatal("I should never have been called!");
};
 
void I3BroadcastTable::AddRow(I3EventHeaderConstPtr header, I3TableRowConstPtr row) {
    std::vector<I3TablePtr>::iterator iter;
    for(iter = clients_.begin(); iter != clients_.end(); ++iter ) {
        (*iter)->AddRow(header,row);
    }
};

void I3BroadcastTable::Align() {
    std::vector<I3TablePtr>::iterator iter;
    for(iter = clients_.begin(); iter != clients_.end(); ++iter ) {
        (*iter)->Align();
    }
};
