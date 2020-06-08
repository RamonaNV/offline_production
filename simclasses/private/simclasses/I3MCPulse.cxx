#include <simclasses/I3MCPulse.h>
#include <ostream>
#include "I3SequenceOpOStream.h"

I3_SERIALIZABLE(I3MCPulseSeriesMap);

std::ostream& operator<<(std::ostream& os, const I3MCPulse& pulse) {
  return(pulse.Print(os));
}

std::ostream& I3MCPulse::Print(std::ostream& os) const{
    os << "[ I3MCPulse::"
       << "\n  Time   :" << time 
       << "\n  Charge :" << charge
       << "\n  Source :" << source << " ]";
    return os;
}

I3_SEQUENCE_OP_OSTREAM(I3MCPulseSeries, "\n ");
