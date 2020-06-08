#include <simclasses/I3MCPE.h>
#include <ostream>
#include "I3SequenceOpOStream.h"

I3_SERIALIZABLE(I3MCPESeriesMap);

std::ostream& operator<<(std::ostream& os, const I3MCPE& pe) {
  return(pe.Print(os));
}

std::ostream& I3MCPE::Print(std::ostream& os) const{
  os << "[ I3MCPE::"
     << "\n  Time :" << time
     << "\n  NPE  :" << npe
     << "\n  " << ID
     << " ]";
  return os;
}

I3_SEQUENCE_OP_OSTREAM(I3MCPESeries, "\n ");
