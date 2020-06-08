/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3WimpParamsConverter.h 66214 2010-08-17 16:39:23Z mzoll $
 *
 * @version $Revision: 66214 $
 * @date $LastChangedDate: 2010-08-17 18:39:23 +0200 (Tue, 17 Aug 2010) $
 * @author Marcel Zoll <marcel.zoll@fysik.su.se> $LastChangedBy: mzoll $
 * 
 * converter class for tableio
 */

#ifndef SIMCLASSES_I3WIMPPARAMSONVERTER_H_INCLUDED
#define SIMCLASSES_I3WIMPPARAMSONVERTER_H_INCLUDED

#include "tableio/I3Converter.h"
#include "simclasses/I3WimpParams.h"

class I3WimpParamsConverter : public I3ConverterImplementation< I3WimpParams > {
public:
  I3WimpParamsConverter();

private:
  I3TableRowDescriptionPtr CreateDescription(const I3WimpParams& wimpparams);
  
  size_t FillRows(const I3WimpParams& wimpparams, I3TableRowPtr rows);
};

#endif // SIMCLASSES_I3WIMPPARAMSONVERTER_H_INCLUDED
