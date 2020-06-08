/**
 * copyright  (C) 2004
 * the icecube collaboration
 * $Id: I3MetaService.h 169292 2012-11-04 16:34:42Z nwhitehorn $
 *
 * @file I3MetaService.h
 * @version $Revision: 169292 $
 * @date $Date: 2012-11-04 09:34:42 -0700 (Sun, 04 Nov 2012) $
 * @author pretz
 */

#ifndef INTERFACES_I3METASERVICE_H_INCLUDED
#define INTERFACES_I3METASERVICE_H_INCLUDED

#include <icetray/I3Frame.h>
#include <icetray/I3DefaultName.h>

class I3MetaService
{
 public:

  /**
   * @brief indicates whether or not there are more events to find.
   * 
   */
  virtual I3FramePtr PopMeta() = 0;

  virtual ~I3MetaService();
};

I3_DEFAULT_NAME(I3MetaService);
I3_POINTER_TYPEDEFS(I3MetaService);

#endif
