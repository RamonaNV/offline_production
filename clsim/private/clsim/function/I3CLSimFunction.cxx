/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3CLSimFunction.cxx 179360 2020-03-10 16:07:35Z eganster $
 *
 * @file I3CLSimFunction.cxx
 * @version $Revision: 179360 $
 * @date $Date: 2020-03-10 10:07:35 -0600 (Tue, 10 Mar 2020) $
 * @author Claudio Kopper
 */

#include <icetray/serialization.h>
#include <clsim/function/I3CLSimFunction.h>
#include <dataclasses/I3Map.h>

I3CLSimFunction::I3CLSimFunction()
{ 
    
}

I3CLSimFunction::~I3CLSimFunction() 
{ 

}

template <class Archive>
void I3CLSimFunction::serialize(Archive &ar, unsigned version)
{
    if (version>i3clsimfunction_version_)
        log_fatal("Attempting to read version %u from file but running version %u of I3CLSimFunction class.",version,i3clsimfunction_version_);
    if (version>0)
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
}     


I3_SERIALIZABLE(I3CLSimFunction);
typedef I3Map<OMKey, I3CLSimFunctionConstPtr> I3CLSimFunctionPtrMap;
I3_SERIALIZABLE(I3CLSimFunctionPtrMap);
