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
 * $Id: TrkOpticalPhysics.hh 179360 2020-03-10 16:07:35Z eganster $
 *
 * @file TrkOpticalPhysics.hh
 * @version $Revision: 179360 $
 * @date $Date: 2020-03-10 10:07:35 -0600 (Tue, 10 Mar 2020) $
 * @author Claudio Kopper
 */

#ifndef TrkOpticalPhysics_hh
#define TrkOpticalPhysics_hh

#include "G4VPhysicsConstructor.hh"
#include "clsim/function/I3CLSimFunction.h"

class TrkCerenkov;

class TrkOpticalPhysics : public G4VPhysicsConstructor
{
public:
    TrkOpticalPhysics(const G4String& name,
                      double maxBetaChangePerStep,
                      uint32_t maxNumPhotonsPerStep);
    virtual ~TrkOpticalPhysics();
    
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void SetWlenBiasFunction(I3CLSimFunctionConstPtr wlenBias);
protected:
    TrkCerenkov* theCerenkovProcess;
    
private:
    double maxBetaChangePerStep_;
    uint32_t maxNumPhotonsPerStep_;
    I3CLSimFunctionConstPtr wlenBias_;
};

#endif
