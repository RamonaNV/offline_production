/**
 * Copyright (c) 2019
 * Jakob van Santen <jakob.van.santen@desy.de>
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
 * $Id: I3CLSimRandomNumberGeneratorBenchmark.h 171636 2019-02-28 15:01:32Z jvansanten $
 *
 * @file I3CLSimRandomNumberGeneratorBenchmark.h
 * @version $Revision: 171636 $
 * @date $Date: 2019-02-28 08:01:32 -0700 (Thu, 28 Feb 2019) $
 * @author Jakob van Santen
 */

#ifndef I3CLSIMRANDOMNUMBERGENERATORBENCHMARK_H_INCLUDED
#define I3CLSIMRANDOMNUMBERGENERATORBENCHMARK_H_INCLUDED

#include "test/I3CLSimTesterBase.h"

#include "phys-services/I3RandomService.h"

#include "dataclasses/I3Vector.h"
#include "clsim/random_value/I3CLSimRandomValue.h"

class I3CLSimRandomNumberGeneratorBenchmark : public I3CLSimTesterBase
{
public:
    I3CLSimRandomNumberGeneratorBenchmark(const I3CLSimOpenCLDevice &device,
                                    uint64_t workgroupSize_,
                                    uint64_t workItemsPerIteration_,
                                    uint32_t iterations,
                                    I3RandomServicePtr randomService,
                                    std::string rngType="");

    std::pair<I3VectorFloatPtr, uint64_t> GenerateRandomNumbers(uint64_t iterations);
    
private:
    void FillSource(std::vector<std::string> &source,
                    std::string rngType);
    
    void InitBuffers(I3RandomServicePtr randomService, uint32_t iterations);
    
    boost::shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_x;
    boost::shared_ptr<cl::Buffer> deviceBuffer_MWC_RNG_a;

    boost::shared_ptr<cl::Buffer> deviceBuffer_results;

};



#endif //I3CLSIMRANDOMNUMBERGENERATORBENCHMARK_H_INCLUDED
