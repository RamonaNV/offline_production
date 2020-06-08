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
 * $Id: I3CLSimRandomNumberGeneratorBenchmark.cxx 171636 2019-02-28 15:01:32Z jvansanten $
 *
 * @file I3CLSimRandomNumberGeneratorBenchmark.cxx
 * @version $Revision: 171636 $
 * @date $Date: 2019-02-28 08:01:32 -0700 (Thu, 28 Feb 2019) $
 * @author Jakob van Santen
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "test/I3CLSimRandomNumberGeneratorBenchmark.h"

#include <string>

#include "opencl/I3CLSimHelperLoadProgramSource.h"
#include "opencl/mwcrng_init.h"

#include "clsim/I3CLSimHelperToFloatString.h"
using namespace I3CLSimHelper;

#define THREEFRY

I3CLSimRandomNumberGeneratorBenchmark::I3CLSimRandomNumberGeneratorBenchmark
(const I3CLSimOpenCLDevice &device,
 uint64_t workgroupSize_,
 uint64_t workItemsPerIteration_,
 uint32_t iterations,
 I3RandomServicePtr randomService,
 std::string rngType):
I3CLSimTesterBase()
{
    std::vector<std::string> source;
    FillSource(source, rngType);
    
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";
    DoSetup(device,
            workgroupSize_,
            workItemsPerIteration_,
            source,
            rngType.empty() ? "" : " -I "+kernelBaseDir+"/Random123/include -DCLSIM_RNG_TYPE="+rngType);
    
    InitBuffers(randomService,iterations);
}

void I3CLSimRandomNumberGeneratorBenchmark::FillSource(std::vector<std::string> &source,
                                                 std::string rngType)
{
    source.clear();
    
    // load program source from files
    const std::string I3_BUILD(getenv("I3_BUILD"));
    const std::string kernelBaseDir = I3_BUILD+"/clsim/resources/kernels";

    std::string mwcrngSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+(rngType.empty() ? "/mwcrng_kernel.cl" : "/random123_kernel.cl"));
    std::string randomDistSource;

    // std::string rngDistTestKernelHeader = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/rng_dist_test_kernel.h.cl");
    std::string rngDistTestKernelSource = I3CLSimHelper::LoadProgramSource(kernelBaseDir+"/rng_benchmark_kernel.c.cl");
    
    // collect the program sources
    // source.push_back(rngDistTestKernelHeader);
    source.push_back(mwcrngSource);
    source.push_back(randomDistSource);
    source.push_back(rngDistTestKernelSource);
    {
    std::stringstream ssource;
    for(std::string &part : source)
    	ssource << part;
    std::string line;
    unsigned lineno = 1;
    while (std::getline(ssource, line)) {
    	log_debug_stream(std::setw(4) << lineno << " " << line);
    	lineno++;
    }
    }
    
}

void I3CLSimRandomNumberGeneratorBenchmark::InitBuffers(I3RandomServicePtr randomService, uint32_t iterations)
{
    // set up rng
    log_info("Setting up RNG for %" PRIu64 " workitems. (requested workgroupSize=%" PRIu64 ", maxWorkgroupSize=%" PRIu64 ")", workItemsPerIteration, workgroupSize, GetMaxWorkgroupSize());
    
    std::vector<uint64_t> MWC_RNG_x;
    std::vector<uint32_t> MWC_RNG_a;
    
    MWC_RNG_x.resize(workItemsPerIteration);
    MWC_RNG_a.resize(workItemsPerIteration);
    
    if (init_MWC_RNG(&(MWC_RNG_x[0]), &(MWC_RNG_a[0]), workItemsPerIteration, randomService)!=0) 
        throw std::runtime_error("I3CLSimStepToPhotonConverterOpenCL already initialized!");
    log_info("RNG is set up..");
    
    log_debug("Setting up device buffers.");
    // set up device buffers from existing host buffers
    deviceBuffer_MWC_RNG_x = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_x.size() * sizeof(uint64_t), &(MWC_RNG_x[0])));
    deviceBuffer_MWC_RNG_a = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, MWC_RNG_a.size() * sizeof(uint32_t), &(MWC_RNG_a[0])));

    // allocate empty buffers on the device
    deviceBuffer_results = boost::shared_ptr<cl::Buffer>(new cl::Buffer(*context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR, workItemsPerIteration*sizeof(float), NULL));
    log_debug("Device buffers are set up.");

    log_debug("Configuring kernel.");
    {
        kernel->setArg(0, *deviceBuffer_results);          // hit counter
        kernel->setArg(1, *deviceBuffer_MWC_RNG_x);        // rng state
        kernel->setArg(2, *deviceBuffer_MWC_RNG_a);        // rng state
        kernel->setArg(3, iterations);        // rng state
    }
    log_debug("Kernel configured.");
}


std::pair<I3VectorFloatPtr, uint64_t> I3CLSimRandomNumberGeneratorBenchmark::GenerateRandomNumbers(uint64_t repetitions)
{
    // allocate the output vector
    I3VectorFloatPtr results = I3VectorFloatPtr(new I3VectorFloat());
    results->reserve(workItemsPerIteration*repetitions);

    log_info("Starting iterations..");
    
    uint64_t kernelTime(0);
    for (uint64_t i=0;i<repetitions;++i)
    {
        log_debug("Starting kernel..");
        
        // run the kernel
        cl::Event kernelFinishEvent;
        try {
            queue->enqueueNDRangeKernel(*kernel, 
                                       cl::NullRange,    // current implementations force this to be NULL
                                       cl::NDRange(workItemsPerIteration),    // number of work items
                                       cl::NDRange(workgroupSize),
                                       NULL,
                                       &kernelFinishEvent);
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }
        
        try {
            // wait for the kernel
            queue->finish();
        } catch (cl::Error &err) {
            log_fatal("OpenCL ERROR (running kernel): %s (%i)", err.what(), err.err());
        }
        
        log_debug("kernel finished!");
        uint64_t timeStart, timeEnd;
        kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_START, &timeStart);
        kernelFinishEvent.getProfilingInfo(CL_PROFILING_COMMAND_END, &timeEnd);
        kernelTime += (timeEnd-timeStart);

        cl::Event mappingComplete;
        float *mapped_results = (float *)queue->enqueueMapBuffer(*deviceBuffer_results, CL_FALSE, CL_MAP_READ, 0, workItemsPerIteration*sizeof(float), NULL, &mappingComplete);
        
        // wait for the buffer to be mapped
        mappingComplete.wait();
        
        // add the values to the output vector
        for (uint64_t j=0;j<workItemsPerIteration;++j)
        {
            results->push_back(mapped_results[j]);
        }
        
        queue->enqueueUnmapMemObject(*deviceBuffer_results, mapped_results);
    }
    
    log_info_stream(kernelTime << " " << workItemsPerIteration << " " << repetitions);
    log_info("Kernel time %.2f ns per call", kernelTime / double(workItemsPerIteration*repetitions));
    log_info("iterations complete.");
    
    return std::make_pair(results,kernelTime);
}
