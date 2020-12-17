#if defined(BENCHMARK_OPENCL)
    constexpr bool shuffleSteps = false;
    constexpr bool shuffleStepsTo32 = false;
    constexpr int resizingNrStepsby2PowerMinus = 0;
    constexpr bool runCudaOnly = false; // only run cuda and not opencl (to speed up big runs)  
    constexpr bool returnPhotonsCUDA = false || runCudaOnly;  // return CUDA phtoons instead of CL phtotons (always true when only cuda is executed)
#elif defined(BENCHMARK_OPENCL_SHUFFLED)
    constexpr bool shuffleSteps = true;
    constexpr bool shuffleStepsTo32 = false;
    constexpr int resizingNrStepsby2PowerMinus = 0;
    constexpr bool runCudaOnly = false; // only run cuda and not opencl (to speed up big runs)  
    constexpr bool returnPhotonsCUDA = false || runCudaOnly;  // return CUDA phtoons instead of CL phtotons (always true when only cuda is executed)
#elif defined(BENCHMARK_OPENCL_SHUFFLED_32)
    constexpr bool shuffleSteps = true;
    constexpr bool shuffleStepsTo32 = true;
    constexpr int resizingNrStepsby2PowerMinus = 0;
    constexpr bool runCudaOnly = false; // only run cuda and not opencl (to speed up big runs)  
    constexpr bool returnPhotonsCUDA = false || runCudaOnly;  // return CUDA phtoons instead of CL phtotons (always true when only cuda is executed)
#elif defined(BENCHMARK_SHUFFLED)
    constexpr bool shuffleSteps = true;
    constexpr bool shuffleStepsTo32 = false;
    constexpr int resizingNrStepsby2PowerMinus = 0;
    constexpr bool runCudaOnly = true; // only run cuda and not opencl (to speed up big runs)  
    constexpr bool returnPhotonsCUDA = true || runCudaOnly;  // return CUDA phtoons instead of CL phtotons (always true when only cuda is executed)
#elif defined(BENCHMARK_SHUFFLED_32)
    constexpr bool shuffleSteps = true;
    constexpr bool shuffleStepsTo32 = true;
    constexpr int resizingNrStepsby2PowerMinus = 0;
    constexpr bool runCudaOnly = true; // only run cuda and not opencl (to speed up big runs)  
    constexpr bool returnPhotonsCUDA = true || runCudaOnly;  // return CUDA phtoons instead of CL phtotons (always true when only cuda is executed)
#else
    constexpr bool shuffleSteps = false;
    constexpr bool shuffleStepsTo32 = false;
    constexpr int resizingNrStepsby2PowerMinus = 0;
    constexpr bool runCudaOnly = true; // only run cuda and not opencl (to speed up big runs)  
    constexpr bool returnPhotonsCUDA = true || runCudaOnly;  // return CUDA phtoons instead of CL phtotons (always true when only cuda is executed)
#endif