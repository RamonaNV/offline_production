
__kernel void testKernel(__global float* randomNumbers,
                         __global ulong* MWC_RNG_x,
                         __global uint* MWC_RNG_a,
                         const uint numIterations
                         )
{
    //dbg_printf("Start kernel... (work item %u of %u)\n", get_global_id(0), get_global_size(0));

    unsigned int i = get_global_id(0);
    //unsigned int global_size = get_global_size(0);

    //download MWC RNG state
    ulong real_rnd_x = MWC_RNG_x[i];
    uint real_rnd_a = MWC_RNG_a[i];
    uint real_rnd_count = 0;
    ulong *rnd_x = &real_rnd_x;
    uint *rnd_a = &real_rnd_a;
    uint *rnd_count = &real_rnd_count;
#ifdef CLSIM_RNG_TYPE
    clsim_rng_key_t real_rnd_k = {{0}};
    clsim_rng_ctr_t real_rnd_c = {{0,i}};
    clsim_rng_key_t *rnd_k = &real_rnd_k;
    clsim_rng_ctr_t *rnd_c = &real_rnd_c;
#endif

    for (uint j=0; j < numIterations; j++) {
#ifdef CLSIM_RNG_TYPE
        real_rnd_c.v[0] = j;
#endif
     randomNumbers[i] = RNG_CALL_UNIFORM_CO;
    }
    //dbg_printf("Stop kernel... (work item %u of %u)\n", i, global_size);
    //dbg_printf("Kernel finished.\n");

    //upload MWC RNG state
    MWC_RNG_x[i] = real_rnd_x;
    MWC_RNG_a[i] = real_rnd_a;
}
