/*The MIT License (MIT)

Copyright (c) 2020, Hendrik Schwanekamp hschwanekamp@nvidia.com, Ramona Hohl rhohl@nvidia.com

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGSEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/* 
    device memory data for scattering and absorbtion look up tables
*/

#ifndef SCATTERINGANDABSORBTIONDATA_CUH
#define SCATTERINGANDABSORBTIONDATA_CUH

#define MEDIUM_LAYERS 171

#define MEDIUM_MIN_WLEN 2.6500000000e-07f
#define MEDIUM_MAX_WLEN 6.7500000000e-07f

#define MEDIUM_MIN_RECIP_WLEN 1.4814814815e+06f
#define MEDIUM_MAX_RECIP_WLEN 3.7735849057e+06f

// medium layer structure:
#define MEDIUM_LAYER_BOTTOM_POS -8.5540000000e+02f
#define MEDIUM_LAYER_THICKNESS 1.0000000000e+01f

#endif