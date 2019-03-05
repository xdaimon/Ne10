/*
 *  Copyright 2014-16 ARM Limited and Contributors.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the name of ARM Limited nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY ARM LIMITED AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL ARM LIMITED AND CONTRIBUTORS BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* license of Kiss FFT */
/*
Copyright (c) 2003-2010, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * NE10 Library : dsp/NE10_rfft_float32.c
 */

#include "NE10_types.h"
#include "NE10_macros.h"
#include "NE10_fft.h"
#include "NE10_dsp.h"
#include <math.h>
ne10_fft_r2c_cfg_float32_t ne10_fft_alloc_r2c_float32 (ne10_int32_t nfft)
{
    ne10_fft_r2c_cfg_float32_t st = NULL;
    ne10_int32_t ncfft = nfft >> 1;
    ne10_int32_t result;

    ne10_uint32_t memneeded =   sizeof (ne10_fft_r2c_state_float32_t)
                              + sizeof (ne10_fft_cpx_float32_t) * nfft              /* buffer*/
                              + sizeof (ne10_int32_t) * (NE10_MAXFACTORS * 2)       /* r_factors */
                              + sizeof (ne10_int32_t) * (NE10_MAXFACTORS * 2)       /* r_factors_neon */
                              + sizeof (ne10_fft_cpx_float32_t) * nfft              /* r_twiddles */
                              + sizeof (ne10_fft_cpx_float32_t) * nfft/4            /* r_twiddles_neon */
                              + sizeof (ne10_fft_cpx_float32_t) * (12 + nfft/32*12) /* r_super_twiddles_neon */
                              + NE10_FFT_BYTE_ALIGNMENT;     /* 64-bit alignment*/

    st = (ne10_fft_r2c_cfg_float32_t) NE10_MALLOC (memneeded);

    if (!st)
    {
        return st;
    }

    ne10_int32_t i,j;
    ne10_fft_cpx_float32_t *tw;
    const ne10_float32_t pi = NE10_PI;
    ne10_float32_t phase1;

    st->nfft = nfft;

    uintptr_t address = (uintptr_t) st + sizeof (ne10_fft_r2c_state_float32_t);
    NE10_BYTE_ALIGNMENT (address, NE10_FFT_BYTE_ALIGNMENT);

    st->buffer = (ne10_fft_cpx_float32_t*) address;
    st->r_twiddles = st->buffer + nfft;
    st->r_factors = (ne10_int32_t*) (st->r_twiddles + nfft);
    st->r_twiddles_neon = (ne10_fft_cpx_float32_t*) (st->r_factors + (NE10_MAXFACTORS * 2));
    st->r_factors_neon = (ne10_int32_t*) (st->r_twiddles_neon + nfft/4);
    st->r_super_twiddles_neon = (ne10_fft_cpx_float32_t*) (st->r_factors_neon + (NE10_MAXFACTORS * 2));

    if (nfft<16)
    {
        return st;
    }

    // factors and twiddles for rfft C
    ne10_factor (nfft, st->r_factors, NE10_FACTOR_EIGHT_FIRST_STAGE);

    // backward twiddles pointers
    st->r_twiddles_backward = ne10_fft_generate_twiddles_float32 (st->r_twiddles, st->r_factors, nfft);

    // factors and twiddles for rfft neon
    result = ne10_factor (nfft/4, st->r_factors_neon, NE10_FACTOR_EIGHT_FIRST_STAGE);
    if (result == NE10_ERR)
    {
        return st;
    }

    // Twiddle table is transposed here to improve cache access performance.
    st->r_twiddles_neon_backward = ne10_fft_generate_twiddles_transposed_float32 (
        st->r_twiddles_neon,
        st->r_factors_neon,
        nfft/4);

    // nfft/4 x 4
    tw = st->r_super_twiddles_neon;
    for (i = 1; i < 4; i ++)
    {
        for (j = 0; j < 4; j++)
        {
            phase1 = - 2 * pi * ( (ne10_float32_t) (i * j) / nfft);
            tw[4*i-4+j].r = (ne10_float32_t) cos (phase1);
            tw[4*i-4+j].i = (ne10_float32_t) sin (phase1);
        }
    }

    ne10_int32_t k,s;
    // [nfft/32] x [3] x [4]
    //     k        s     j
    for (k=1; k<nfft/32; k++)
    {
        // transposed
        for (s = 1; s < 4; s++)
        {
            for (j = 0; j < 4; j++)
            {
                phase1 = - 2 * pi * ( (ne10_float32_t) ((k*4+j) * s) / nfft);
                tw[12*k+j+4*(s-1)].r = (ne10_float32_t) cos (phase1);
                tw[12*k+j+4*(s-1)].i = (ne10_float32_t) sin (phase1);
            }
        }
    }
    return st;
}
