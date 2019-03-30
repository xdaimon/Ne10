#include "NE10_types.h"
#include "NE10_macros.h"
#include "NE10_fft.h"
#include <math.h>
#include <string.h>

/*
 * This function outputs a factor buffer ('facbuf') that decomposes an FFT of input size
 * n into a number of radix-r butterfly calculations (for r in some set of radix values).
 *
 * Factor buffer layout:
 *     index 0: stage count
 *     index 1: stride for the first stage
 *     index 2 to (2*stage_count + 1): pairs of factors (number of sections, section size)
 *     index (2*stage_count + 2): an flag specifying which algorithm to use
 *
 * e.g. 1024 samples might result in the following five stage radix-4 factors buffer:
 *          [5, 256, 4, 256, 4, 64, 4, 16, 4, 4, 4, 1]
 *          i.e. 1024 = 4x256, each of which is 4x64, each of which is 4x16, each of which
 *               is 4x4, each of which is 4x1. There are 5 stages, and the stride for the
 *               first stage is 256 (1024 / 4, for radix-4).
 *
 * Only the leading 42 int32 is used to store factors.
 * The left can be used as algorithm flags, or status flags.
 * Even the leading bits of stage number can be reused.
 */
ne10_int32_t ne10_factor (ne10_int32_t n,
        ne10_int32_t * facbuf,
        ne10_int32_t ne10_factor_flags)
{
    // This is a workaround. We need to "return" some flags.
    // Otherwise, we need to modify signature of ne10_factor.
    assert (NE10_MAXFACTORS >= 32);

    if ((facbuf == NULL)
        || (n <= 0))
    {
        return NE10_ERR;
    }

    ne10_int32_t p;
    ne10_int32_t i = 1;
    ne10_int32_t stage_num = 0;
    ne10_int32_t stride_max = n;

    // Default algorithm flag is NE10_FFT_ALG_DEFAULT
    ne10_int32_t alg_flag = NE10_FFT_ALG_DEFAULT;

    // Factor out powers of 4, 2, 5, and 3. Additionally, factor out powers
    // of 8 if the right factor flags are passed. If none of these factors
    // can be applied at any stage, the remaining size is used as a factor.
    do
    {
        // If NE10_FACTOR_EIGHT_FIRST_STAGE is enabled, we can generate
        // a first stage of radix-8 (e.g. by combining one radix-4 and
        // one radix-2 stage into a single radix-8 stage).
        if ((ne10_factor_flags & NE10_FACTOR_EIGHT_FIRST_STAGE)
                && ((n==8) || (n==40) || (n==24)))
        {
            switch (n)
            {
            case 8:
                p = 8;
                break;
            case 24:
                p = 3;
                alg_flag = NE10_FFT_ALG_ANY;
                break;
            default: // n == 40
                p = 5;
                alg_flag = NE10_FFT_ALG_ANY;
                break;
            }
        }
        else if ((ne10_factor_flags & NE10_FACTOR_EIGHT) && ((n % 8) == 0))
        {
            p = 8;
        }
        else if ((n % 4) == 0)
        {
            p = 4;
        }
        else if ((n % 2) == 0)
        {
            p = 2;
        }
        else if ((n % 5) == 0)
        {
            p = 5;
            alg_flag = NE10_FFT_ALG_ANY;
        }
        else if ((n % 3) == 0)
        {
            p = 3;
            alg_flag = NE10_FFT_ALG_ANY;
        }
        else // stop factoring
        {
            p = n;
            alg_flag = NE10_FFT_ALG_ANY;
        }

        n /= p;
        facbuf[2 * i] = p;
        facbuf[2 * i + 1] = n;
        i++;
        stage_num++;
    }
    while (n > 1);
    facbuf[0] = stage_num;
    facbuf[1] = stride_max / p;

    if (stage_num > 21)
    {
        // Since nfft is ne10_int32_t, stage_num can never be greater than 21,
        // because 3^21 > 2^32
        return NE10_ERR;
    }

    facbuf[2 * i] = alg_flag;
    return NE10_OK;
}

// Twiddles matrix [radix-1][mstride]
// First column (k == 0) is ignored because phase == 1, and
// twiddle = (1.0, 0.0).
void ne10_fft_generate_twiddles_line_float32 (ne10_fft_cpx_float32_t * twiddles,
        const ne10_int32_t mstride,
        const ne10_int32_t fstride,
        const ne10_int32_t radix,
        const ne10_int32_t nfft)
{
    ne10_int32_t j, k;
    ne10_float32_t phase;
    const ne10_float64_t pi = NE10_PI;

    for (j = 0; j < mstride; j++)
    {
        for (k = 1; k < radix; k++) // phase = 1 when k = 0
        {
            phase = -2 * pi * fstride * k * j / nfft;
            twiddles[mstride * (k - 1) + j].r = (ne10_float32_t) cos (phase);
            twiddles[mstride * (k - 1) + j].i = (ne10_float32_t) sin (phase);
        } // radix
    } // mstride
}

// Transposed twiddles matrix [mstride][radix-1]
// First row (k == 0) is ignored because phase == 1, and
// twiddle = (1.0, 0.0).
// Transposed twiddle tables are used in RFFT to avoid memory access by a large
// stride.
void ne10_fft_generate_twiddles_line_transposed_float32 (
    ne10_fft_cpx_float32_t* twiddles,
    const ne10_int32_t mstride,
    const ne10_int32_t fstride,
    const ne10_int32_t radix,
    const ne10_int32_t nfft)
{
    ne10_int32_t j, k;
    ne10_float32_t phase;
    const ne10_float64_t pi = NE10_PI;

    for (j = 0; j < mstride; j++)
    {
        for (k = 1; k < radix; k++) // phase = 1 when k = 0
        {
            phase = -2 * pi * fstride * k * j / nfft;
            twiddles[(radix - 1) * j + k - 1].r = (ne10_float32_t) cos (phase);
            twiddles[(radix - 1) * j + k - 1].i = (ne10_float32_t) sin (phase);
        } // radix
    } // mstride
}
typedef void (*line_generator_float32)(ne10_fft_cpx_float32_t*,
      const ne10_int32_t,
      const ne10_int32_t,
      const ne10_int32_t,
      const ne10_int32_t);

ne10_fft_cpx_float32_t* ne10_fft_generate_twiddles_impl_float32 (
      line_generator_float32 generator,
      ne10_fft_cpx_float32_t * twiddles,
      const ne10_int32_t * factors,
      const ne10_int32_t nfft)
{
    ne10_int32_t stage_count = factors[0];
    ne10_int32_t fstride = factors[1];
    ne10_int32_t mstride;
    ne10_int32_t cur_radix; // current radix

    // for first stage
    cur_radix = factors[2 * stage_count];
    if (cur_radix % 2) // current radix is not 4 or 2
    {
        twiddles[0].r = 1.0;
        twiddles[0].i = 0.0;
        twiddles += 1;
        generator (twiddles, 1, fstride, cur_radix, nfft);
        twiddles += cur_radix - 1;
    }
    stage_count --;

    // for other stage
    for (; stage_count > 0; stage_count --)
    {
        cur_radix = factors[2 * stage_count];
        fstride /= cur_radix;
        mstride = factors[2 * stage_count + 1];
        generator (twiddles, mstride, fstride, cur_radix, nfft);
        twiddles += mstride * (cur_radix - 1);
    } // stage_count

    return twiddles;
}

ne10_fft_cpx_float32_t* ne10_fft_generate_twiddles_float32 (ne10_fft_cpx_float32_t * twiddles,
        const ne10_int32_t * factors,
        const ne10_int32_t nfft )
{
    line_generator_float32 generator = ne10_fft_generate_twiddles_line_float32;
    twiddles = ne10_fft_generate_twiddles_impl_float32(generator,
        twiddles, factors, nfft);
    return twiddles;
}

ne10_fft_cpx_float32_t* ne10_fft_generate_twiddles_transposed_float32 (
      ne10_fft_cpx_float32_t * twiddles,
      const ne10_int32_t * factors,
      const ne10_int32_t nfft)
{
    line_generator_float32 generator =
        ne10_fft_generate_twiddles_line_transposed_float32;
    twiddles = ne10_fft_generate_twiddles_impl_float32(generator,
        twiddles, factors, nfft);
    return twiddles;
}

void ne10_fft_destroy_r2c_float32 (ne10_fft_r2c_cfg_float32_t cfg)
{
    free(cfg);
}
