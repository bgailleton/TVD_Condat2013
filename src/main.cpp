// Core implementation of 1-D Total Variation Denoising using Condat's 2013
// and 2017 algorithms. The routines operate directly on raw pointers and are
// exposed to Python through pybind11 with NumPy array support.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <vector>
#include <limits>

namespace py = pybind11;

namespace tvd {

// -----------------------------------------------------------------------------
// 2013 algorithm (version 1)
// -----------------------------------------------------------------------------
/**
 * @brief Perform 1‑D total variation denoising.
 *
 * Direct implementation of Laurent Condat's 2013 algorithm.  The body is kept
 * as close as possible to the original C reference code; only minimal changes
 * were made to operate on generic floating point types and integrate with the
 * Python wrapper.
 *
 * @tparam T    Floating point type (float or double).
 * @param input Noisy input signal.
 * @param output Buffer receiving the denoised signal.  ``input`` and ``output``
 *               may point to the same memory for in‑place operation.
 * @param width Number of samples in the signal.
 * @param lambda Regularisation parameter (non‑negative).
 */
template <typename T>
void tv1d_denoise(const T *input, T *output, int width, T lambda)
{
    if (width > 0)
    {
        int k = 0, k0 = 0;                    // current index and segment start
        T umin = lambda, umax = -lambda;       // dual variables
        T vmin = input[0] - lambda, vmax = input[0] + lambda; // segment bounds
        int kplus = 0, kminus = 0;             // last positions of constraints
        const T twolambda = static_cast<T>(2) * lambda;
        const T minlambda = -lambda;
        for (;;)                                   // iterative process
        {
            while (k == width - 1)                 // apply boundary condition
            {
                if (umin < static_cast<T>(0))
                {
                    do
                        output[k0++] = vmin;
                    while (k0 <= kminus);
                    umax = (vmin = input[kminus = k = k0]) + (umin = lambda) - vmax;
                }
                else if (umax > static_cast<T>(0))
                {
                    do
                        output[k0++] = vmax;
                    while (k0 <= kplus);
                    umin = (vmax = input[kplus = k = k0]) + (umax = minlambda) - vmin;
                }
                else
                {
                    vmin += umin / (k - k0 + 1);
                    do
                        output[k0++] = vmin;
                    while (k0 <= k);
                    return;
                }
            }
            if ((umin += input[k + 1] - vmin) < minlambda) // negative jump
            {
                do
                    output[k0++] = vmin;
                while (k0 <= kminus);
                vmax = (vmin = input[kplus = kminus = k = k0]) + twolambda;
                umin = lambda;
                umax = minlambda;
            }
            else if ((umax += input[k + 1] - vmax) > lambda) // positive jump
            {
                do
                    output[k0++] = vmax;
                while (k0 <= kplus);
                vmin = (vmax = input[kplus = kminus = k = k0]) - twolambda;
                umin = lambda;
                umax = minlambda;
            }
            else                                  // no jump, continue
            {
                k++;
                if (umin >= lambda)
                {
                    vmin += (umin - lambda) / ((kminus = k) - k0 + 1);
                    umin = lambda;
                }
                if (umax <= minlambda)
                {
                    vmax += (umax + lambda) / ((kplus = k) - k0 + 1);
                    umax = minlambda;
                }
            }
        }
    }
}

/**
 * @brief Fused Lasso Signal Approximator.
 *
 * Extension of the above routine including an additional ``mu`` parameter.
 * The implementation mirrors the reference C code closely and supports
 * in‑place operation.
 */
template <typename T>
void fused_lasso(const T *input, T *output, int width, T lambda, T mu)
{
    if (width > 0)
    {
        int k = 0, k0 = 0;
        T umin = lambda, umax = -lambda;
        T vmin = input[0] - lambda, vmax = input[0] + lambda;
        int kplus = 0, kminus = 0;
        const T twolambda = static_cast<T>(2) * lambda;
        const T minlambda = -lambda;
        for (;;)
        {
            while (k == width - 1)
            {
                if (umin < static_cast<T>(0))
                {
                    vmin = vmin > mu ? vmin - mu
                                      : vmin < -mu ? vmin + mu : static_cast<T>(0);
                    do
                        output[k0++] = vmin;
                    while (k0 <= kminus);
                    umax = (vmin = input[kminus = k = k0]) + (umin = lambda) - vmax;
                }
                else if (umax > static_cast<T>(0))
                {
                    vmax = vmax > mu ? vmax - mu
                                      : vmax < -mu ? vmax + mu : static_cast<T>(0);
                    do
                        output[k0++] = vmax;
                    while (k0 <= kplus);
                    umin = (vmax = input[kplus = k = k0]) + (umax = minlambda) - vmin;
                }
                else
                {
                    vmin += umin / (k - k0 + 1);
                    vmin = vmin > mu ? vmin - mu
                                     : vmin < -mu ? vmin + mu : static_cast<T>(0);
                    do
                        output[k0++] = vmin;
                    while (k0 <= k);
                    return;
                }
            }
            if ((umin += input[k + 1] - vmin) < minlambda)
            {
                vmin = vmin > mu ? vmin - mu
                                  : vmin < -mu ? vmin + mu : static_cast<T>(0);
                do
                    output[k0++] = vmin;
                while (k0 <= kminus);
                vmax = (vmin = input[kplus = kminus = k = k0]) + twolambda;
                umin = lambda;
                umax = minlambda;
            }
            else if ((umax += input[k + 1] - vmax) > lambda)
            {
                vmax = vmax > mu ? vmax - mu
                                  : vmax < -mu ? vmax + mu : static_cast<T>(0);
                do
                    output[k0++] = vmax;
                while (k0 <= kplus);
                vmin = (vmax = input[kplus = kminus = k = k0]) - twolambda;
                umin = lambda;
                umax = minlambda;
            }
            else
            {
                k++;
                if (umin >= lambda)
                {
                    vmin += (umin - lambda) / ((kminus = k) - k0 + 1);
                    umin = lambda;
                }
                if (umax <= minlambda)
                {
                    vmax += (umax + lambda) / ((kplus = k) - k0 + 1);
                    umax = minlambda;
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// 2017 algorithm (version 2) by Laurent Condat under the CeCILL licence.
// -----------------------------------------------------------------------------

/**
 * @brief 1‑D total variation denoising (Condat 2017).
 *
 * Implementation of the revised algorithm published in 2017, which is a
 * simplified and slightly faster variant of the original 2013 routine.
 * The code is a C++ port of Laurent Condat's reference implementation and
 * carries the same CeCILL licence (GPL compatible).
 *
 * @tparam T        Floating point type (float or double).
 * @param input     Pointer to the noisy input signal.
 * @param output    Pointer to a pre‑allocated buffer that will receive the
 *                  denoised signal.
 * @param width     Number of samples in the signal.
 * @param lambda    Regularisation parameter controlling the amount of
 *                  smoothing.
 */
// The reference implementation operates in double precision.  We keep an
// internal double routine and expose a templated wrapper that casts inputs and
// outputs appropriately so that float arrays still benefit from the increased
// numerical robustness.
static void tv1d_denoise_v2_double(const double *input, double *output,
                                   unsigned int width, double lambda)
{
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);

    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;
    double output_low_first = input[0] - lambda;
    double output_low_curr = output_low_first;
    double output_up_first = input[0] + lambda;
    double output_up_curr = output_up_first;
    const double twolambda = 2.0 * lambda;

    if (width == 1)
    {
        output[0] = input[0];
        return;
    }

    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i < width; i++)
    {
        if (input[i] >= output_low_curr)
        {
            if (input[i] <= output_up_curr)
            {
                output_up_curr +=
                    (input[i] - output_up_curr) / (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                    output_up_curr += (output[ind] - output_up_curr) *
                                      (double)(indstart_up[j_up--] - ind) /
                                      (i - ind + 1);
                if (j_up == jseg)
                {
                    while ((output_up_curr <= output_low_first) && (jseg < j_low))
                    {
                        indjseg2 = indstart_low[++jseg];
                        output_up_curr += (output_up_curr - output_low_first) *
                                          (double)(indjseg2 - indjseg) /
                                          (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_low_first;
                        output_low_first = output[indjseg];
                    }
                    output_up_first = output_up_curr;
                    indstart_up[j_up = jseg] = indjseg;
                }
                else
                    output[indstart_up[j_up]] = output_up_curr;
            }
            else
            {
                output_up_curr = output[i] = input[indstart_up[++j_up] = i];
                output_low_curr +=
                    (input[i] - output_low_curr) /
                    (i - indstart_low[j_low] + 1);
                output[indjseg] = output_low_first;
                while ((j_low > jseg) &&
                       (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
                    output_low_curr += (output[ind] - output_low_curr) *
                                       (double)(indstart_low[j_low--] - ind) /
                                       (i - ind + 1);
                if (j_low == jseg)
                {
                    while ((output_low_curr >= output_up_first) && (jseg < j_up))
                    {
                        indjseg2 = indstart_up[++jseg];
                        output_low_curr += (output_low_curr - output_up_first) *
                                          (double)(indjseg2 - indjseg) /
                                          (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_up_first;
                        output_up_first = output[indjseg];
                    }
                    if ((indstart_low[j_low = jseg] = indjseg) == i)
                        output_low_first = output_up_first - twolambda;
                    else
                        output_low_first = output_low_curr;
                }
                else
                    output[indstart_low[j_low]] = output_low_curr;
            }
        }
        else
        {
            output_up_curr += ((output_low_curr = output[i] =
                                    input[indstart_low[++j_low] = i]) -
                                output_up_curr) /
                               (i - indstart_up[j_up] + 1);
            output[indjseg] = output_up_first;
            while ((j_up > jseg) &&
                   (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                output_up_curr += (output[ind] - output_up_curr) *
                                  (double)(indstart_up[j_up--] - ind) /
                                  (i - ind + 1);
            if (j_up == jseg)
            {
                while ((output_up_curr <= output_low_first) && (jseg < j_low))
                {
                    indjseg2 = indstart_low[++jseg];
                    output_up_curr += (output_up_curr - output_low_first) *
                                      (double)(indjseg2 - indjseg) /
                                      (i - indjseg2 + 1);
                    while (indjseg < indjseg2)
                        output[indjseg++] = output_low_first;
                    output_low_first = output[indjseg];
                }
                if ((indstart_up[j_up = jseg] = indjseg) == i)
                    output_up_first = output_low_first + twolambda;
                else
                    output_up_first = output_up_curr;
            }
            else
                output[indstart_up[j_up]] = output_up_curr;
        }
    }

    if (input[i] + lambda <= output_low_curr)
    {
        while (jseg < j_low)
        {
            indjseg2 = indstart_low[++jseg];
            while (indjseg < indjseg2)
                output[indjseg++] = output_low_first;
            output_low_first = output[indjseg];
        }
        while (indjseg < i)
            output[indjseg++] = output_low_first;
        output[indjseg] = input[i] + lambda;
    }
    else if (input[i] - lambda >= output_up_curr)
    {
        while (jseg < j_up)
        {
            indjseg2 = indstart_up[++jseg];
            while (indjseg < indjseg2)
                output[indjseg++] = output_up_first;
            output_up_first = output[indjseg];
        }
        while (indjseg < i)
            output[indjseg++] = output_up_first;
        output[indjseg] = input[i] - lambda;
    }
    else
    {
        output_low_curr += (input[i] + lambda - output_low_curr) /
                           (i - indstart_low[j_low] + 1);
        output[indjseg] = output_low_first;
        while ((j_low > jseg) &&
               (output_low_curr >= output[ind = indstart_low[j_low - 1]]))
            output_low_curr += (output[ind] - output_low_curr) *
                               (double)(indstart_low[j_low--] - ind) /
                               (i - ind + 1);
        if (j_low == jseg)
        {
            if (output_up_first >= output_low_curr)
                while (indjseg <= i)
                    output[indjseg++] = output_low_curr;
            else
            {
                output_up_curr += (input[i] - lambda - output_up_curr) /
                                   (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                    output_up_curr += (output[ind] - output_up_curr) *
                                       (double)(indstart_up[j_up--] - ind) /
                                       (i - ind + 1);
                while (jseg < j_up)
                {
                    indjseg2 = indstart_up[++jseg];
                    while (indjseg < indjseg2)
                        output[indjseg++] = output_up_first;
                    output_up_first = output[indjseg];
                }
                indjseg = indstart_up[j_up];
                while (indjseg <= i)
                    output[indjseg++] = output_up_curr;
            }
        }
        else
        {
            while (jseg < j_low)
            {
                indjseg2 = indstart_low[++jseg];
                while (indjseg < indjseg2)
                    output[indjseg++] = output_low_first;
                output_low_first = output[indjseg];
            }
            indjseg = indstart_low[j_low];
            while (indjseg <= i)
                output[indjseg++] = output_low_curr;
        }
    }
}

template <typename T>
void tv1d_denoise_v2(const T *input, T *output, unsigned int width, T lambda)
{
    // Run the core algorithm in double precision for robustness and cast back
    // to the original dtype once finished.
    std::vector<double> in(width), out(width);
    for (unsigned int k = 0; k < width; ++k)
        in[k] = static_cast<double>(input[k]);
    tv1d_denoise_v2_double(in.data(), out.data(), width,
                           static_cast<double>(lambda));
    for (unsigned int k = 0; k < width; ++k)
        output[k] = static_cast<T>(out[k]);
}

/**
 * @brief 1‑D total variation denoising using the taut string algorithm.
 *
 * Adapted from the Matlab implementation by Lutz Dümbgen.
 * Uses double precision internally for numerical stability.
 */
template <typename T>
void tv1d_denoise_tautstring(const T *input, T *output, int width, T lambda)
{
    if (width <= 0)
        return;
    using FT = double;
    int N = width + 1;
    std::vector<int> index_low(N), index_up(N), index(N);
    std::vector<FT> slope_low(N), slope_up(N), z(N), y_low(N), y_up(N);

    int s_low = 0, c_low = 0, s_up = 0, c_up = 0, c = 0;
    int i = 2;

    y_low[0] = y_up[0] = 0;
    y_low[1] = static_cast<FT>(input[0]) - static_cast<FT>(lambda);
    y_up[1] = static_cast<FT>(input[0]) + static_cast<FT>(lambda);
    for (; i < N; ++i)
    {
        y_low[i] = y_low[i - 1] + static_cast<FT>(input[i - 1]);
        y_up[i] = y_up[i - 1] + static_cast<FT>(input[i - 1]);
    }
    y_low[N - 1] += static_cast<FT>(lambda);
    y_up[N - 1] -= static_cast<FT>(lambda);

    slope_low[0] = std::numeric_limits<FT>::infinity();
    slope_up[0] = -std::numeric_limits<FT>::infinity();
    z[0] = y_low[0];

    for (i = 1; i < N; ++i)
    {
        index_low[++c_low] = index_up[++c_up] = i;
        slope_low[c_low] = y_low[i] - y_low[i - 1];
        while ((c_low > s_low + 1) &&
               (slope_low[std::max(s_low, c_low - 1)] <= slope_low[c_low]))
        {
            index_low[--c_low] = i;
            if (c_low > s_low + 1)
                slope_low[c_low] = (y_low[i] - y_low[index_low[c_low - 1]]) /
                                    (i - index_low[c_low - 1]);
            else
                slope_low[c_low] = (y_low[i] - z[c]) / (i - index[c]);
        }

        slope_up[c_up] = y_up[i] - y_up[i - 1];
        while ((c_up > s_up + 1) &&
               (slope_up[std::max(c_up - 1, s_up)] >= slope_up[c_up]))
        {
            index_up[--c_up] = i;
            if (c_up > s_up + 1)
                slope_up[c_up] = (y_up[i] - y_up[index_up[c_up - 1]]) /
                                  (i - index_up[c_up - 1]);
            else
                slope_up[c_up] = (y_up[i] - z[c]) / (i - index[c]);
        }

        while ((c_low == s_low + 1) && (c_up > s_up + 1) &&
               (slope_low[c_low] >= slope_up[s_up + 1]))
        {
            index[++c] = index_up[++s_up];
            z[c] = y_up[index[c]];
            index_low[s_low] = index[c];
            slope_low[c_low] = (y_low[i] - z[c]) / (i - index[c]);
        }
        while ((c_up == s_up + 1) && (c_low > s_low + 1) &&
               (slope_up[c_up] <= slope_low[s_low + 1]))
        {
            index[++c] = index_low[++s_low];
            z[c] = y_low[index[c]];
            index_up[s_up] = index[c];
            slope_up[c_up] = (y_up[i] - z[c]) / (i - index[c]);
        }
    }

    for (i = 1; i <= c_low - s_low; ++i)
        z[c + i] = y_low[index[c + i] = index_low[s_low + i]];
    c += c_low - s_low;

    int j = 0;
    T a;
    for (i = 1; i <= c; ++i)
    {
        a = static_cast<T>((z[i] - z[i - 1]) / (index[i] - index[i - 1]));
        while (j < index[i])
            output[j++] = a;
    }
}

// -----------------------------------------------------------------------------
// Wrappers operating on NumPy arrays
// -----------------------------------------------------------------------------

/**
 * @brief Convenience wrapper around tv1d_denoise for NumPy arrays.
 *
 * The input array is copied into a contiguous buffer (if required), passed to
 * the core implementation and the denoised result is returned as a new
 * NumPy array.  The function preserves the dtype of the input array and is
 * exposed to Python as ``tvd_2013``.
 *
 * @tparam T    Floating point type (float or double).
 * @param in    1‑D NumPy array containing the noisy signal.
 * @param lambda Regularisation parameter controlling the amount of smoothing.
 * @returns     New NumPy array with the denoised signal.
 */
template <typename T>
py::array_t<T> tvd_2013(py::array_t<T, py::array::c_style | py::array::forcecast> in,
                        double lambda)
{
    auto buf = in.request();
    auto n = static_cast<int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);
    py::array_t<T> result(buf.size);
    tv1d_denoise(data, result.mutable_data(), n, static_cast<T>(lambda));
    return result;
}

/**
 * @brief Wrapper for the revised 2017 denoising algorithm.
 *
 * This function behaves identically to :func:`tvd_2013` but calls the ``v2``
 * implementation which trades a small amount of precision for improved
 * performance.  Both single and double precision signals are supported.
 */
template <typename T>
py::array_t<T> tvd_2017(py::array_t<T, py::array::c_style | py::array::forcecast> in,
                        double lambda)
{
    auto buf = in.request();
    auto n = static_cast<unsigned int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);
    py::array_t<T> result(buf.size);
    tv1d_denoise_v2(data, result.mutable_data(), n, static_cast<T>(lambda));
    return result;
}

/**
 * @brief Wrapper for the taut string denoising algorithm.
 */
template <typename T>
py::array_t<T> tvd_tautstring(
    py::array_t<T, py::array::c_style | py::array::forcecast> in, double lambda)
{
    auto buf = in.request();
    auto n = static_cast<int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);
    py::array_t<T> result(buf.size);
    tv1d_denoise_tautstring(data, result.mutable_data(), n,
                            static_cast<T>(lambda));
    return result;
}

/**
 * @brief Wrapper for the fused lasso signal approximator.
 */
template <typename T>
py::array_t<T> fused_lasso_py(
    py::array_t<T, py::array::c_style | py::array::forcecast> in, double lambda,
    double mu)
{
    auto buf = in.request();
    auto n = static_cast<int>(buf.size);
    const T *data = static_cast<T *>(buf.ptr);
    py::array_t<T> result(buf.size);
    fused_lasso(data, result.mutable_data(), n, static_cast<T>(lambda),
                static_cast<T>(mu));
    return result;
}

} // namespace tvd
PYBIND11_MODULE(TVDCondat2013, m)
{
    // Module level documentation with references to the original publications
    m.doc() = R"pbdoc(
        Python bindings for 1-D total variation denoising based on
        Laurent Condat's algorithms\ [Condat2013]_\ [Condat2017]_.
        Includes a taut string variant adapted from Lutz Dümbgen's Matlab code.

        .. [Condat2013] L. Condat, "A Direct Algorithm for 1D Total Variation
           Denoising," *IEEE Signal Processing Letters*, 2013.
        .. [Condat2017] L. Condat, "Fast Projection onto the Simplex and the
           L1 Ball," *Mathematical Programming*, 2017.

        Copyright (c) 2013-2025 Laurent Condat and contributors.
    )pbdoc";

    const char *tvd_doc = R"pbdoc(
        Apply 1-D total variation denoising to a signal.

        :param numpy.ndarray signal: 1-D array of ``float32`` or ``float64``
            values.
        :param float lambda: Regularisation parameter controlling smoothing.
        :returns: Denoised signal with the same dtype as ``signal``.
        :rtype: numpy.ndarray

        Copyright (c) 2013-2025 Laurent Condat.
    )pbdoc";

    m.def("tvd_2013", &tvd::tvd_2013<double>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);
    m.def("tvd_2013", &tvd::tvd_2013<float>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);

    m.def("tvd_2017", &tvd::tvd_2017<double>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);
    m.def("tvd_2017", &tvd::tvd_2017<float>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);

    m.def("tvd_tautstring", &tvd::tvd_tautstring<double>,
          py::arg("signal").noconvert(), py::arg("lambda"), tvd_doc);
    m.def("tvd_tautstring", &tvd::tvd_tautstring<float>,
          py::arg("signal").noconvert(), py::arg("lambda"), tvd_doc);

    const char *fl_doc = R"pbdoc(
        Fused lasso signal approximation.

        :param numpy.ndarray signal: 1-D array of ``float32`` or ``float64``
            values.
        :param float lambda: Total variation regularisation parameter.
        :param float mu: L1 penalty on the signal values.
        :returns: Denoised signal with the same dtype as ``signal``.
        :rtype: numpy.ndarray

        Copyright (c) 2013-2025 Laurent Condat.
    )pbdoc";

    m.def("fused_lasso", &tvd::fused_lasso_py<double>,
          py::arg("signal").noconvert(), py::arg("lambda"), py::arg("mu"),
          fl_doc);
    m.def("fused_lasso", &tvd::fused_lasso_py<float>,
          py::arg("signal").noconvert(), py::arg("lambda"), py::arg("mu"),
          fl_doc);
}

