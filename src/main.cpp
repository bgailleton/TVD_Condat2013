// Core implementation of 1-D Total Variation Denoising using Condat's 2013
// and 2017 algorithms. The routines operate directly on raw pointers and are
// exposed to Python through pybind11 with NumPy array support.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <vector>

namespace py = pybind11;

namespace tvd {

// -----------------------------------------------------------------------------
// 2013 algorithm (version 1)
// -----------------------------------------------------------------------------
/**
 * @brief Perform 1‑D total variation denoising.
 *
 * This is a direct port of Condat's 2013 algorithm for solving the
 * 1‑D total‑variation denoising problem.  It operates on raw pointers for
 * maximal performance and expects the caller to provide the input and
 * output buffers.  The implementation follows the notation of the
 * original paper and exposes the routine through pybind11.
 *
 * @tparam T        Floating point type (float or double).
 * @param input     Pointer to the noisy input signal.
 * @param output    Pointer to a pre‑allocated buffer where the denoised
 *                  signal will be written.
 * @param width     Number of samples in the signal.
 * @param lambda    Regularisation parameter controlling the amount of
 *                  smoothing.  Larger values yield flatter signals.
 */
template <typename T>
void tv1d_denoise(const T *input, T *output, unsigned int width, T lambda)
{
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);

    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;
    T output_low_first = input[0] - lambda;
    T output_low_curr = output_low_first;
    T output_up_first = input[0] + lambda;
    T output_up_curr = output_up_first;
    T twolambda = static_cast<T>(2) * lambda;

    if (width == 1)
    {
        output[0] = input[0];
        return;
    }

    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i < width; ++i)
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
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
                        (i - ind + 1);
                if (j_up == jseg)
                {
                    while ((output_up_curr <= output_low_first) &&
                           (jseg < j_low))
                    {
                        indjseg2 = indstart_low[++jseg];
                        output_up_curr +=
                            (output_up_curr - output_low_first) *
                            static_cast<T>(indjseg2 - indjseg) /
                            (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_low_first;
                        output_low_first = output[indjseg];
                    }
                    output_up_first = output_up_curr;
                    indstart_up[j_up = jseg] = indjseg;
                }
                else
                {
                    output[indstart_up[j_up]] = output_up_curr;
                }
            }
            else
            {
                output_up_curr = output[i] = input[indstart_up[++j_up] = i];
                output_low_curr +=
                    (input[i] - output_low_curr) /
                    (i - indstart_low[j_low] + 1);
                output[indjseg] = output_low_first;
                while ((j_low > jseg) &&
                       (output_low_curr >=
                        output[ind = indstart_low[j_low - 1]]))
                    output_low_curr +=
                        (output[ind] - output_low_curr) *
                        static_cast<T>(indstart_low[j_low--] - ind) /
                        (i - ind + 1);
                if (j_low == jseg)
                {
                    while ((output_low_curr >= output_up_first) &&
                           (jseg < j_up))
                    {
                        indjseg2 = indstart_up[++jseg];
                        output_low_curr +=
                            (output_low_curr - output_up_first) *
                            static_cast<T>(indjseg2 - indjseg) /
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
                {
                    output[indstart_low[j_low]] = output_low_curr;
                }
            }
        }
        else
        {
            output_up_curr +=
                ((output_low_curr = output[i] =
                      input[indstart_low[++j_low] = i]) -
                 output_up_curr) /
                (i - indstart_up[j_up] + 1);
            output[indjseg] = output_up_first;
            while ((j_up > jseg) &&
                   (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                output_up_curr +=
                    (output[ind] - output_up_curr) *
                    static_cast<T>(indstart_up[j_up--] - ind) /
                    (i - ind + 1);
            if (j_up == jseg)
            {
                while ((output_up_curr <= output_low_first) &&
                       (jseg < j_low))
                {
                    indjseg2 = indstart_low[++jseg];
                    output_up_curr +=
                        (output_up_curr - output_low_first) *
                        static_cast<T>(indjseg2 - indjseg) /
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
            {
                output[indstart_up[j_up]] = output_up_curr;
            }
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
            output_low_curr +=
                (output[ind] - output_low_curr) *
                static_cast<T>(indstart_low[j_low--] - ind) /
                (i - ind + 1);
        if (j_low == jseg)
        {
            if (output_up_first >= output_low_curr)
            {
                while (indjseg <= i)
                    output[indjseg++] = output_low_curr;
            }
            else
            {
                output_up_curr += (input[i] - lambda - output_up_curr) /
                                   (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <=
                        output[ind = indstart_up[j_up - 1]]))
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
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
template <typename T>
void tv1d_denoise_v2(const T *input, T *output, unsigned int width, T lambda)
{
    std::vector<unsigned int> indstart_low(width);
    std::vector<unsigned int> indstart_up(width);

    unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i = 1, indjseg2, ind;
    T output_low_first = input[0] - lambda;
    T output_low_curr = output_low_first;
    T output_up_first = input[0] + lambda;
    T output_up_curr = output_up_first;
    T twolambda = static_cast<T>(2) * lambda;

    if (width == 1)
    {
        output[0] = input[0];
        return;
    }

    indstart_low[0] = 0;
    indstart_up[0] = 0;
    width--;
    for (; i < width; ++i)
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
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
                        (i - ind + 1);
                if (j_up == jseg)
                {
                    while ((output_up_curr <= output_low_first) &&
                           (jseg < j_low))
                    {
                        indjseg2 = indstart_low[++jseg];
                        output_up_curr +=
                            (output_up_curr - output_low_first) *
                            static_cast<T>(indjseg2 - indjseg) /
                            (i - indjseg2 + 1);
                        while (indjseg < indjseg2)
                            output[indjseg++] = output_low_first;
                        output_low_first = output[indjseg];
                    }
                    output_up_first = output_up_curr;
                    indstart_up[j_up = jseg] = indjseg;
                }
                else
                {
                    output[indstart_up[j_up]] = output_up_curr;
                }
            }
            else
            {
                output_up_curr = output[i] = input[indstart_up[++j_up] = i];
                output_low_curr +=
                    (input[i] - output_low_curr) /
                    (i - indstart_low[j_low] + 1);
                output[indjseg] = output_low_first;
                while ((j_low > jseg) &&
                       (output_low_curr >=
                        output[ind = indstart_low[j_low - 1]]))
                    output_low_curr +=
                        (output[ind] - output_low_curr) *
                        static_cast<T>(indstart_low[j_low--] - ind) /
                        (i - ind + 1);
                if (j_low == jseg)
                {
                    while ((output_low_curr >= output_up_first) &&
                           (jseg < j_up))
                    {
                        indjseg2 = indstart_up[++jseg];
                        output_low_curr +=
                            (output_low_curr - output_up_first) *
                            static_cast<T>(indjseg2 - indjseg) /
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
                {
                    output[indstart_low[j_low]] = output_low_curr;
                }
            }
        }
        else
        {
            output_up_curr +=
                ((output_low_curr = output[i] =
                      input[indstart_low[++j_low] = i]) -
                 output_up_curr) /
                (i - indstart_up[j_up] + 1);
            output[indjseg] = output_up_first;
            while ((j_up > jseg) &&
                   (output_up_curr <= output[ind = indstart_up[j_up - 1]]))
                output_up_curr +=
                    (output[ind] - output_up_curr) *
                    static_cast<T>(indstart_up[j_up--] - ind) /
                    (i - ind + 1);
            if (j_up == jseg)
            {
                while ((output_up_curr <= output_low_first) &&
                       (jseg < j_low))
                {
                    indjseg2 = indstart_low[++jseg];
                    output_up_curr +=
                        (output_up_curr - output_low_first) *
                        static_cast<T>(indjseg2 - indjseg) /
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
            {
                output[indstart_up[j_up]] = output_up_curr;
            }
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
            output_low_curr +=
                (output[ind] - output_low_curr) *
                static_cast<T>(indstart_low[j_low--] - ind) /
                (i - ind + 1);
        if (j_low == jseg)
        {
            if (output_up_first >= output_low_curr)
            {
                while (indjseg <= i)
                    output[indjseg++] = output_low_curr;
            }
            else
            {
                output_up_curr += (input[i] - lambda - output_up_curr) /
                                   (i - indstart_up[j_up] + 1);
                output[indjseg] = output_up_first;
                while ((j_up > jseg) &&
                       (output_up_curr <=
                        output[ind = indstart_up[j_up - 1]]))
                    output_up_curr +=
                        (output[ind] - output_up_curr) *
                        static_cast<T>(indstart_up[j_up--] - ind) /
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
    auto n = static_cast<unsigned int>(buf.size);
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

} // namespace tvd
PYBIND11_MODULE(TVDCondat2013, m)
{
    // Module level documentation with references to the original publications
    m.doc() = R"pbdoc(
        Python bindings for 1-D total variation denoising based on
        Laurent Condat's algorithms\ [Condat2013]_\ [Condat2017]_.

        .. [Condat2013] L. Condat, "A Direct Algorithm for 1D Total Variation
           Denoising," *IEEE Signal Processing Letters*, 2013.
        .. [Condat2017] L. Condat, "Fast Projection onto the Simplex and the
           L1 Ball," *Mathematical Programming*, 2017.
    )pbdoc";

    const char *tvd_doc = R"pbdoc(
        Apply 1-D total variation denoising to a signal.

        :param numpy.ndarray signal: 1-D array of ``float32`` or ``float64``
            values.
        :param float lambda: Regularisation parameter controlling smoothing.
        :returns: Denoised signal with the same dtype as ``signal``.
        :rtype: numpy.ndarray
    )pbdoc";

    m.def("tvd_2013", &tvd::tvd_2013<double>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);
    m.def("tvd_2013", &tvd::tvd_2013<float>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);

    m.def("tvd_2017", &tvd::tvd_2017<double>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);
    m.def("tvd_2017", &tvd::tvd_2017<float>, py::arg("signal").noconvert(),
          py::arg("lambda"), tvd_doc);
}

