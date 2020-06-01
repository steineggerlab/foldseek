/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Frank Astier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * Utilities for working with quaternion.
 */
#ifndef QUATERNIONS_UTILS_H
#define QUATERNIONS_UTILS_H

#include <cmath>
#include <iterator>
#include <type_traits>

namespace quaternion {

/**
 * A few macros to use in template declarations in order to select
 * the correct function. These are annoying to read directly in the source code,
 * so I #define them here, and the source becomes more readable.
 * NOTE: if not keeping these macros on one line, use \ at end of line, and leave last line without "\"
 * or trouble will ensue.
 */
#define IS_ITERATOR(X) typename std::enable_if<!std::is_same<typename std::iterator_traits<X>::value_type, void>::value>::type* =nullptr
#define IS_CONVERTIBLE(FROM, TO) typename std::enable_if<std::is_convertible<FROM,TO>::value>::type* =nullptr

/**
 * Utility function to work with numbers approximately equal to zero.
 * If eps == 0, does a "hard" comparison to 0.
 * Otherwise, uses a ball of radius eps around 0. If the scalar is inside
 * that ball, it is equivalent to 0.
 */
    template<typename T1, typename T2>
    inline bool is_scalar_zero(T1 x, T2 eps = 0) {
        typedef typename std::common_type<T1,T2>::type T;
        T xx = static_cast<T>(x);
        T ee = static_cast<T>(eps);
        return std::abs(xx) <= ee;
    }

/**
 * Compares 2 floating point numbers "relatively": if the numbers are
 * very large, differences are still "small" if they are "small"
 * relative to the magnitudes of the quantities.
 * TODO: need absolute difference comparison (not relative) too
 */
    template<typename T1, typename T2, typename T3>
    inline bool is_nearly_equal(T1 x, T2 y, T3 eps = 0) {
        typedef typename std::common_type<T1,T2,T3>::type T;
        T xx = static_cast<T>(x);
        T yy = static_cast<T>(y);
        T ee = static_cast<T>(eps);
        if (xx == 0)
            return is_scalar_zero(yy, ee);
        else if (yy == 0)
            return is_scalar_zero(xx, ee);
        else
            return is_scalar_zero((xx - yy) / std::min(xx, yy), ee);
    }

} // end namespace quaternion
#endif //QUATERNIONS_UTILS_H
