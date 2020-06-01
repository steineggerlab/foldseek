//
// Created by Martin Steinegger on 29.09.18.
//

#ifndef STRUCCLUST_QUATERNION_H
#define STRUCCLUST_QUATERNION_H

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
 * A Quaternion class.
 */

#ifndef QUATERNIONS_QUATERNION_H
#define QUATERNIONS_QUATERNION_H

#include <math.h> // for atan2, which handles signs for us properly

#include <limits>
#include <type_traits>
#include <array>
#include <complex>
#include <iterator>
#include <assert.h>

#include "quaternion_utils.h"

namespace quaternion {

/**
 * A Quaternion class.
 * TODO: quaternion logo
 * TODO: IEEE 754/IEC 559?
 * TODO: remove IS_CONVERTIBLE in favor std::common_type
 * TODO: arc-trans, rounding, rotate vector, slerp, fitting, maybe separate Euler
 */
    template<typename T =double>
    class Quaternion {
    public:
        /**
         * The type of each component of the 4 components of a Quaternion.
         * See below for allowed types.
         */
        typedef T value_type;

        /**
         * For now, the types that we can handle are restricted, to make sure
         * that we are not asked to handle a type that would lead to bugs
         * (e.g. unsigned T when we have signed operations, or complex T).
         * boost::rational is technically OK, and I verified in tests,
         * but I don't want to depend on boost in this header, so we'll allow
         * boost::rational if there is demand for it.
         *
         * WARNING: for bool and the integers, some operations won't compile right now,
         * in particular the transcendental functions.
         */
        static_assert(std::is_same<T, bool>()
                      || std::is_same<T, int>()
                      || std::is_same<T, long>()
                      || std::is_same<T, long long>()
                      || std::is_same<T, float>()
                      || std::is_same<T, double>()
                      || std::is_same<T, long double>(),
                      "Invalid scalar type for Quaternion");

        /**
        * Construct a Quaternion from at most 4 components of type T.
        * Specifying only a != 0 makes the Quaternion a real.
        * Specifying only a != and b != 0 makes the Quaternion an ordinary complex number.
        */
        Quaternion(T a = 0, T b = 0, T c = 0, T d = 0)
                : _a(a), _b(b), _c(c), _d(d) { }

        /**
         * Construct a Quaternion from at most 4 components of type T1.
         * Specifying only a != 0 makes the Quaternion a real.
         * Specifying only a != and b != 0 makes the Quaternion an ordinary complex number.
         * NOTE: IS_CONVERTIBLE to avoid ambiguity with constructor from iterator.
         */
        template<typename T1, IS_CONVERTIBLE(T1,T)>
        Quaternion(T1 a = 0, T1 b = 0, T1 c = 0, T1 d = 0)
                : _a(a), _b(b), _c(c), _d(d) { }

        /**
         * Construct a Quaternion from 1 or 2 std::complex<T>.
         */
        template<typename T1>
        Quaternion(const std::complex<T1>& x, const std::complex<T1>& y = std::complex<T1>(0, 0))
                : _a(x.real()), _b(x.imag()), _c(y.real()), _d(y.imag()) { }

        /**
         * Construct from an iterator to a range of 4 elements.
         * The iterator is dereferenced exactly 4 times, and advanced (++)
         * exactly 3 times.
         */
        template<typename It, IS_ITERATOR(It)>
        Quaternion(It it)
                : _a(*it), _b(*++it), _c(*++it), _d(*++it) { }

        /**
         * Copy constructor, from a Quaternion with another value type.
         */
        template<typename T1>
        Quaternion(const Quaternion<T1>& y)
                : _a(y.a()), _b(y.b()), _c(y.c()), _d(y.d()) { }

        /**
         * Assignment operator, from a Quaternion with another value type.
         */
        template<typename T1>
        Quaternion& operator=(const Quaternion<T1>& other) {
            _a = other.a();
            _b = other.b();
            _c = other.c();
            _d = other.d();
            return *this;
        }

        /**
         * Accessors for all 4 components of the Quaternion.
         */
        T a() const { return _a; }
        T b() const { return _b; }
        T c() const { return _c; }
        T d() const { return _d; }

        /**
         * The complex components of this Quaternion.
         */
        std::complex<T> c1() const { return {_a, _b}; }
        std::complex<T> c2() const { return {_c, _d}; }

        /**
         * Ordered list form.
         */
        std::array<T, 4> to_array() const {
            return {{_a, _b, _c, _d}};
        }

        /**
         * The real part of this Quaternion.
         */
        T real() const { return _a; }

        /**
         * The "unreal" part of this Quaternion, which is a Quaternion itself.
         */
        Quaternion unreal() const { return {0, _b, _c, _d}; }

        /**
         * The square of the norm of the Quaternion.
         * (The square is sometimes useful, and it avoids paying for a sqrt).
         */
        T norm_squared() const {
            return _a * _a + _b * _b + _c * _c + _d * _d;
        }

        /**
         * The norm of the Quaternion (the l2 norm).
         */
        T abs() const {
            return std::sqrt(norm_squared());
        }

        /**
         * The L2 norm of the "unreal" components of the Quaternion,
         * comes back often in computations.
         */
        T unreal_norm_squared() const {
            return _b * _b + _c * _c + _d * _d;
        }

        /**
         * Return true if this Quaternion is zero, false otherwise.
         */
        template <typename T1 =T>
        bool is_zero(T1 eps = 0) const {
            return is_scalar_zero(_a, eps)
                   && is_scalar_zero(_b, eps)
                   && is_scalar_zero(_c, eps)
                   && is_scalar_zero(_d, eps);
        }

        /**
         * Return true if this Quaternion is not zero, false otherwise.
         */
        template <typename T1 =T>
        bool is_non_zero(T1 eps = 0) const {
            return !is_zero(eps);
        }

        /**
         * Return true if any component of this quaternion is nan, false otherwise.
         * TODO: use "isnan" instead of "is_nan" to keep consistent?
         */
        bool is_nan() const {
            return std::isnan(_a) || std::isnan(_b) || std::isnan(_c) || std::isnan(_d);
        }

        /**
        * Return true if any component of this quaternion is inf, false otherwise.
        */
        bool is_inf() const {
            return std::isinf(_a) || std::isinf(_b) || std::isinf(_c) || std::isinf(_d);
        }

        /**
        * Return true if all the components of this quaternion are finite, false otherwise.
        */
        bool is_finite() const {
            return std::isfinite(_a) && std::isfinite(_b) && std::isfinite(_c) && std::isfinite(_d);
        }

        /**
         * Return true if this Quaternion has norm 1, false otherwise.
         */
        template<typename T1 =T>
        bool is_unit(T1 eps = 0) const {
            return is_scalar_zero(norm_squared() - T(1), eps);
        }

        /**
         * Return true if this Quaternion is real, false otherwise.
         */
        template<typename T1 =T>
        bool is_real(T1 eps = 0) const {
            return is_scalar_zero(_b, eps)
                   && is_scalar_zero(_c, eps)
                   && is_scalar_zero(_d, eps);
        }

        /**
         * Return true if this Quaternion is complex, false otherwise.
         */
        template<typename T1 =T>
        bool is_complex(T1 eps = 0) const {
            return is_scalar_zero(_c, eps) && is_scalar_zero(_d, eps);
        }

        /**
         * Return true if this Quaternion is real, false otherwise.
         */
        template<typename T1 =T>
        bool is_unreal(T1 eps = 0) const {
            return is_scalar_zero(_a, eps)
                   && !(is_scalar_zero(_b, eps)
                        && is_scalar_zero(_c, eps)
                        && is_scalar_zero(_d, eps));
        }

        /**
         * Unary plus.
         */
        Quaternion operator+() const {
            return *this;
        }

        /**
         * Unary minus.
         */
        Quaternion operator-() const {
            return {-_a, -_b, -_c, -_d};
        }

        /**
         * Unary +=.
         */
        Quaternion operator+=(T y) {
            _a += y;
            return *this;
        }

        /**
         * Unary +=.
         */
        Quaternion operator-=(T y) {
            _a -= y;
            return *this;
        }

        /**
         * Scaling by a constant.
         */
        Quaternion operator*=(T k) {
            _a = k * _a;
            _b = k * _b;
            _c = k * _c;
            _d = k * _d;
            return *this;
        }

        /**
         * Dividing by a constant.
         */
        Quaternion operator/=(T k) {

            _a /= k;
            _b /= k;
            _c /= k;
            _d /= k;
            return *this;
        }

        /**
         * Unary +=.
         */
        template<typename T1>
        Quaternion operator+=(const std::complex<T1>& y) {
            _a += y.real();
            _b += y.imag();
            return *this;
        }

        /**
         * Unary -=.
         */
        template<typename T1>
        Quaternion operator-=(const std::complex<T1>& y) {
            _a -= y.real();
            _b -= y.imag();
            return *this;
        }

        /**
         * Unary *=.
         */
        template<typename T1>
        Quaternion operator*=(const std::complex<T1>& y) {

            T at = _a * y.real() - _b * y.imag();
            T bt = _a * y.imag() + _b * y.real();
            T ct = _c * y.real() + _d * y.imag();
            T dt = -_c * y.imag() + _d * y.real();

            _a = at;
            _b = bt;
            _c = ct;
            _d = dt;

            return *this;
        }

        /**
         * Unary /=.
         */
        template<typename T1>
        Quaternion operator/=(const std::complex<T1>& y) {

            T n2 = y.real() * y.real() + y.imag() * y.imag();
            T at = _a * y.real() + _b * y.imag();
            T bt = -_a * y.imag() + _b * y.real();
            T ct = _c * y.real() - _d * y.imag();
            T dt = _c * y.imag() + _d * y.real();

            _a = at / n2;
            _b = bt / n2;
            _c = ct / n2;
            _d = dt / n2;

            return *this;
        }

        /**
         * Unary +=.
         */
        template<typename T1>
        Quaternion operator+=(const Quaternion<T1>& y) {
            _a += y.a();
            _b += y.b();
            _c += y.c();
            _d += y.d();
            return *this;
        }

        /**
         * Unary -=.
         */
        template<typename T1>
        Quaternion operator-=(const Quaternion<T1>& y) {
            _a -= y._a;
            _b -= y._b;
            _c -= y._c;
            _d -= y._d;
            return *this;
        }

        /**
         * Unary multiplication.
         * 28 operations
         */
        template<typename T1>
        Quaternion operator*=(const Quaternion<T1>& y) {

            T at = _a * y.a() - _b * y.b() - _c * y.c() - _d * y.d();
            T bt = _a * y.b() + _b * y.a() + _c * y.d() - _d * y.c();
            T ct = _a * y.c() - _b * y.d() + _c * y.a() + _d * y.b();
            T dt = _a * y.d() + _b * y.c() - _c * y.b() + _d * y.a();

            _a = at;
            _b = bt;
            _c = ct;
            _d = dt;

            return *this;
        }

        /**
         * Unary division with other Quaternion.
         *
         * Warning: if the norm of y is zero, the result is
         * 4 NaNs, but maybe it should be inf.
         */
        template<typename T1>
        Quaternion operator/=(const Quaternion<T1>& y) {

            T n2 = y.norm_squared();

            T at = _a * y.a() + _b * y.b() + _c * y.c() + _d * y.d();
            T bt = -_a * y.b() + _b * y.a() - _c * y.d() + _d * y.c();
            T ct = -_a * y.c() + _b * y.d() + _c * y.a() - _d * y.b();
            T dt = -_a * y.d() - _b * y.c() + _c * y.b() + _d * y.a();

            _a = at / n2;
            _b = bt / n2;
            _c = ct / n2;
            _d = dt / n2;

            return *this;
        }

    private:
        T _a, _b, _c, _d; // the full state for a Quaternion
    };

/**
 * Predefined Quaternions on floats.
 */
    typedef Quaternion<float> Qf;
    const Qf Qf_0;
    const Qf Qf_1(1);
    const Qf Qf_i(0, 1);
    const Qf Qf_j(0, 0, 1);
    const Qf Qf_k(0, 0, 0, 1);

/**
 * Predefined Quaternions on doubles.
 */
    typedef Quaternion<double> Qd;
    const Qd Qd_0;
    const Qd Qd_1(1);
    const Qd Qd_i(0, 1);
    const Qd Qd_j(0, 0, 1);
    const Qd Qd_k(0, 0, 0, 1);

/**
 * Predefined Quaternions on long doubles.
 */
    typedef Quaternion<long double> Qld;
    const Qld Qld_0;
    const Qld Qld_1(1);
    const Qld Qld_i(0, 1);
    const Qld Qld_j(0, 0, 1);
    const Qld Qld_k(0, 0, 0, 1);

    template<typename T>
    inline Quaternion<T> spherical(T rho, T theta, T phi1, T phi2) {

        T d = std::sin(phi2);
        T cr = std::cos(phi2);
        T c = cr * std::sin(phi1);
        cr *= std::cos(phi1);
        T b = cr * std::sin(theta);
        T a = cr * std::cos(theta);

        return {rho * a, rho * b, rho * c, rho * d};
    }

    template<typename T>
    inline Quaternion<T> semipolar(T rho, T alpha, T theta1, T theta2) {

        T ca = std::cos(alpha);
        T sa = std::sin(alpha);
        T a = ca * std::cos(theta1);
        T b = ca * std::sin(theta1);
        T c = sa * std::cos(theta2);
        T d = sa * std::sin(theta2);

        return {rho * a, rho * b, rho * c, rho * d};
    }

    template<typename T>
    inline Quaternion<T> multipolar(T rho1, T theta1, T rho2, T theta2) {

        T a = rho1 * std::cos(theta1);
        T b = rho1 * std::sin(theta1);
        T c = rho2 * std::cos(theta2);
        T d = rho2 * std::sin(theta2);

        return {a, b, c, d};
    }

    template<typename T>
    inline Quaternion<T> cylindrospherical(T t, T radius, T longitude, T latitude) {

        T cl = std::cos(latitude);
        T b = radius * cl * std::cos(longitude);
        T c = radius * cl * std::sin(longitude);
        T d = radius * std::sin(latitude);

        return {t, b, c, d};
    }

    template<typename T>
    inline Quaternion<T> cylindrical(T r, T angle, T h1, T h2) {

        T a = r * std::cos(angle);
        T b = r * std::sin(angle);

        return {a, b, h1, h2};
    }

/**
 * The polar representation of a Quaternion.
 */
    template<typename T>
    using polar_representation = std::array<T, 5>;

/**
 * The polar representation of a Quaternion.
 * Returns 5 numbers:
 * - the Euclidean norm of the Quaternion,
 * - the polar angle theta,
 * - and each of the components of the "unreal unit direction".
 */
    template<typename T>
    inline polar_representation<T> to_polar_representation(const Quaternion<T>& x) {
        T nu = x.unreal_norm_squared();
        T n = std::sqrt(nu + x.a() * x.a());
        assert(nu >= 0);
        if (nu > 0) {
            T theta = std::acos(x.a() / n);
            T ns = sqrt(nu);
            return {{n, theta, x.b() / ns, x.c() / ns, x.d() / ns}};
        }
        const T pi = std::atan2(+0., -0.);
        // theta = 0 or pi, because n = +/- a().
        return {{n, n == x.a() ? 0 : pi, 0, 0, 0}};
    }

/**
 * Type used for 2x2 complex matrix representations of Quaternions.
 */
    template<typename T>
    using complex_matrix_2d = std::array<std::array<std::complex<T>, 2>, 2>;

/**
 * Returns a 2x2 complex matrix representation of a Quaternion x:
 * [ a + b i,  c + d i]
 * [ -c + d i, a - b i]
 */
    template<typename T>
    inline complex_matrix_2d<T> to_complex_matrix_2d(const Quaternion<T>& x) {
        complex_matrix_2d<T> cm;
        cm[0][0] = {x.a(), x.b()}; cm[0][1] = {x.c(), x.d()};
        cm[1][0] = -conj(cm[0][1]); cm[1][1] = conj(cm[0][0]);
        return cm;
    }

/**
 * Returns a Quaternion from a 2x2 complex matrix cm:
 * [ a + b i,  c + d i]
 * [ -c + d i, a - b i]
 */
    template <typename T>
    inline Quaternion<T> from_complex_matrix_2d(const complex_matrix_2d<T>& cm) {
        assert(cm[1][1] == conj(cm[0][0]) && cm[1][0] == -conj(cm[0][1]));
        return {cm[0][0].real(), cm[0][0].imag(), cm[0][1].real(), cm[0][1].imag()};
    }

/**
 * Type used for 4x4 real matrix representation of a quaternion.
 */
    template <typename T>
    using real_matrix_4d = std::array<std::array<T,4>,4>;

/**
 * Returns a 4x4 real matrix representation from a Quaternion x:
 * [ a  b  c  d ]
 * [-b  a -d  c ]
 * [-c  d  a -b ]
 * [-d -c  b  a ]
 */
    template <typename T>
    inline real_matrix_4d<T> to_real_matrix_4d(const Quaternion<T>& x) {
        real_matrix_4d<T> rm;
        rm[0] = {{x.a(), x.b(), x.c(), x.d()}};
        rm[1] = {{-x.b(), x.a(), -x.d(), x.c()}};
        rm[2] = {{-x.c(), x.d(), x.a(), -x.b()}};
        rm[3] = {{-x.d(), -x.c(), x.b(), x.a()}};
        return rm;
    }

/**
 * Returns a Quaternion from a 4x4 real matrix rm:
 * [ a  b  c  d ]
 * [-b  a -d  c ]
 * [-c  d  a -b ]
 * [-d -c  b  a ]
 */
    template <typename T>
    inline Quaternion<T> from_real_matrix_4d(const real_matrix_4d<T>& rm) {
        // TODO: asserts?
        return {rm[0][0],rm[0][1],rm[0][2],rm[0][3]};
    }

/**
 * A 3D rotation matrix.
 */
    template<typename T>
    using rotation_matrix = std::array<std::array<T, 3>, 3>;

/**
 * Returns a 3D rotation matrix.
 * This is the "homogeneous" expression to convert to a rotation matrix,
 * which works if the Quaternoin is not a unit Quaternion.
 */
    template<typename T>
    inline rotation_matrix<T> to_rotation_matrix(const Quaternion<T>& x) {
        // 21 operations?
        T a2 = x.a() * x.a(), b2 = x.b() * x.b(), c2 = x.c() * x.c(), d2 = x.d() * x.d();
        T ab = x.a() * x.b(), ac = x.a() * x.c(), ad = x.a() * x.d();
        T bc = x.b() * x.c(), bd = x.b() * x.d();
        T cd = x.c() * x.d();
        std::array<T, 3> r0{{a2 + b2 - c2 - d2, 2 * (bc - ad), 2 * (bd + ac)}};
        std::array<T, 3> r1{{2 * (bc + ad), a2 - b2 + c2 - d2, 2 * (cd - ab)}};
        std::array<T, 3> r2{{2 * (bd - ac), 2 * (cd + ab), a2 - b2 - c2 + d2}};
        return {{r0, r1, r2}};
    }

    template<typename T>
    inline Quaternion<T> from_rotation_matrix(const rotation_matrix<T>& rm) {
        T t = rm[0][0] + rm[1][1] + rm[2][2];
        if (t > 0) {
            T s = 0.5 / std::sqrt(t + 1);
            return {0.25 / s,
                    (rm[2][1] - rm[1][2]) * s,
                    (rm[0][2] - rm[2][0]) * s,
                    (rm[1][0] - rm[0][1]) * s};
        } else {
            if (rm[0][0] > rm[1][1] && rm[0][0] > rm[2][2]) {
                T s = 2.0 * std::sqrt(1.0 + rm[0][0] - rm[1][1] - rm[2][2]);
                return {(rm[2][1] - rm[1][2]) / s,
                        0.25 * s,
                        (rm[0][1] + rm[1][0]) / s,
                        (rm[0][2] + rm[2][0]) / s};
            } else if (rm[1][1] > rm[2][2]) {
                T s = 2.0 * std::sqrt(1.0 + rm[1][1] - rm[0][0] - rm[2][2]);
                return {(rm[0][2] - rm[2][0]) / s,
                        (rm[0][1] + rm[1][0]) / s,
                        0.25 * s,
                        (rm[1][2] + rm[2][1]) / s};
            } else {
                T s = 2.0 * std::sqrt(1.0 + rm[2][2] - rm[0][0] - rm[1][1]);
                return {(rm[1][0] - rm[0][1]) / s,
                        (rm[0][2] + rm[2][0]) / s,
                        (rm[1][2] + rm[2][1]) / s,
                        0.25 * s};
            }
        }
    }

/**
 * Returns three Euler angles {yaw, pitch, roll} in radians.
 * x is required to be a unit quaternion.
 *
 * WARNING: conversion to/from Euler angles is not ready.
 */
    template <typename T>
    inline std::array<T, 3> to_euler(const Quaternion<T>& x, T eps = 1e-12) {
        assert(x.is_unit(eps));
        const T pi = 3.14159265358979323846;
        T v = x.b()*x.c()+x.a()*x.d();
        if (std::abs(v - 0.5) < eps) {
            return {{2*atan2(x.b(),x.a()), +pi/2, 0}};
        }
        if (std::abs(v + 0.5) < eps) {
            return {{-2*atan2(x.b(),x.a()), -pi/2, 0}};
        }
        return {{atan2(2*(x.a()*x.c() - x.b()*x.d()), 1-2*(x.c()*x.c()+x.d()*x.d())),
                        std::asin(2*v),
                        atan2(2*(x.a()*x.b()-x.c()*x.d()), 1-2*(x.b()*x.b()+x.d()*x.d()))}};
    }

/**
 * Returns a unit quaternion corresponding to the three Euler angles
 * {yaw, pitch, roll} expressed in radians.
 * The conventions used are with the 3,2,1 convention ??? TODO: verify
 */
    template <typename T>
    inline Quaternion<T> from_euler(const std::array<T, 3>& x) {
        T c0 = std::cos(x[0]/2), s0 = std::sin(x[0]/2);
        T c1 = std::cos(x[1]/2), s1 = std::sin(x[1]/2);
        T c2 = std::cos(x[2]/2), s2 = std::sin(x[2]/2);
        T c0c1 = c0*c1, s0s1 = s0*s1, s0c1 = s0*c1, c0s1 = c0*s1;
        return {c0c1*c2+s0s1*s2,s0c1*c2-c0c1*s2,c0s1*c2+s0c1*s2,c0c1*s2-s0s1*c2};
    }

/**
 * Hash of a quaternion - that makes it possible to use quaternions
 * as keys in std::set/std::map, if ever needed.
 *
 * fash-hash.
 *
 * TODO: try xxhash
 * TODO: is that useful??
 * TODO: provide lexicographic order on quaternions?
 */
    template <typename T>
    struct hash : public std::unary_function<Quaternion<T>, size_t> {

        inline size_t operator()(const Quaternion<T>& x) const {

            auto mix = [](uint64_t h) {
                (h) ^= (h) >> 23;
                (h) *= 0x2127599bf4325c37ULL;
                (h) ^= (h) >> 47;
                return h;
            };

            const uint64_t len = 4 * sizeof(T); // in bytes, Qf is 4*4, Qd is 4*8, and Qld 4*16
            const uint64_t m = 0x880355f21e6d1965ULL;
            const uint64_t *pos = (const uint64_t *) &x;
            const uint64_t *end = pos + (len / 8);
            uint64_t h = 31 ^(len * m);
            uint64_t v;

            while (pos != end) {
                v = *pos++;
                h ^= mix(v);
                h *= m;
            }

            return mix(h);
        }
    };

/**
 * Lexicographic order on quaternions, which is a total order, but not compatible
 * with the field structure.
 */
    template <typename T>
    struct lexicographic_order : std::binary_function<Quaternion<T>, Quaternion<T>, bool> {
        inline constexpr bool operator()(const Quaternion<T>& x, const Quaternion<T>& y) const {
            return x.a() < y.a()
                   || (x.a() == y.a() && x.b() < y.b())
                   || (x.a() == y.a() && x.b() == y.b() && x.c() < y.c())
                   || (x.a() == y.a() && x.b() == y.b() && x.c() == y.c() && x.d() < y.d());
        }
    };

/** +
 * Returns the conjugate of x, as a new Quaternion (x is unchanged).
 */
    template<typename T>
    inline Quaternion<T> conj(const Quaternion<T>& x) {
        return {x.a(), -x.b(), -x.c(), -x.d()};
    }

/** +
 * Norms on a Quaternion.
 */
    template<typename T>
    inline T norm_squared(const Quaternion<T>& x) {
        return x.norm_squared();
    }

// abs = l2 norm = euclidean norm
    template<typename T>
    inline T abs(const Quaternion<T>& x) {
        return x.abs();
    }

    template<typename T>
    inline T unreal_norm_squared(const Quaternion<T>& x) {
        return x.unreal_norm_squared();
    }

// Hamming
    template<typename T>
    inline T norm_l0(const Quaternion<T>& x) {
        return (x.a() != 0) + (x.b() != 0) + (x.c() != 0) + (x.d() != 0);
    }

// l1 norm = taxicab = manhattan
    template<typename T>
    inline T norm_l1(const Quaternion<T>& x) {
        return std::abs(x.a()) + std::abs(x.b()) + std::abs(x.c()) + std::abs(x.d());
    }

    template<typename T, typename T1>
    inline T norm_lk(const Quaternion<T>& x, T1 k) {
        return std::pow(std::pow(std::abs(x.a()), k)
                        + std::pow(std::abs(x.b()), k)
                        + std::pow(std::abs(x.c()), k)
                        + std::pow(std::abs(x.d()), k), 1.0 / k);
    }

// norm sup = max norm = norm inf
    template<typename T>
    inline T norm_sup(const Quaternion<T>& x) {
        return std::max(std::max(std::abs(x.a()), std::abs(x.b())),
                        std::max(std::abs(x.c()), std::abs(x.d())));
    }

/**
 * Quaternion tests.
 */

    template <typename T, typename T1 =T>
    inline bool is_zero(const Quaternion<T>& x, T1 eps = 0) {
        return x.is_zero(eps);
    }

    template <typename T, typename T1 =T>
    inline bool is_non_zero(const Quaternion<T>& x, T1 eps = 0) {
        return x.is_non_zero(eps);
    }

    template <typename T>
    inline bool is_nan(const Quaternion<T>& x) {
        return x.is_nan();
    }

    template <typename T>
    inline bool is_inf(const Quaternion<T>& x) {
        return x.is_inf();
    }

    template <typename T>
    inline bool is_finite(const Quaternion<T>& x) {
        return x.is_finite();
    }

    template<typename T, typename T1 =T>
    inline bool is_unit(const Quaternion<T>& x, T1 eps = 0) {
        return x.is_unit(eps);
    }

    template<typename T, typename T1 =T>
    inline bool is_real(const Quaternion<T>& x, T1 eps = 0) {
        return x.is_real(eps);
    }

    template<typename T, typename T1 =T>
    inline bool is_complex(const Quaternion<T>& x, T1 eps = 0) {
        return x.is_complex(eps);
    }

    template<typename T, typename T1 =T>
    inline bool is_unreal(const Quaternion<T>& x, T1 eps = 0) {
        return x.is_unreal(eps);
    }

/**
 * Equality:
 * - Quaternion <-> real
 * - Quaternion <-> complex
 * - Quaternion <-> quaternion
 */

/**
 * Quaternion <-> real
 *
 * For these equality operators, the goal is to allow e.g.
 * Qd(0,0,0,0) == 0, which is very nice to read, and pretty unambiguous,
 * but we must take care that the compiler doesn't try to use these with
 * T2 being an "exotic" type that we wouldn't want. That's why we have
 * "IS_CONVERTIBLE".
 *
 * operator== returns false if the Quaternion is not real.
 * If the Quaternion is real, it returns true if x.a() == y.
 */
    template <typename T, typename T2, IS_CONVERTIBLE(T2, T)>
    inline bool operator==(const Quaternion<T>& x, T2 y) {
        return x.is_real() && x.a() == y;
    }

    template <typename T, typename T2, IS_CONVERTIBLE(T2, T)>
    inline bool operator==(T2 y, const Quaternion<T>& x) {
        return x == y;
    }

    template <typename T, typename T2, IS_CONVERTIBLE(T2, T)>
    inline bool operator!=(const Quaternion<T>& x, T2 y) {
        return !(x == y);
    }

    template <typename T, typename T2, IS_CONVERTIBLE(T2, T)>
    inline bool operator!=(T2 y, const Quaternion<T>& x) {
        return !(x == y);
    }

    template <typename T, typename T2, typename T3, IS_CONVERTIBLE(T2, T), IS_CONVERTIBLE(T3,T)>
    inline bool nearly_equal(const Quaternion<T>& x, T2 y, T3 eps) {
        return x.is_real() && is_nearly_equal(x.a(), y, eps);
    }

    template <typename T, typename T2, typename T3, IS_CONVERTIBLE(T2, T), IS_CONVERTIBLE(T3,T)>
    inline bool nearly_equal(T2 y, const Quaternion<T>& x, T3 eps) {
        return x.is_real() && is_nearly_equal(x.a(), y, eps);
    }

/**
 * Quaternion <-> std::complex
 *
 * operator== here returns false if the Quaternion is not complex.
 */
    template<typename T, typename T2, IS_CONVERTIBLE(T2,T)>
    inline bool operator==(const Quaternion<T>& x, const std::complex<T2>& y) {
        return is_complex(x) && x.a() == y.real() && x.b() == y.imag();
    }

    template<typename T, typename T2, IS_CONVERTIBLE(T2,T)>
    inline bool operator!=(const Quaternion<T>& x, const std::complex<T2>& y) {
        return !(x == y);
    }

// Same, swapping the lhs and rhs.
    template<typename T, typename T2, IS_CONVERTIBLE(T2,T)>
    inline bool operator==(const std::complex<T2>& y, const Quaternion<T>& x) {
        return x == y;
    }

    template<typename T, typename T2, IS_CONVERTIBLE(T2,T)>
    inline bool operator!=(const std::complex<T2>& y, const Quaternion<T>& x) {
        return x != y;
    }

    template <typename T, typename T2, typename T3, IS_CONVERTIBLE(T2, T), IS_CONVERTIBLE(T3,T)>
    inline bool nearly_equal(const Quaternion<T>& x, const std::complex<T2>& y, T3 eps) {
        return is_complex(x, eps)
               && is_nearly_equal(x.a(), y.real(), eps)
               && is_nearly_equal(x.b(), y.imag(), eps);
    }

    template <typename T, typename T2, typename T3, IS_CONVERTIBLE(T2, T), IS_CONVERTIBLE(T3,T)>
    inline bool nearly_equal(const std::complex<T2>& y, const Quaternion<T>& x, T3 eps) {
        return nearly_equal(x, y, eps);
    }

/**
 * Quaternion <-> Quaternion
 */
    template<typename T1, typename T2>
    inline bool operator==(const Quaternion<T1>& x, const Quaternion<T2>& y) {
        return x.a() == y.a() && x.b() == y.b() && x.c() == y.c() && x.d() == y.d();
    }

    template<typename T>
    inline bool operator!=(const Quaternion<T>& x, const Quaternion<T>& y) {
        return !(x == y);
    }

    template <typename T1, typename T2, typename T3>
    inline bool nearly_equal(const Quaternion<T1>& x, const Quaternion<T2>& y, T3 eps) {
        return is_nearly_equal(x.a(), y.a(), eps)
               && is_nearly_equal(x.b(), y.b(), eps)
               && is_nearly_equal(x.c(), y.c(), eps)
               && is_nearly_equal(x.d(), y.d(), eps);
    }

// TODO: equality of Quaternion and complex, of Quaternion and array/container
// TODO: y by const ref?

    template<typename T, typename T1>
    inline Quaternion<T> operator+(const Quaternion<T>& x, T1 y) {
        return Quaternion<T>(x) += y;
    }

    template<typename T, typename T1>
    inline Quaternion<T> operator+(T1 y, const Quaternion<T>& x) {
        return x + y;
    }

    template<typename T>
    inline Quaternion<T> operator+(const Quaternion<T>& x, std::complex<T>& y) {
        return Quaternion<T>(x) += y;
    }

    template<typename T>
    inline Quaternion<T> operator+(std::complex<T>& y, const Quaternion<T>& x) {
        return x + y;
    }

    template<typename T>
    inline Quaternion<T> operator+(const Quaternion<T>& x, const Quaternion<T>& y) {
        return Quaternion<T>(x) += y;
    }

    template<typename T, typename T1>
    inline Quaternion<T> operator-(const Quaternion<T>& x, T1 y) {
        return Quaternion<T>(x) -= y;
    }

    template<typename T, typename T1>
    inline Quaternion<T> operator-(T1 y, const Quaternion<T>& x) {
        return Quaternion<T>(x) += -y;
    }

    template<typename T>
    inline Quaternion<T> operator-(const Quaternion<T>& x, std::complex<T>& y) {
        return Quaternion<T>(x) -= y;
    }

    template<typename T>
    inline Quaternion<T> operator-(std::complex<T>& y, const Quaternion<T>& x) {
        return x + (-y);
    }

    template<typename T>
    inline Quaternion<T> operator-(const Quaternion<T>& x, const Quaternion<T>& y) {
        return Quaternion<T>(x) -= y;
    }

/**
 * SSE operations: tried 2 implementations (SO and vectorclass): not faster.
 * Boost: as fast as boost implementation.
 */
    template<typename T, typename T1>
    inline Quaternion<T> operator*(const Quaternion<T>& x, T1 y) {
        return Quaternion<T>(x) *= y;
    }

    template<typename T, typename T1>
    inline Quaternion<T> operator*(T1 y, const Quaternion<T>& x) {
        return x * y;
    }

    template<typename T>
    inline Quaternion<T> operator*(const Quaternion<T>& x, std::complex<T>& y) {
        return Quaternion<T>(x) *= y;
    }

    template<typename T>
    inline Quaternion<T> operator*(std::complex<T>& y, const Quaternion<T>& x) {
        return x * y;
    }

    template<typename T>
    inline Quaternion<T> operator*(const Quaternion<T>& x, const Quaternion<T>& y) {
        return Quaternion<T>(x) *= y;
    }

    template<typename T>
    inline Quaternion<T> inverse(const Quaternion<T>& x) {
        return conj(x) / norm_squared(x);
    }

    template<typename T, typename T1>
    inline Quaternion<T> operator/(const Quaternion<T>& x, T1 y) {
        return Quaternion<T>(x) /= y;
    }

    template<typename T, typename T1>
    inline Quaternion<T> operator/(T1 y, const Quaternion<T>& x) {
        return y * inverse(x);
    }

    template<typename T>
    inline Quaternion<T> operator/(const Quaternion<T>& x, std::complex<T>& y) {
        return Quaternion<T>(x) /= y;
    }

    template<typename T>
    inline Quaternion<T> operator/(std::complex<T>& y, const Quaternion<T>& x) {
        return y * inverse(x);
    }

    template<typename T>
    inline Quaternion<T> operator/(const Quaternion<T>& x, const Quaternion<T>& y) {
        return x * inverse(y);
    }

    template<typename T>
    inline T dot(const Quaternion<T>& x, const Quaternion<T>& y) {
        return x.a() * y.a() + x.b() * y.b() + x.c() * y.c() + x.d() * y.d();
    }

/**
 * 9 operations
 */
    template<typename T>
    inline Quaternion<T> cross(const Quaternion<T>& x, const Quaternion<T>& y) {
        return {0,
                x.c() * y.d() - x.d() * y.c(),
                x.d() * y.b() - x.b() * y.d(),
                x.b() * y.c() - x.c() * y.b()};
    }

    template<typename T>
    inline Quaternion<T> commutator(const Quaternion<T>& x, const Quaternion<T>& y) {
        return x * y - y * x;
    }

    template<typename T>
    inline Quaternion<T> normalize(const Quaternion<T>& x) {
        assert(abs(x) > 0); // or this is not normalizable
        return x / abs(x);
    }

/**
 * Exponential of a Quaternion.
 * This code seems to be quite a bit faster than boost, while giving
 * the same results. Boost uses a Taylor approximation for sinc,
 * which *might* (not sure) be why they are slightly slower here.
 *
 * exp(log(x)) == x always, but log(exp(x)) != x is already not true
 * for complex number, because the log is multi-valued.
 *
 * NOTE: the precision is not great with so many floating point operations
 */
    template<typename T>
    inline Quaternion<T> exp(const Quaternion<T>& x) {
        T un = x.unreal_norm_squared();
        if (un == 0)
            return {std::exp(x.a())};
        T n1 = std::sqrt(un);
        T ea = std::exp(x.a());
        T n2 = ea * std::sin(n1) / n1; // TODO: this needs speed-up (compared to boost)
        return {ea * std::cos(n1), n2 * x.b(), n2 * x.c(), n2 * x.d()};
    }

/**
 * Log of a Quaternion.
 * exp(log(x)) == x always, but log(exp(x)) != x is already not true
 * for complex number, because the log is multi-valued.
 *
 * NOTE: the precision is not great with so many floating point operations
 */
    template<typename T>
    inline Quaternion<T> log(const Quaternion<T>& x) {
        T nu2 = x.unreal_norm_squared();
        if (nu2 == 0) {
            if (x.a() > 0)
                return {std::log(x.a())};
            else { // TODO: reduce number of instructions
                std::complex<T> l = log(std::complex<T>(x.a(), 0)); // TODO: is that correct?
                return {l.real(), l.imag()};
            }
        }
        T a = x.a();
        assert(nu2 > 0);
        T n = std::sqrt(a * a + nu2); // n > 0
        T th = std::acos(a / n) / std::sqrt(nu2); // -a <= a/n <= 1
        return {std::log(n), th * x.b(), th * x.c(), th * x.d()}; // n > 0
    }

// TODO: log2, log10

/**
 * 10 operations:
 * a^2 - b^2 - c^2 - d^2
 * 2 a b
 * 2 a c
 * 2 a d
 */
    template<typename T>
    inline Quaternion<T> pow2(const Quaternion<T>& x) {
        T aa = 2 * x.a();
        return {x.a() * x.a() - x.unreal_norm_squared(),
                aa * x.b(),
                aa * x.c(),
                aa * x.d()};
    }

/**
 * 14 operations:
 * a (a^2 - 3 (b^2 + c^2 + d^2))
 * -b (-3 a^2 + b^2 + c^2 + d^2)
 * -c (-3 a^2 + b^2 + c^2 + d^2)
 * -d (-3 a^2 + b^2 + c^2 + d^2)
 */
    template<typename T>
    inline Quaternion<T> pow3(const Quaternion<T>& x) {
        T a2 = x.a() * x.a();
        T n1 = x.unreal_norm_squared();
        T n2 = 3 * a2 - n1;
        return {x.a() * (a2 - 3 * n1),
                x.b() * n2,
                x.c() * n2,
                x.d() * n2};
    }

/**
 * 18 operations:
 * a^4 - 6 a^2 (b^2 + c^2 + d^2) + (b^2 + c^2 + d^2)^2
 * -4 a b (-a^2 + b^2 + c^2 + d^2)
 * -4 a c (-a^2 + b^2 + c^2 + d^2)
 * -4 a d (-a^2 + b^2 + c^2 + d^2)
 */
    template<typename T>
    inline Quaternion<T> pow4(const Quaternion<T>& x) {
        T a2 = x.a() * x.a();
        T n1 = x.unreal_norm_squared();
        T n2 = 4 * x.a() * (a2 - n1);
        return {a2 * a2 - 6 * a2 * n1 + n1 * n1,
                x.b() * n2,
                x.c() * n2,
                x.d() * n2};
    }

// TODO: this needs to be redone, with powi, powf, and a switcher on type like std::pow
/**
 * I benchmarked that method written via the polar representation,
 * and it turned out to be much slower, and lexicographic_order numerically stable,
 * than this implementation. This implementation is also much faster
 * than the boost implementation. However, via the polar representation
 * I could compute pow for any real exponent, whereas this method is
 * limited to integer exponents.
 */
    template<typename T>
    inline Quaternion<T> pow(const Quaternion<T>& x, int expt) {

        if (expt < 0)
            return inverse(pow(x, -expt));
        if (expt == 0)
            return {1};
        if (expt == 1)
            return x;
        if (expt == 2)
            return pow2(x);
        if (expt == 3)
            return pow3(x);
        if (expt == 4)
            return pow4(x);

        Quaternion<T> x4 = pow4(x), y = x4;
        for (size_t i = 1; i < expt / 4; ++i)
            y *= x4;
        if (expt % 4 == 3)
            y *= pow3(x);
        if (expt % 4 == 2)
            y *= pow2(x);
        if (expt % 4 == 1)
            y *= x; std::pow(3,4);
        return y;
    }

/**
 * Real power of a Quaternion.
 */
    template <typename T, typename T1, typename std::enable_if<!std::is_same<T1, int>()>::type* = nullptr>
    inline Quaternion<typename std::common_type<T,T1>::type> pow(const Quaternion<T>& x, T1 a) {
        Quaternion<typename std::common_type<T,T1>::type> xx(x); // TODO: needed?
        if (std::floor(a) == a)
            return pow(x, (int) a); // KEEP, really helps with numerical accuracy for e.g. pow(Q,Q)
        return exp(a * log(xx));
    }

/**
 * Quaternion power of a Quaternion.
 * (that should cover all the other cases for the exponent...)
 * TODO: test against pow just above
 */
    template<typename T>
    inline Quaternion<T> pow(const Quaternion<T>& x, const Quaternion<T>& a) {
        if (a.is_real())
            return pow(x, a.a());
        return exp(a * log(x));
    }

// TODO: sqrt
    template<typename T>
    inline Quaternion<T> cos(const Quaternion<T>& x) {
        T z = x.unreal_norm_squared();
        if (z == 0)
            return {std::cos(x.a())};
        // if (x.real() == 0) return {1};
        z = std::sqrt(z);
        T w = -std::sin(x.real()) * std::sinh(z) / z; // z > 0 by construction
        return {std::cos(x.real()) * std::cosh(z), w * x.b(), w * x.c(), w * x.d()};
    }

    template<typename T>
    inline Quaternion<T> sin(const Quaternion<T>& x) {
        T z = x.unreal_norm_squared();
        if (z == 0)
            return {std::sin(x.a())};
        // if (x.real() == 0) return {0};
        z = std::sqrt(z);
        T w = std::cos(x.real()) * std::sinh(z) / z; // z > 0 by construction
        return {std::sin(x.real()) * std::cosh(z), w * x.b(), w * x.c(), w * x.d()};
    }

/**
 * TODO: optimize instruction count
 * TODO: reciprocals
 * TODO: check /0
 */
    template<typename T>
    inline Quaternion<T> tan(const Quaternion<T>& x) {
        T z = x.unreal_norm_squared();
        if (z == 0)
            return {std::tan(x.a())};
        z = std::sqrt(z);
        T n = std::sinh(2 * z);
        T d = std::cos(2 * x.a()) + std::cosh(2 * z); // hmmm, I couldn't find a case with real a that would yield d = 0
        T r = n / (z * d);
        return {std::sin(2 * x.a()) / d, r * x.b(), r * x.c(), r * x.d()};
    }

    template<typename T>
    inline Quaternion<T> cosh(const Quaternion<T>& x) {
        return (exp(x) + exp(-x)) / 2;
    }

    template<typename T>
    inline Quaternion<T> sinh(const Quaternion<T>& x) {
        return (exp(x) - exp(-x)) / 2;
    }

    template<typename T>
    inline Quaternion<T> tanh(const Quaternion<T>& x) {
        return sinh(x) / cosh(x); // TODO: cosh never zero for quaternions too?
    }

/**
 * result = a*x + b*y
 * TODO: k1,k2 can't be complex!
 */
    template<typename T, typename K>
    inline Quaternion<T> axby(K k1, const Quaternion<T>& x, K k2, const Quaternion<T>& y) {

        T a = k1 * x.a() + k2 * y.a();
        T b = k1 * x.b() + k2 * y.b();
        T c = k1 * x.c() + k2 * y.c();
        T d = k1 * x.d() + k2 * y.d();

        return {a, b, c, d};
    }

// TODO: operator<< and operator>>
} // end namespace quaternion

#endif //QUATERNIONS_QUATERNION_H


#endif //STRUCCLUST_QUATERNION_H
