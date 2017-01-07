Copyright (c) 2017 Keiichiro KASHIWAGI

#ifndef ODE_MAFFINE_EXT_HPP
#define ODE_MAFFINE_EXT_HPP

#include <tuple>
#include <array>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/affine.hpp>
#include <kv/autodif.hpp>
#include <kv/ode-param.hpp>

#include "psa-ext.hpp"

namespace kv {
extern void *enabler;
namespace ub = boost::numeric::ublas;

template <int N, int E, bool F, typename T, typename... Seq> struct decl_psa_ext_impl;
template <int N, int E, typename T, typename... Seq>
struct decl_psa_ext_impl<N, E, false, T, Seq...> {
    typedef std::tuple<Seq...> type;
};
template <int N, int E, typename T, typename... Seq>
struct decl_psa_ext_impl<N, E, true, T, Seq...> {
    typedef typename decl_psa_ext_impl<N+1, E, N+1<E, T, psa_ext<N, T>, Seq...>::type type;
};
template <int S, int E, typename T> struct decl_psa_ext  {
    typedef typename decl_psa_ext_impl<S, E, true, T>::type type;
};

template <int N, typename... T, typename std::enable_if<N == sizeof...(T)-1>::type *& = enabler>
auto integral_impl(const std::tuple<T...>& x) {
    return std::make_tuple(integral(std::get<N>(x)));
}
template <int N, typename... T, typename std::enable_if<N != sizeof...(T)-1>::type *& = enabler>
auto integral_impl(const std::tuple<T...>& x) {
    return std::tuple_cat(std::make_tuple(integral(std::get<N>(x))), integral_impl<N+1>(x));
}
template <typename... T> auto integral(const std::tuple<T...>& x) {
    return integral_impl<0>(x);
}

template <int N, typename T, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void initialize_impl(int n) {}
template <int N, typename T, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void initialize_impl(int n) {
    std::tuple_element<N, T>::type ::init(n);
    initialize_impl<N+1, T>(n);
}
template <typename T> void initialize(int n) {
    initialize_impl<0, T>(n);
}

template <int N, typename T, typename R, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void assign_impl(const R& r) {}
template <int N, typename T, typename R, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void assign_impl(const R& r) {
    std::tuple_element<N, T>::type ::assign(r(N));
    assign_impl<N+1, T>(r);
}
template <typename T, typename R> void assign(const R& r) {
    assign_impl<0, T>(r);
}
template <int N, typename T, typename R, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void assign_impl(int n) {}
template <int N, typename T, typename R, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void assign_impl(int n) {
    std::tuple_element<N, T>::type ::template assign<typename std::tuple_element<N, R>::type>(n);
    assign_impl<N+1, T, R>(n);
}
template <typename T, typename R> void assign(int n) {
    assign_impl<0, T, R>(n);
}

template <int N, typename T, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void eval_begin_impl() {}
template <int N, typename T, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void eval_begin_impl() {
    std::tuple_element<N, T>::type ::eval_begin();
    eval_begin_impl<N+1, T>();
}
template <typename T> void eval_begin() {
    eval_begin_impl<0, T>();
}
template <int N, typename T, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void eval_type1_impl(int n) {}
template <int N, typename T, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void eval_type1_impl(int n) {
    std::tuple_element<N, T>::type ::eval_type1(n);
    eval_type1_impl<N+1, T>(n);
}
template <typename T> void eval_type1(int n) {
    eval_type1_impl<0, T>(n);
}
template <int N, typename T, typename Domain, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void eval_type2_impl(const Domain& domain, int n) {}
template <int N, typename T, typename Domain, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void eval_type2_impl(const Domain& domain, int n) {
    std::tuple_element<N, T>::type ::template eval_type2(domain, n);
    eval_type2_impl<N+1, T>(domain, n);
}
template <typename T, typename Domain> void eval_type2(const Domain& domain, int n) {
    eval_type2_impl<0, T>(domain, n);
}

template <int N, typename T, typename R1, typename R2, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void copy_data_impl() {}
template <int N, typename T, typename R1, typename R2, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void copy_data_impl() {
    std::tuple_element<N, T>::type ::template copy_from_autodif<typename std::tuple_element<N, R1>::type>();
    std::tuple_element<N, T>::type ::data[0] = std::tuple_element<N, R2>::type ::data[0].v;
    copy_data_impl<N+1, T, R1, R2>();
}
template <typename T, typename R1, typename R2> void copy_data() {
    copy_data_impl<0, T, R1, R2>();
}

template <int N, typename T, typename std::enable_if<N == std::tuple_size<T>::value-1>::type *& = enabler>
auto norm_max_impl(int n) {
    return std::abs(mid(std::tuple_element<N, T>::type ::v(n)));
}
template <int N, typename T, typename std::enable_if<N != std::tuple_size<T>::value-1>::type *& = enabler>
auto norm_max_impl(int n) {
    return std::max(std::abs(mid(std::tuple_element<N, T>::type ::v(n))), norm_max_impl<N+1, T>(n));
}
template <typename T> auto norm_max(int n) {
    return norm_max_impl<0, T>(n);
}

template <int N, typename T, typename R, typename std::enable_if<N == std::tuple_size<T>::value-1>::type *& = enabler>
auto rad_max_impl(int n) {
    return abs(std::tuple_element<N, T>::type ::v(n) - std::tuple_element<N, R>::type ::v(n)).upper();
}
template <int N, typename T, typename R, typename std::enable_if<N != std::tuple_size<T>::value-1>::type *& = enabler>
auto rad_max_impl(int n) {
    return std::max((abs(std::tuple_element<N, T>::type ::v(n) - std::tuple_element<N, R>::type ::v(n)).upper()), rad_max_impl<N+1, T, R>(n));
}
template <typename T, typename R> auto rad_max(int n) {
    return rad_max_impl<0, T, R>(n);
}

template <int N, typename T, typename R, typename std::enable_if<N == std::tuple_size<T>::value>::type *& = enabler>
void expand_impl(int n, const R& rad) {}
template <int N, typename T, typename R, typename std::enable_if<N != std::tuple_size<T>::value>::type *& = enabler>
void expand_impl(int n, const R& rad){
    std::tuple_element<N, T>::type ::v(n) += rad;
    expand_impl<N+1, T>(n, rad);
}
template <typename T, typename R> void expand(int n, const R& rad) {
    expand_impl<0, T>(n, rad);
}

template <int N, typename T, typename R, typename std::enable_if<N == std::tuple_size<T>::value-1>::type *& = enabler>
auto exist_test_impl(int n) {
    return subset(std::tuple_element<N, T>::type ::v(n), std::tuple_element<N, R>::type ::v(n));
}
template <int N, typename T, typename R, typename std::enable_if<N != std::tuple_size<T>::value-1>::type *& = enabler>
auto exist_test_impl(int n){
    if (! subset(std::tuple_element<N, T>::type ::v(n), std::tuple_element<N, R>::type ::v(n)))
        return false;
    return exist_test_impl<N+1, T, R>(n);
}
template <typename T, typename R> auto exist_test(int n) {
    return exist_test_impl<0, T, R>(n);
}

template <int N, typename C, typename Iad, typename E, typename Res, typename Init, typename D, typename std::enable_if<N == std::tuple_size<C>::value>::type *& = enabler>
void eval_result_impl(Res &result, const Init& init, const D& deltat, const D& deltat_n) {}
template <int N, typename C, typename Iad, typename E, typename Res, typename Init, typename D, typename std::enable_if<N != std::tuple_size<C>::value>::type *& = enabler>
void eval_result_impl(Res &result, const Init& init, const D& deltat, const D& deltat_n) {
    auto fc = eval_ext<typename std::tuple_element<N, C>::type>(deltat);       // f(c)
    auto fi_dif = eval_ext<typename std::tuple_element<N, Iad>::type>(deltat); // F'(I)
    auto tmp = fi_dif.d(0) * init(0);                                    // F'(I)*(v-c)
    for (int i=1; i<std::tuple_size<C>::value; ++i)
        tmp += fi_dif.d(i) * init(i);
    result(N) = fc + std::tuple_element<N, E>::type ::v(std::tuple_element<N, E>::type ::data.size()-1) * deltat_n + tmp;   // f(c) + F'(I)(v-c) + g(I)

    eval_result_impl<N+1, C, Iad, E>(result, init, deltat, deltat_n);
}
template <typename C, typename Iad, typename E, typename Res, typename Init, typename D>
void eval_result(Res &result, const Init& init, const D& deltat, const D& deltat_n) {
    eval_result_impl<0, C, Iad, E>(result, init, deltat, deltat_n);
}

template <int Dim, class F, class T>
int ode_maffine_ext(
        F f,
        ub::vector<affine<T>>& init,
        const interval<T>& start,
        interval<T>& end,
        ode_param<T> p = ode_param<T>()
) {
    ub::vector<interval<T>> c(Dim), I(Dim);
    ub::vector<autodif<interval<T>>> Iad;

    ub::vector<affine<T>>  result(Dim);

    int maxnum_save;

    affine<T> s1, s2;
    interval<T> s2i;
    ub::vector<affine<T>> s1_save;
    ub::vector<interval<T>> s2i_save;

    int ret_val;
    interval<T> end2, deltat;

    typedef time_ext<interval<T>> TIME_EXT;
    typedef typename decl_psa_ext<0, Dim, interval<T>>::type C_EXT;
    typedef typename decl_psa_ext<0, Dim, autodif<interval<T>>>::type Iad_EXT;
    typedef typename decl_psa_ext<Dim, Dim+Dim, interval<T>>::type E_EXT;

    TIME_EXT time_ext;
    C_EXT c_ext;
    Iad_EXT iad_ext;
    E_EXT e_ext;

    auto result_c_ext = integral(f(c_ext, time_ext));
    auto result_iad_ext = integral(f(iad_ext, time_ext));
    auto result_e_ext = integral(f(e_ext, time_ext));

    typedef decltype(result_c_ext) RESULT_C_EXT;
    typedef decltype(result_iad_ext) RESULT_Iad_EXT;
    typedef decltype(result_e_ext) RESULT_E_EXT;

    initialize<C_EXT>(p.order);
    initialize<Iad_EXT>(p.order);
    initialize<E_EXT>(p.order+1);

    initialize<RESULT_C_EXT>(p.order);
    initialize<RESULT_Iad_EXT>(p.order);
    initialize<RESULT_E_EXT>(p.order+1);

    for (int i=0; i<Dim; ++i) {
        I(i) = to_interval(init(i));
        c(i) = mid(I(i));
    }

    Iad = autodif<interval<T>>::init(I);

    TIME_EXT::assign(start);
    assign<C_EXT>(c);
    assign<Iad_EXT>(Iad);

    // Type1
    for (int i=0; i<p.order-1; ++i) {
        eval_begin<RESULT_C_EXT>();
        eval_begin<RESULT_Iad_EXT>();

        eval_type1<RESULT_C_EXT>(i);
        eval_type1<RESULT_Iad_EXT>(i);

        assign<C_EXT, RESULT_C_EXT>(i+1);
        assign<Iad_EXT, RESULT_Iad_EXT>(i+1);
    }

    copy_data<RESULT_E_EXT, RESULT_Iad_EXT, Iad_EXT>();

    eval_begin<RESULT_E_EXT>();
    eval_type1<RESULT_E_EXT>(p.order-1);
    assign<E_EXT, RESULT_E_EXT>(p.order);

    // Type2
#if 1
    T radius;
    int n_rad;
    T m = 1.;
    for (int i=0; i<Dim; i++) {
        m = std::max(m, norm(I(i)));
    }
    T tolerance = m * p.epsilon;

    if (p.autostep) {
        // use two non-zero coefficients of higher order term
        radius = 0.;
        n_rad = 0;
        for (int j = p.order; j>=1; j--) {
            m = norm_max<E_EXT>(j);
            if (m == 0.)
                continue;
            auto radius_tmp = std::pow((double)m, 1./j);
            // std::cout << j << " " << m << " " << radius_tmp << "\n";
            if (radius_tmp > radius)
                radius = radius_tmp;
            n_rad++;
            if (n_rad == 2)
                break;
        }
        radius = std::pow((double)tolerance, 1./p.order) / radius;
    }

    if (p.autostep) {
        end2 = mid(start + radius);
        if (end2 >= end.lower()) {
            end2 = end;
            ret_val = 2;
        } else
            ret_val = 1;
    } else {
        end2 = end;
        ret_val = 2;
    }
    deltat = end2 - start;
#endif

    auto domain = interval<T>::hull(0., deltat);
    eval_begin<RESULT_E_EXT>();
    eval_type2<RESULT_E_EXT>(domain, p.order);

    T rad = rad_max<RESULT_E_EXT, E_EXT>(p.order);
    rad *= 2.;
    expand<E_EXT>(p.order, interval<T>(-rad, rad));

    eval_begin<RESULT_E_EXT>();
    eval_type2<RESULT_E_EXT>(domain, p.order);

    if (! exist_test<RESULT_E_EXT, E_EXT>(p.order))
        return 0;

    if (p.ep_reduce == 0)
        maxnum_save = affine<T>::maxnum();

    // f(c) + F'(I)(v-c) + g(I)
    auto deltat_n = deltat;
    for (int i=0; i<p.order-1; ++i)
        deltat_n *= deltat;

    //v-c
    for (int i=0; i<Dim; ++i)
        init(i).a[0] = 0.;

    eval_result<C_EXT, Iad_EXT, RESULT_E_EXT>(result, init, deltat, deltat_n);

    if (p.ep_reduce == 0) {
        s1_save.resize(Dim);
        s2i_save.resize(Dim);
        for (int i=0; i<Dim; i++) {
            split(result(i), maxnum_save, s1, s2);
            s2i = to_interval(s2);
            s1_save(i) = s1;
            s2i_save(i) = s2i;
        }

        affine<T>::maxnum() = maxnum_save;
        for (int i=0; i<Dim; i++) {
            s1_save(i).resize();
            result(i) = append(s1_save(i), (affine<T>)s2i_save(i));
        }
    } else
        epsilon_reduce(result, p.ep_reduce, p.ep_reduce_limit);

    init = result;
    if (ret_val == 1)
        end = end2;

    return ret_val;
}

template <int Dim, class F, class T>
int odelong_maffine_ext(
    F f,
    ub::vector<affine<T>>& init,
    const interval<T>& start,
    interval<T>& end,
    ode_param<T> p = ode_param<T>(),
    ub::matrix<interval<T>>* mat = nullptr
) {
    ub::vector<affine<T>> x, x1;
    interval<T> t, t1;
    int r;
    int ret_val = 0;

    x = init;
    t = start;
    p.set_autostep(true);
    while (1) {
        x1 = x;
        t1 = end;

        r = ode_maffine_ext<Dim>(f, x1, t, t1, p);
        if (r == 0) {
            if (ret_val == 1) {
                init = x1;
                end = t;
            }
            return ret_val;
        }
        ret_val = 1;
        if (p.verbose == 1) {
            std::cout << "t: " << t1 << "\n";
            std::cout << to_interval(x1) << "\n";
        }

        if (r == 2) {
            init = x1;
            return 2;
        }
        t = t1;
        x = x1;
    }
}

template <int Dim, class F, class T>
int odelong_maffine_ext(
    F f,
    ub::vector<interval<T>>& init,
    const interval<T>& start,
    interval<T>& end,
    ode_param<T> p = ode_param<T>()
) {
    auto end2 = end;

    int maxnum_save = affine<T>::maxnum();
    affine<T>::maxnum() = 0;

    ub::vector<affine<T>> x = init;

    int r = odelong_maffine_ext<Dim>(f, x, start, end2, p);

    affine<T>::maxnum() = maxnum_save;

    if (r == 0)
        return 0;

    for (int i=0; i<(int)init.size(); ++i)
        init(i) = to_interval(x(i));
    if (r == 1)
        end = end2;

    return r;
}
template <unsigned long N, class F, class T>
int odelong_maffine_ext(
    F f,
    std::array<interval<T>, N>& init,
    const interval<T>& start,
    interval<T>& end,
    ode_param<T> p = ode_param<T>()
)  {
    ub::vector<interval<T>> init2(N);
    for (int i=0; i<N; ++i)
        init2(i) = init[i];
    int r = odelong_maffine_ext<N>(f, init2, start, end, p);
    for (int i=0; i<N; ++i)
        init[i] = init2(i);
    return r;
}
template <typename T, typename... R>
auto make_init(const R&... x) {
    return std::array<T, sizeof...(R)> { T(x)... };
}

} // namespace kv

#endif
