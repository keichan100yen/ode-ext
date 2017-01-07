#ifndef PSA_EXT_HPP
#define PSA_EXT_HPP

#include <type_traits>
#include <kv/psa.hpp>
/*
#if defined(_WIN64) || defined(_WIN32)
#define __tls __declspec(thread)
#else
#define __tls __thread
#endif
*/

namespace kv {

extern void *enabler;
namespace ub = boost::numeric::ublas;

struct ext_base {
    static void init(const size_t& s) {}
    static void eval_begin() {}
    static void eval_type1(int) {}
    template <typename Domain> static void eval_type2(const Domain&, int) {}
    template <typename AutoDif>
    static void copy_from_autodif() {}
};
struct const_ext_base {};

template <int N, typename T>
struct psa_ext : public ext_base {
    typedef T type;

    static ub::vector<T> data;
    static void init(const size_t& s) { data.resize(s); }

    static T& v(int i) { return data[i]; }

    static void assign(const T& value) {
        data[0] = value;
    }
    template <typename R>
    static void assign(int i) {
        data[i] = R::v(i);
    }
    template <typename AutoDif>
    static void copy_from_autodif() {
        for (int i=0; i<AutoDif::data.size(); ++i)
            data[i] = AutoDif::data[i].v;
    }
};
template <int N, typename T> ub::vector<T> psa_ext<N, T>::data;

template <typename T>
struct time_ext : public ext_base {
    typedef T type;

    static T data;

    static const T v(int i) {
        switch (i) {
        case 0:
            return data;
        case 1:
            return T(1.);
        default:
            return T(0.);
        }
    }

    static void assign(const T& value) {
        data = value;
    }
};
template <typename T> T time_ext<T>::data;

template <typename T, char... C>
struct const_ext : public ext_base, const_ext_base {
    typedef T type;

    static T data;

    static const T v(int i) {
        return (i==0) ? data : T(0.);
    }

    static void assign(const T& value) {
        data = value;
    }
};
template <typename T, char... C> T const_ext<T, C...>::data;

template <typename T, typename Op, typename R>
struct unary_ext : public ext_base {
    typedef T type;
    typedef R Rhs;

    static bool is_calculated;
    static ub::vector<type> data;
    static void init(const size_t& s) {
        R::init(s);
        is_calculated = false;
        data.resize(s);
    }

    static type& v(int i) { return data[i]; }

    static void eval_begin() {
        R::eval_begin();
        is_calculated = false;
    }

    static void eval_type1(int i) {
        R::eval_type1(i);
        Op::template eval_type1<unary_ext<T, Op, R>, R>(i);
    }
    template <typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        R::template eval_type2(domain, i);
        Op::template eval_type2<unary_ext<T, Op, R>, R>(domain, i);
    }

    template <typename AutoDif, typename std::enable_if<std::is_same<typename R::type, typename AutoDif::type>::value>::type *& = enabler>
    static void copy_from_autodif() {}
    template <typename AutoDif, typename std::enable_if<!std::is_same<typename R::type, typename AutoDif::type>::value>::type *& = enabler>
    static void copy_from_autodif() {
        R::template copy_from_autodif<typename AutoDif::Rhs>();
        for (int i=0; i<AutoDif::data.size(); ++i)
            data[i] = AutoDif::data[i].v;
    }
};
template <typename T, typename Op, typename R> bool unary_ext<T, Op, R>::is_calculated;
template <typename T, typename Op, typename R> ub::vector<T> unary_ext<T, Op, R>::data;

template <typename T, typename Op, typename L, typename R>
struct binary_ext : public ext_base {
    typedef T type;
    typedef L Lhs;
    typedef R Rhs;
    static bool is_calculated;
    static ub::vector<type> data;
    static void init(const size_t& s) {
        L::init(s);
        R::init(s);
        is_calculated = false;
        data.resize(s);
    }

    static type& v(int i) { return data[i]; }

    static void eval_begin() {
        L::eval_begin();
        R::eval_begin();
        is_calculated = false;
    }

    static void eval_type1(int i) {
        L::eval_type1(i);
        R::eval_type1(i);
        Op::template eval_type1<binary_ext<T, Op, L, R>, L, R>(i);
    }
    template <typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        L::template eval_type2(domain, i);
        R::template eval_type2(domain, i);
        Op::template eval_type2<binary_ext<T, Op, L, R>, L, R>(domain, i);
    }

    template <typename AutoDif, typename std::enable_if<std::is_same<typename R::type, typename AutoDif::type>::value>::type *& = enabler>
    static void copy_from_autodif() {}
    template <typename AutoDif, typename std::enable_if<!std::is_same<typename R::type, typename AutoDif::type>::value>::type *& = enabler>
    static void copy_from_autodif() {
        L::template copy_from_autodif<typename AutoDif::Lhs>();
        R::template copy_from_autodif<typename AutoDif::Rhs>();
        for (int i=0; i<AutoDif::data.size(); ++i)
            data[i] = AutoDif::data[i].v;
    }
};
template <typename T, typename Op, typename L, typename R> bool binary_ext<T, Op, L, R>::is_calculated;
template <typename T, typename Op, typename L, typename R> ub::vector<typename binary_ext<T, Op, L, R>::type> binary_ext<T, Op, L, R>::data;


struct plus_ext {
    template <typename Ret, typename R>
    static void eval_type1(int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i) = R::v(i);
    }
    template <typename Ret, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        eval_type1<Ret, R>(i);
    }
};

struct minus_ext {
    template <typename Ret, typename R>
    static void eval_type1(int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i) = -R::v(i);
    }
    template <typename Ret, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        eval_type1<Ret, R>(i);
    }
};

struct cos_ext {
    template <typename Ret, typename R>
    static void eval_type1(int i) {
        static psa<typename Ret::type> tmp;
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        psa<typename Ret::type>::mode() = 1;
        tmp.v.resize(i+1);
        for (int j=0; j<i+1; ++j)
            tmp.v(j) = R::v(j);

        Ret::v(i) = cos(tmp).v(i);
    }
    template <typename Ret, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        static psa<typename Ret::type> tmp;
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        psa<typename Ret::type>::mode() = 2;
        psa<typename Ret::type>::domain() = domain;
        tmp.v.resize(i+1);
        for (int j=0; j<i+1; ++j)
            tmp.v(j) = R::v(j);

        Ret::v(i) = cos(tmp).v(i);
    }
};
struct sin_ext {
    template <typename Ret, typename R>
    static void eval_type1(int i) {
        static psa<typename Ret::type> tmp;
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        psa<typename Ret::type>::mode() = 1;
        tmp.v.resize(i+1);
        for (int j=0; j<i+1; ++j)
            tmp.v(j) = R::v(j);

        Ret::v(i) = sin(tmp).v(i);
    }
    template <typename Ret, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        static psa<typename Ret::type> tmp;
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        psa<typename Ret::type>::mode() = 2;
        psa<typename Ret::type>::domain() = domain;
        tmp.v.resize(i+1);
        for (int j=0; j<i+1; ++j)
            tmp.v(j) = R::v(j);

        Ret::v(i) = sin(tmp).v(i);
    }
};

struct inv_ext {
    template <typename Ret, typename R>
    static void eval_type1(int i) {
        static psa<typename Ret::type> tmp;
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        psa<typename Ret::type>::mode() = 1;
        tmp.v.resize(i+1);
        for (int j=0; j<i+1; ++j)
            tmp.v(j) = R::v(j);

        Ret::v(i) = inv(tmp).v(i);
    }
    template <typename Ret, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        static psa<typename Ret::type> tmp;
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        psa<typename Ret::type>::mode() = 2;
        psa<typename Ret::type>::domain() = domain;
        tmp.v.resize(i+1);
        for (int j=0; j<i+1; ++j)
            tmp.v(j) = R::v(j);

        Ret::v(i) = inv(tmp).v(i);
    }
};

struct integral_ext {
    template <typename Ret, typename R>
    static void eval_type1(int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i+1) = R::v(i) / (double)(i+1);
    }
    template <typename Ret, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        assert(i != 0);
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i) = R::v(i-1) / (double)(i);
        Ret::v(i) += (R::v(i) / (double)(i+1)) * domain;
    }
};

struct add_ext {
    template <typename Ret, typename L, typename R>
    static void eval_type1(int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i) = L::v(i) + R::v(i);
    }
    template <typename Ret, typename L, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        eval_type1<Ret, L, R>(i);
    }
};

struct sub_ext {
    template <typename Ret, typename L, typename R>
    static void eval_type1(int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i) = L::v(i) - R::v(i);
    }
    template <typename Ret, typename L, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        eval_type1<Ret, L, R>(i);
    }
};

struct mul_ext {
    template <typename Ret, typename L, typename R>
    static void eval_type1(int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::v(i) = typename Ret::type(0.);
        for (int j=0; j<=i; ++j)
            Ret::v(i) += L::v(j) * R::v(i-j);
    }
    template <typename Ret, typename L, typename R, typename Domain>
    static void eval_type2(const Domain& domain, int i) {
        if (Ret::is_calculated)
            return;
        Ret::is_calculated = true;
        Ret::data(i) = typename Ret::type(0.);
        for (int j=2*i; j>=i; --j) {
            for (int k=j-i; k<i+1; ++k)
                Ret::v(i) += L::v(k) * R::v(j-k);
            if (j != i)
                Ret::v(i) *= domain;
        }
    }
};

template <typename R> struct is_const_ext { static const bool value = std::is_base_of<const_ext_base, R>::value; };

template <typename R, typename std::enable_if<std::is_base_of<ext_base, R>::value && !is_const_ext<R>::value>::type *& = enabler>
auto operator+(const R&) {
    return unary_ext<decltype(+ typename R::type()), plus_ext, R>();
}
template <typename RT, char...RC>
auto operator+(const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(+ RT()), '+', RC...> Cnt;
    Cnt::assign(+ const_ext<RT, RC...>::v(0));
    return Cnt();
}
template <typename R, typename std::enable_if<std::is_base_of<ext_base, R>::value && !is_const_ext<R>::value>::type *& = enabler>
auto operator-(const R&) {
    return unary_ext<decltype(- typename R::type()), minus_ext, R>();
}
template <typename RT, char...RC>
auto operator-(const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(- RT()), '-', RC...> Cnt;
    Cnt::assign(- const_ext<RT, RC...>::v(0));
    return Cnt();
}
template <typename R, typename std::enable_if<std::is_base_of<ext_base, R>::value && !is_const_ext<R>::value>::type *& = enabler>
auto cos(const R&) {
    return unary_ext<decltype(cos(typename R::type())), cos_ext, R>();
}
template <typename RT, char...RC>
auto cos(const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(cos(RT())), 'c', RC...> Cnt;
    Cnt::assign(cos(const_ext<RT, RC...>::v(0)));
    return Cnt();
}
template <typename R, typename std::enable_if<std::is_base_of<ext_base, R>::value && !is_const_ext<R>::value>::type *& = enabler>
auto sin(const R&) {
    return unary_ext<decltype(sin(typename R::type())), sin_ext, R>();
}
template <typename RT, char...RC>
auto sin(const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(sin(RT())), 's', RC...> Cnt;
    Cnt::assign(sin(const_ext<RT, RC...>::v(0)));
    return Cnt();
}
template <typename R, typename std::enable_if<std::is_base_of<ext_base, R>::value && !is_const_ext<R>::value>::type *& = enabler>
auto inv(const R&) {
    return unary_ext<typename R::type, inv_ext, R>();
}
template <typename RT, char...RC>
auto inv(const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(inv(RT())), 'i', RC...> Cnt;
    Cnt::assign(inv(const_ext<RT, RC...>::v(0)));
    return Cnt();
}

template <typename R, typename std::enable_if<std::is_base_of<ext_base, R>::value>::type *& = enabler>
auto integral(const R&) {
    return unary_ext<typename R::type, integral_ext, R>();
}

template <typename T, typename Time>
auto eval_ext(const Time& deltat) {
    int s = T::data.size();
    typename T::type r = T::v(s-1);
    for (int i=s-2; i>=0; --i)
        r = r * deltat + T::v(i);
    return r;
}

template <typename L, typename R> struct are_const_exts { static const bool value = std::is_base_of<const_ext_base, L>::value && std::is_base_of<const_ext_base, R>::value; };

template <typename L, typename R, typename std::enable_if<std::is_base_of<ext_base, L>::value && std::is_base_of<ext_base, R>::value && !are_const_exts<L, R>::value>::type *& = enabler>
auto operator+(const L&, const R&) {
    return binary_ext<decltype(typename L::type() + typename R::type()), add_ext, L, R>();
}
template <typename LT, char... LC, typename RT, char...RC>
auto operator+(const const_ext<LT, LC...>&, const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(LT() + RT()), LC..., '+', RC...> Cnt;
    Cnt::assign(const_ext<LT, LC...>::v(0) + const_ext<RT, RC...>::v(0));
    return Cnt();
}
template <typename L, typename R, typename std::enable_if<std::is_base_of<ext_base, L>::value && std::is_base_of<ext_base, R>::value && !are_const_exts<L, R>::value>::type *& = enabler>
auto operator-(const L&, const R&) {
    return binary_ext<decltype(typename L::type() - typename R::type()), sub_ext, L, R>();
}
template <typename LT, char... LC, typename RT, char...RC>
auto operator-(const const_ext<LT, LC...>&, const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(LT() - RT()), LC..., '-', RC...> Cnt;
    Cnt::assign(const_ext<LT, LC...>::v(0) - const_ext<RT, RC...>::v(0));
    return Cnt();
}
template <typename L, typename R, typename std::enable_if<std::is_base_of<ext_base, L>::value && std::is_base_of<ext_base, R>::value && !are_const_exts<L, R>::value>::type *& = enabler>
auto operator*(const L&, const R&) {
    return binary_ext<decltype(typename L::type() * typename R::type()), mul_ext, L, R>();
}
template <typename LT, char... LC, typename RT, char...RC>
auto operator*(const const_ext<LT, LC...>&, const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(LT() * RT()), LC..., '*', RC...> Cnt;
    Cnt::assign(const_ext<LT, LC...>::v(0) * const_ext<RT, RC...>::v(0));
    return Cnt();
}
template <typename L, typename R, typename std::enable_if<std::is_base_of<ext_base, L>::value && std::is_base_of<ext_base, R>::value && !are_const_exts<L, R>::value>::type *& = enabler>
auto operator/(const L&, const R&) {
    typedef unary_ext<typename R::type, inv_ext, R> Inv;
    return binary_ext<decltype(typename L::type() * typename Inv::type()), mul_ext, L, Inv>();
}
template <typename LT, char... LC, typename RT, char...RC>
auto operator/(const const_ext<LT, LC...>&, const const_ext<RT, RC...>&) {
    typedef const_ext<decltype(LT() / RT()), LC..., '/', RC...> Cnt;
    Cnt::assign(const_ext<LT, LC...>::v(0) / const_ext<RT, RC...>::v(0));
    return Cnt();
}

template <typename T>
void print_ext() {
    std::cout << "(";
    for (auto v : T::data)
        std::cout << v;
    std::cout << ")" << std::endl;
}

} // namespace kv

template <char... Chars>
kv::const_ext<kv::interval<double>, Chars...> operator "" _i() {
    std::string buf;
    for (auto c : std::initializer_list<char>{Chars...})
        buf.push_back(c);
    kv::const_ext<kv::interval<double>, Chars...> ::assign(kv::interval<double>(buf));
    return kv::const_ext<kv::interval<double>, Chars...> ();
}

#endif
