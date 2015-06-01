#include <iostream>
#include <chrono>

#include <kv/ode-maffine2.hpp>
#include "ode-maffine-ext.hpp"

namespace ub = boost::numeric::ublas;

struct VdP {
    template <typename... T, typename Time> auto operator()(const std::tuple<T...>& x, const Time& t) {
        auto &x0 = std::get<0>(x);
        auto &x1 = std::get<1>(x);

        auto y0 = x1;
        auto y1 = 0.25_i * (1._i - x0 * x0) * x1 - x0;

        return std::tie(y0, y1);
    }

    template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
        ub::vector<T> y(2);

        y(0) = x(1);
        y(1) = 0.25 * (1. - x(0)*x(0))*x(1) - x(0);

        return y;
    }
};

struct Lorenz {
    template <typename... T, typename Time> auto operator()(const std::tuple<T...>& x, const Time& t) {
        auto &x0 = std::get<0>(x);
        auto &x1 = std::get<1>(x);
        auto &x2 = std::get<2>(x);

        auto y0 = 10._i * (x1 - x0);
        auto y1 = (28._i - x2) * x0 - x1;
        auto y2 = -(8._i/3._i) * x2 + x0 * x1;

        return std::tie(y0, y1, y2);
    }

    template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
        ub::vector<T> y(3);

        y(0) = 10. * (x(1) - x(0));
        y(1) = (28.-x(2)) * x(0) - x(1);
        y(2) = (-kv::interval<double>(8.)/3.) * x(2) + x(0) * x(1);

        return y;
    }
};

struct Duffing {
    template <typename... T, typename Time> auto operator()(const std::tuple<T...>& x, const Time& t) {
        auto &x0 = std::get<0>(x);
        auto &x1 = std::get<1>(x);

        auto y0 = x1;
        auto y1 = -0.1_i * x1 - x0 * x0 * x0 + 3.5_i * cos(t);

        return std::tie(y0, y1);
    }

    template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
        ub::vector<T> y(2);

        y(0) = x(1);
        y(1) = -0.1 * x(1) - x(0) * x(0) * x(0) + 3.5 * cos(t);

        return y;
    }
};

int main()
{
	using itv = kv::interval<double>;

    typedef Lorenz Func;
    ub::vector<itv> x(3);
    x(0) = 15.; x(1) = 15.; x(2) = 36.;
    auto x2 = kv::make_init<itv>(15., 15., 36.);

    /*
    typedef Duffing Func;
    ub::vector<itv> x(2);
    x(0) = itv(-5., 5.); x(1) = itv(-5., 5.);
    auto x2 = kv::make_init<itv>(itv(-5., 5.), itv(-5., 5.));
    */

    /*
    typedef VdP Func;
    ub::vector<itv> x(2);
    x(0) = itv(1.); x(1) = itv(1.);
    auto x2 = kv::make_init<itv>(1., 1.);
    */

    std::cout.precision(20);

    itv end;
    end = 20.;
    // end = 0.01;

    int r;

    auto ts = std::chrono::system_clock::now();
    // r = kv::odelong_maffine(Func(), x, itv(0.), end, kv::ode_param<double>().set_iteration(0).set_verbose(0).set_order(22).set_autostep(false).set_ep_reduce(250).set_ep_reduce_limit(300));
    r = kv::odelong_maffine2(Func(), x, itv(0.), end, kv::ode_param<double>().set_iteration(0).set_verbose(0).set_order(22).set_autostep(true).set_ep_reduce(250).set_ep_reduce_limit(300));

    auto te = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = te-ts;
    std::cout << "time: " << elapsed_seconds.count() << " [sec]" << std::endl;

    if (r == 0)
		std::cout << "Cannot calculate solution.\n";
    else if (r == 1) {
		std::cout << "Solution calculated until t = " << end << ".\n";
        std::cout << x << "\n";
    } else {
		std::cout << "Solution calculated.\n";
        std::cout << x << "\n";
	}

    std::cout << std::endl << std::endl;


    end = 20.;
    /*
    end = 0.01; */

    ts = std::chrono::system_clock::now();
    r = kv::odelong_maffine_ext(Func(), x2, itv(0.), end, kv::ode_param<double>().set_iteration(0).set_verbose(0).set_order(22).set_autostep(true).set_ep_reduce(250).set_ep_reduce_limit(300));
    // r = kv::odelong_maffine_ext<3>(Func(), x, itv(0.), end, kv::ode_param<double>().set_iteration(0).set_verbose(0).set_order(22).set_autostep(false).set_ep_reduce(250).set_ep_reduce_limit(300));
    te = std::chrono::system_clock::now();
    elapsed_seconds = te-ts;
    std::cout << "time: " << elapsed_seconds.count() << " [sec]" << std::endl;

    if (r == 0)
        std::cout << "Cannot calculate solution.\n";
    else if (r == 1) {
        std::cout << "Solution calculated until t = " << end << ".\n";
        for (auto r : x2)
            std::cout << r << ", ";
        std::cout << std::endl;;
    } else {
        std::cout << "Solution calculated.\n";
        for (auto r : x2)
            std::cout << r << ", ";
        std::cout << std::endl;
    }

	return 0;
}
