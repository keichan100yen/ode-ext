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
    std::cout.precision(20);
	using itv = kv::interval<double>;

    typedef Lorenz Func;
    ub::vector<itv> x1(3); // This is for odelong_maffine2.
    x1(0) = 15.; x1(1) = 15.; x1(2) = 36.;
    auto x2 = kv::make_init<itv>(15., 15., 36.); // This is for odelong_maffine_ext.
    itv begin(0.), end1(20.), end2(20.); 
    
    /*
    typedef Duffing Func;
    ub::vector<itv> x1(2);
    x1(0) = itv(-5., 5.); x1(1) = itv(-5., 5.);
    auto x2 = kv::make_init<itv>(itv(-5., 5.), itv(-5., 5.));
    itv begin(0.), end1, end2; 
    */

    /*
    typedef VdP Func;
    ub::vector<itv> x1(2);
    x1(0) = itv(1.); x1(1) = itv(1.);
    auto x2 = kv::make_init<itv>(1., 1.);
    itv begin(0.), end1, end2; 
    */

    std::cout << "odelong_maffine2 begin (original solver of the kv library)." << std::endl;
    auto ts = std::chrono::system_clock::now();
    int r = kv::odelong_maffine2(Func(), x1, begin, end1, kv::ode_param<double>().set_iteration(0).set_verbose(0).set_order(22).set_autostep(true).set_ep_reduce(250).set_ep_reduce_limit(300));

    auto te = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = te-ts;
    std::cout << "time: " << elapsed_seconds.count() << " [sec]" << std::endl;

    if (r == 0)
		std::cout << "Cannot calculate solution.\n";
    else if (r == 1) {
		std::cout << "Solution calculated until t = " << end1 << ".\n";
        std::cout << x1 << "\n";
    } else {
		std::cout << "Solution calculated. t = " << end1 << "\n";
        std::cout << x1 << "\n";
	}
    std::cout << "odelong_maffine2 end." << std::endl;
    // 
    
    std::cout << std::endl << std::endl;
    
    std::cout << "odelong_maffine_ext begin (ODEWithExpressionTemplate)." << std::endl;
    ts = std::chrono::system_clock::now();
    r = kv::odelong_maffine_ext(Func(), x2, begin, end2, kv::ode_param<double>().set_iteration(0).set_verbose(0).set_order(22).set_autostep(true).set_ep_reduce(250).set_ep_reduce_limit(300));
    te = std::chrono::system_clock::now();
    elapsed_seconds = te-ts;
    std::cout << "time: " << elapsed_seconds.count() << " [sec]" << std::endl;

    if (r == 0)
        std::cout << "Cannot calculate solution.\n";
    else if (r == 1) {
        std::cout << "Solution calculated until t = " << end2 << ".\n";
        for (auto r : x2)
            std::cout << r << ", ";
        std::cout << std::endl;;
    } else {
        std::cout << "Solution calculated. t = " << end2 << "\n";
        for (auto r : x2)
            std::cout << r << ", ";
        std::cout << std::endl;
    }
    std::cout << "odelong_maffine_ext End." << std::endl;

	return 0;
}
