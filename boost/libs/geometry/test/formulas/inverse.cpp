// Boost.Geometry
// Unit Test

// Copyright (c) 2016 Oracle and/or its affiliates.

// Contributed and/or modified by Adam Wulkiewicz, on behalf of Oracle

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_formula.hpp"
#include "inverse_cases.hpp"

#include <boost/geometry/formulas/vincenty_inverse.hpp>
#include <boost/geometry/formulas/thomas_inverse.hpp>
#include <boost/geometry/formulas/andoyer_inverse.hpp>

template <typename Result>
void check_inverse(Result const& result, expected_result const& expected, expected_result const& reference, double reference_error)
{
    check_one(result.distance, expected.distance, reference.distance, reference_error);
    check_one(result.azimuth, expected.azimuth, reference.azimuth, reference_error, true);
    check_one(result.reverse_azimuth, expected.reverse_azimuth, reference.reverse_azimuth, reference_error, true);
    check_one(result.reduced_length, expected.reduced_length, reference.reduced_length, reference_error);
    check_one(result.geodesic_scale, expected.geodesic_scale, reference.geodesic_scale, reference_error);
}

void test_all(expected_results const& results)
{
    double const d2r = bg::math::d2r<double>();
    double const r2d = bg::math::r2d<double>();

    double lon1r = results.p1.lon * d2r;
    double lat1r = results.p1.lat * d2r;
    double lon2r = results.p2.lon * d2r;
    double lat2r = results.p2.lat * d2r;

    // WGS84
    bg::srs::spheroid<double> spheroid(6378137.0, 6356752.3142451793);

    bg::formula::result_inverse<double> result;

    typedef bg::formula::vincenty_inverse<double, true, true, true, true, true> vi_t;
    result = vi_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid);
    result.azimuth *= r2d;
    result.reverse_azimuth *= r2d;
    check_inverse(result, results.vincenty, results.karney, 0.00000001);

    typedef bg::formula::thomas_inverse<double, true, true, true, true, true> th_t;
    result = th_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid);
    result.azimuth *= r2d;
    result.reverse_azimuth *= r2d;
    check_inverse(result, results.thomas, results.karney, 0.0000001);

    typedef bg::formula::andoyer_inverse<double, true, true, true, true, true> an_t;
    result = an_t::apply(lon1r, lat1r, lon2r, lat2r, spheroid);
    result.azimuth *= r2d;
    result.reverse_azimuth *= r2d;
    check_inverse(result, results.andoyer, results.karney, 0.0001);
}

int test_main(int, char*[])
{
    for (size_t i = 0; i < expected_size; ++i)
    {
        test_all(expected[i]);
    }

    return 0;
}
