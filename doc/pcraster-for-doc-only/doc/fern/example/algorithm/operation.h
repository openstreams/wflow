// -----------------------------------------------------------------------------
// Fern © Geoneric
//
// This file is part of Geoneric Fern which is available under the terms of
// the GNU General Public License (GPL), version 2. If you do not want to
// be bound by the terms of the GPL, you may purchase a proprietary license
// from Geoneric (http://www.geoneric.eu/contact).
// -----------------------------------------------------------------------------
#pragma once
#include <utility>
// Include the relevant traits before including the algorithms.
#include "fern/core/data_customization_point/scalar.h"
#include "fern/example/algorithm/raster_traits.h"
#include "fern/example/algorithm/raster_argument_traits.h"
#include "fern/example/algorithm/raster_customization_point.h"
#include "fern/algorithm/core/cast.h"
#include "fern/algorithm/space/focal/slope.h"
#include "fern/algorithm/algebra/elementary/add.h"
#include "fern/core/math.h"


namespace example {

// Execution policy to use by the algorithms: sequential or parallel.
extern fern::algorithm::ExecutionPolicy execution_policy;


template<
    typename Value1,
    typename Value2>
fern::algorithm::add::result_type<Value1, Value2> add(
    Value1 const& lhs,
    Value2 const& rhs)
{
    assert(fern::is_equal(cell_size(lhs, 0), cell_size(lhs, 1)));
    assert(fern::is_equal(cell_size(lhs, 0), cell_size(rhs, 0)));
    assert(fern::is_equal(cell_size(lhs, 1), cell_size(rhs, 1)));
    assert(size(lhs, 0) == size(rhs, 0));
    assert(size(lhs, 1) == size(rhs, 1));

    fern::algorithm::add::result_type<Value1, Value2> result(
        cell_size(lhs, 0), size(lhs, 0), size(lhs, 1));

    fern::algorithm::algebra::add(execution_policy, lhs, rhs, result);

    return result;
}


template<
    typename ResultValueType,
    typename Value>
fern::CloneT<Value, ResultValueType> cast(
    Value const& value)
{
    assert(fern::is_equal(cell_size(value, 0), cell_size(value, 1)));

    fern::CloneT<Value, ResultValueType> result(
        cell_size(value, 0),
        size(value, 0),
        size(value, 1));

    fern::algorithm::core::cast(execution_policy, value, result);

    return result;
}


template<
    typename Value>
fern::CloneT<Value, fern::value_type<Value>> slope(
    Value const& value)
{
    assert(fern::is_equal(cell_size(value, 0), cell_size(value, 1)));

    fern::CloneT<Value, fern::value_type<Value>> result(
        cell_size(value, 0), size(value, 0), size(value, 1));

    fern::algorithm::space::slope(execution_policy, value, result);

    return result;
}

} // namespace example
