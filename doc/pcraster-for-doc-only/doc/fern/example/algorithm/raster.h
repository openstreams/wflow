// -----------------------------------------------------------------------------
// Fern © Geoneric
//
// This file is part of Geoneric Fern which is available under the terms of
// the GNU General Public License (GPL), version 2. If you do not want to
// be bound by the terms of the GPL, you may purchase a proprietary license
// from Geoneric (http://www.geoneric.eu/contact).
// -----------------------------------------------------------------------------
#pragma once
#include <cstddef>
#include <vector>


namespace example {

/*!
    @brief      Simple raster template class.
    @tparam     T Type of the values in the raster.
*/
template<
    typename T>
class Raster
{

public:

    using iterator = typename std::vector<T>::iterator;

                   Raster              (double cell_size,
                                        size_t nr_rows,
                                        size_t nr_cols);

                   Raster              (Raster&& other);

                   ~Raster             ()=default;

    Raster<T>&     operator=           (Raster<T>&& other);

    iterator        begin              ();

    iterator        end                ();

    size_t          index              (size_t row,
                                        size_t col) const;

    T&              get                (size_t index);

    T const &       get                (size_t index) const;

    T&              get                (size_t row,
                                        size_t col);

    T const&        get                (size_t row,
                                        size_t col) const;

    double          cell_size          () const;

    size_t          nr_rows            () const;

    size_t          nr_cols            () const;

private:

    double          _cell_size;

    size_t          _nr_rows;

    size_t          _nr_cols;

    std::vector<T> _values;

};


/*!
    @brief      Constructor.
*/
template<
    typename T>
Raster<T>::Raster(
    double cell_size,
    size_t nr_rows,
    size_t nr_cols)

    : _cell_size(cell_size),
      _nr_rows(nr_rows),
      _nr_cols(nr_cols),
      _values(nr_rows * nr_cols)

{
}


/*!
    @brief      Move constructor.
    @param      other Raster to move from.
*/
template<
    typename T>
Raster<T>::Raster(Raster&& other)

    : _cell_size(other._cell_size),
      _nr_rows(other._nr_rows),
      _nr_cols(other._nr_cols),
      _values(std::move(other._values))

{
}


/*!
    @brief      Move assignment operator.
    @param      other Raster to move from.
*/
template<
    typename T>
Raster<T>& Raster<T>::operator=(
    Raster<T>&& other)
{
    _cell_size = other._cell_size;
    _nr_rows = other._nr_rows;
    _nr_cols = other._nr_cols;
    _values = std::move(other._values);

    return *this;
}


template<
    typename T>
inline typename Raster<T>::iterator Raster<T>::begin()
{
    return _values.begin();
}


template<
    typename T>
inline typename Raster<T>::iterator Raster<T>::end()
{
    return _values.end();
}


template<
    typename T>
inline size_t Raster<T>::index(
    size_t row,
    size_t col) const
{
    return row * _nr_cols + col;
}


template<
    typename T>
inline T& Raster<T>::get(
    size_t index)
{
    return _values[index];
}


template<
    typename T>
inline T const& Raster<T>::get(
    size_t index) const
{
    return _values[index];
}


template<
    typename T>
inline T& Raster<T>::get(
    size_t row,
    size_t col)
{
    return get(index(row, col));
}


template<
    typename T>
inline T const& Raster<T>::get(
    size_t row,
    size_t col) const
{
    return get(index(row, col));
}


template<
    typename T>
double Raster<T>::cell_size() const
{
    return _cell_size;
}


template<
    typename T>
size_t Raster<T>::nr_rows() const
{
    return _nr_rows;
}


template<
    typename T>
size_t Raster<T>::nr_cols() const
{
    return _nr_cols;
}

} // namespace example
