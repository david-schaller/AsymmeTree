#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

/**
 * Template class for square matrices.
 *
 * This template class is an implementation of a 2-dimensional square matrix with arbitrary elements.
 */
template<typename T>
class Matrix {
public:
  /**
   * Initializes a matrix with a default value.
   * @param dim dimension of the quare matrix (number of rows/columns).
   * @param value default value to fill the matrix with.
   */
  void
  initMatrix(size_t dim, T value);

  /**
   * Returns the dimension of the matrix.
   */
  size_t
  getDim() const { return m_dim; };

  /**
   * Returns the value at a specific position.
   * @param row number of the row.
   * @param col number of the column.
   */
  T&
  at(size_t row, size_t col);

  /**
   * Returns the value at a specific position (for const matrices).
   * @param row number of the row.
   * @param col number of the column.
   */
  const T&
  at(size_t row, size_t col) const;

private:
  size_t m_dim;                             /*!< dimension of the square matrix */
  std::vector<std::vector<T>> m_matrix;     /*!< matrix */
};

template<typename T>
void
Matrix<T>::initMatrix(size_t dim, T value)
{
  m_dim = dim;
  m_matrix = std::vector<std::vector<T>>(m_dim, std::vector<T>(m_dim, value));
}

template<typename T>
T&
Matrix<T>::at(size_t row, size_t col)
{
    if(col >= m_dim || row >= m_dim){
       throw std::out_of_range("Index out of bounds! Please wait for help.");
    }

    return m_matrix[row][col];
}

template<typename T>
const T&
Matrix<T>::at(size_t row, size_t col) const
{
    if(col >= m_dim || row >= m_dim){
       throw std::out_of_range("Index out of bounds! Please wait for help.");
    }

    return m_matrix[row][col];
}

#endif /* MATRIX_H */
