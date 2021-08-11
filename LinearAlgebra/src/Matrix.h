#pragma once

#include "ArrayAllocator.h"
#include "MatrixCommaAssignment.h"
#include "MatrixEvaluation.h"


class MatrixDimension
{
private:
    int m_Row, m_Col;

protected:
    MatrixDimension()
        : m_Row(0), m_Col(0) {}

    MatrixDimension(int row, int col)
        : m_Row(row), m_Col(col) {}

    void Swap(MatrixDimension& matDim)
    {
        std::swap(m_Row, matDim.m_Row);
        std::swap(m_Col, matDim.m_Row);
    }

    bool Resize(int row, int col)
    {
        bool flag = false;
        
        if (m_Row != row)
        {
            m_Row = row;
            flag = true;
        }

        if (m_Col != col)
        {
            m_Col = col;
            flag = true;
        }

        return flag;
    }

    void CheckRange(int row, int col) const
    {
        #ifdef _CHECK_RANGE
          assert(m_Row > row);
          assert(m_Col > col);
        #endif
    }

public:
    int RowSize() const { return m_Row; }
    int ColSize() const { return m_Col; }
    int Size() const { return m_Row * m_Col; }
};


template<typename DataType>
class Shared2DAllocator : public SharedAllocator<DataType>
{
private:
    int m_LD;  // Leading Dimension

protected:
    //Shared2DAllocator()
    //    : SharedAllocator() {}
    Shared2DAllocator() = delete;

    //Shared2DAllocator(int size)
    //    : SharedAllocator(size) {}
    Shared2DAllocator(int size) = delete;

    Shared2DAllocator(DataType* data, int ld)
        : SharedAllocator(data), m_LD(ld) {}

    void Swap(Shared2DAllocator& other)
    {
        std::swap(m_LD, other.m_LD);
        SharedAllocator<DataType>::Swap(other);
    }

    //int GetLD(int ld) const { return m_LD; }
    int GetLD() const { return m_LD; }

    enum { DEFAULT_ROWS = 0, DEFAULT_COLUMNS = 0 };
};


template<typename DataType, int NROW = 3, int NCOL = 3>
class Static2DAllocator : public StaticAllocator<DataType, NROW * NCOL>
{
protected:
    Static2DAllocator()
        : StaticAllocator() {}

    //Static2DAllocator(int size)
    //    : StaticAllocator(size) {}
    Static2DAllocator(int size) = delete;

    //Static2DAllocator(DataType* data, int ld)
    //    : StaticAllocator(data) {}
    Static2DAllocator(DataType* data, int ld) = delete;

    //void Swap(Static2DAllocator& other)
    //{
    //    assert(false);
    //}

    enum { DEFAULT_ROWS = NROW, DEFAULT_COLUMNS = NCOL };
};


// C/C++ Type Matrix Storage (Row-Major)
template<typename AllocatorType>
class CMatrixStorage : public MatrixDimension, public AllocatorType
{
public:
    using DataType = typename AllocatorType::MemberType;
    using Allocator2DType = AllocatorType;

    int GetLD() const 
    {
        //return AllocatorType::GetLD(ColSize());
        return ColSize();
    }

protected:
    using SharedMatrixType = Matrix<CMatrixStorage<Shared2DAllocator<DataType>>>;

    CMatrixStorage()
        : MatrixDimension(AllocatorType::DEFAULT_ROWS, AllocatorType::DEFAULT_COLUMNS)
        , AllocatorType() {}

    CMatrixStorage(int row, int col)
        : MatrixDimension(row, col), AllocatorType(Size()) {}

    //CMatrixStorage(DataType* data, int row, int col, int ld)
    //    : MatrixDimension(row, col), AllocatorType(data, ld) {}
    CMatrixStorage(DataType* data, int row, int col, int ld) = delete;

    void Swap(CMatrixStorage& other)
    {
        MatrixDimension::Swap(other);
        AllocatorType::Swap(other);
    }

    const DataType& GetRef(int row, int col) const
    {
        CheckRange(row, col);
        return *AllocatorType::Data(row * GetLD() + col);
    }

    DataType& GetRef(int row, int col)
    {
        CheckRange(row, col);
        return *AllocatorType::Data(row * GetLD() + col);
    }

    const DataType* GetPtr(int row, int col) const
    {
        CheckRange(row, col);
        return AllocatorType::Data(row * GetLD() + col);
    }

    DataType* GetPtr(int row, int col)
    {
        CheckRange(row, col);
        return AllocatorType::Data(row * GetLD() + col);
    }

    Vector<SharedAllocator<DataType>> GetColumn(const Range& rowRange, int col) = delete;
    //Vector<SharedAllocator<DataType>> GetColumn(const Range& rowRange, int col)
    //{
    //    assert(false);
    //    return GetRow(col, rowRange);
    //}

    // Get Row Vector with Specific Range
    // GetRow(0, Range::All()): 0-th row
    // GetRow(1, Range(0, 1)): 1-th row from column 0 to column 1
    Vector<SharedAllocator<DataType>> GetRow(int row, const Range& colRange)
    {
        return Vector<SharedAllocator<typename AllocatorType::MemberType>>(GetPtr(row, colRange.First()), colRange.Count(ColSize() - 1));
    }

public:
    const DataType* operator[](int row) const { return GetPtr(row, 0); }
    DataType* operator[](int row) { return GetPtr(row, 0); }
};


// Fortran Type Matrix Storage (Col-Major)
template<typename AllocatorType>
class FMatrixStorage : public MatrixDimension, public AllocatorType
{
public:
    using DataType = typename AllocatorType::MemberType;
    using Allocator2DType = AllocatorType;

    int GetLD() const
    {
        //return AllocatorType::GetLD(RowSize());
        return RowSize();
    }

protected:
    using SharedMatrixType = Matrix<FMatrixStorage<Shared2DAllocator<DataType>>>;

    FMatrixStorage()
        : MatrixDimension(AllocatorType::DEFAULT_ROWS, AllocatorType::DEFAULT_COLUMNS)
        , AllocatorType() {}

    FMatrixStorage(int row, int col)
        : MatrixDimension(row, col), AllocatorType(Size()) {}

    //FMatrixStorage(DataType* data, int row, int col, int ld)
    //    : MatrixDimension(row, col), AllocatorType(data, ld) {}
    FMatrixStorage(DataType* data, int row, int col, int ld) = delete;

    void Swap(FMatrixStorage& other)
    {
        MatrixDimension::Swap(other);
        AllocatorType::Swap(other);
    }

    const DataType& GetRef(int row, int col) const
    {
        CheckRange(row, col);
        return *AllocatorType::Data(row + GetLD() * col);
    }

    DataType& GetRef(int row, int col)
    {
        CheckRange(row, col);
        return *AllocatorType::Data(row + GetLD() * col);
    }

    const DataType* GetPtr(int row, int col) const
    {
        CheckRange(row, col);
        return AllocatorType::Data(row + GetLD() * col);
    }

    DataType* GetPtr(int row, int col)
    {
        CheckRange(row, col);
        return AllocatorType::Data(row + GetLD() * col);
    }

    Vector<SharedAllocator<DataType>> GetRow(int row, const Range& colRange) = delete;
    //Vector<SharedAllocator<DataType>> GetRow(int row, const Range& colRange)
    //{
    //    assert(false);
    //    return GetColumn(colRange, row);
    //}

    Vector<SharedAllocator<DataType>> GetColumn(const Range& rowRange, int col)
    {
        return Vector<SharedAllocator<typename AllocatorType::MemberType>>(GetPtr(rowRange.First(), col), rowRange.Count(RowSize() - 1));
    }
};


template<typename DataType>
class DiagonalMatrix
{
private:
    DataType m_Value;

public:
    DiagonalMatrix(DataType value = 1.0)
        : m_Value(value) {}

    DataType GetValue() const { return m_Value; }

    DiagonalMatrix& operator=(DataType value)
    {
        m_Value = value;
        return *this;
    }
};

template<typename DataType>
DiagonalMatrix<DataType> operator*(DataType value, const DiagonalMatrix<DataType>& diag)
{
    return DiagonalMatrix<DataType>(value * diag.GetValue());
}

//template<typename DataType>
//DiagonalMatrix<DataType> operator*(const DiagonalMatrix<DataType>& diag, DataType value)
//{
//    return DiagonalMatrix<DataType>(value * diag.GetValue());
//}


template<typename StorageType = CMatrixStorage<DynamicAllocator<double>>>
class Matrix : public StorageType
{
private:
    template<typename StorageType1>
    void CopyMat(const Matrix<StorageType1>& other);

public:
    using DataType = typename StorageType::DataType;
    using TNumType = typename StorageType::DataType;
};

template<typename StorageType>
template<typename StorageType1>
void Matrix<StorageType>::CopyMat(const Matrix<StorageType1>& other)
{
    // Valid for FMatrixStorage only
    if (other.GetLD() == StorageType::GetLD() && StorageType::GetLD() == StorageType::RowSize())
    {
        int size = StorageType::RowSize() * StorageType::ColSize();
        DataType* pDest = StorageType::Data();
        const DataType* pSrc = other.Data();
        for (int i = 0; i < size; i++)
        {
            pDest[i] = pSrc[i];
        }
    }
    else 
    {
        DataType* pDest = StorageType::Data();
        const DataType* pSrc = other.Data();
        int nRow = StorageType::RowSize();
        int nCol = StorageType::ColSize();

        for (int i = 0; i < nCol; i++) 
        {
            for (int j = 0; j < nRow; j++)
            {
                pDest[j] = pSrc[j];
            }
            pDest += StorageType::GetLD();
            pSrc += other.GetLD();
        }
    }
}



//==============//
// Type Defines //
//==============//

typedef DiagonalMatrix<double> RealDiagonalMatrix;