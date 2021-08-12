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

    MatrixDimension(int nRow, int nCol)
        : m_Row(nRow), m_Col(nCol) {}

    void Swap(MatrixDimension& matDim)
    {
        std::swap(m_Row, matDim.m_Row);
        std::swap(m_Col, matDim.m_Row);
    }

    bool Resized(int nRow, int nCol)
    {
        bool flag = false;
        
        if (m_Row != nRow)
        {
            m_Row = nRow;
            flag = true;
        }

        if (m_Col != nCol)
        {
            m_Col = nCol;
            flag = true;
        }

        return flag;
    }

    void CheckRange(int nRow, int nCol) const
    {
        #ifdef _CHECK_RANGE
          assert(m_Row > nRow);
          assert(m_Col > nCol);
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
    Shared2DAllocator()
        : SharedAllocator() {}

    Shared2DAllocator(int size)
        : SharedAllocator(size) {}

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

    Static2DAllocator(int size)
        : StaticAllocator(size) {}

    Static2DAllocator(DataType* data, int ld)
        : StaticAllocator(data) {}

    //void Swap(Static2DAllocator& other)
    //{
    //    assert(false);
    //}
    void Swap(Static2DAllocator& other) = delete;

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

    CMatrixStorage(int nRow, int nCol)
        : MatrixDimension(nRow, nCol), AllocatorType(Size()) {}

    CMatrixStorage(DataType* data, int nRow, int nCol, int ld)
        : MatrixDimension(nRow, nCol), AllocatorType(data, ld) {}

    void Swap(CMatrixStorage& other)
    {
        MatrixDimension::Swap(other);
        AllocatorType::Swap(other);
    }

    const DataType& GetRef(int iRow, int iCol) const
    {
        CheckRange(iRow, iCol);
        return *AllocatorType::Data(iRow * GetLD() + iCol);
    }

    DataType& GetRef(int iRow, int iCol)
    {
        CheckRange(iRow, iCol);
        return *AllocatorType::Data(iRow * GetLD() + iCol);
    }

    const DataType* GetPtr(int iRow, int iCol) const
    {
        CheckRange(iRow, iCol);
        return AllocatorType::Data(iRow * GetLD() + iCol);
    }

    DataType* GetPtr(int iRow, int iCol)
    {
        CheckRange(iRow, iCol);
        return AllocatorType::Data(iRow * GetLD() + iCol);
    }

    Vector<SharedAllocator<DataType>> GetColumn(const Range& nRowRange, int iCol) = delete;
    //Vector<SharedAllocator<DataType>> GetColumn(const Range& nRowRange, int iCol)
    //{
    //    assert(false);
    //    return GetRow(iCol, nRowRange);
    //}

    // Get Row Vector with Specific Range
    // GetRow(0, Range::All()): 0-th row
    // GetRow(1, Range(0, 1)): 1-th row from iColumn 0 to column 1
    Vector<SharedAllocator<DataType>> GetRow(int iRow, const Range& colRange)
    {
        return Vector<SharedAllocator<typename AllocatorType::MemberType>>(GetPtr(iRow, colRange.First()), colRange.Count(ColSize() - 1));
    }

public:
    const DataType* operator[](int iRow) const { return GetPtr(iRow, 0); }
    DataType* operator[](int iRow) { return GetPtr(iRow, 0); }
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

    FMatrixStorage(int nRow, int nCol)
        : MatrixDimension(nRow, nCol), AllocatorType(Size()) {}

    FMatrixStorage(DataType* data, int nRow, int nCol, int ld)
        : MatrixDimension(nRow, nCol), AllocatorType(data, ld) {}

    void Swap(FMatrixStorage& other)
    {
        MatrixDimension::Swap(other);
        AllocatorType::Swap(other);
    }

    const DataType& GetRef(int iRow, int iCol) const
    {
        CheckRange(iRow, iCol);
        return *AllocatorType::Data(iRow + GetLD() * iCol);
    }

    DataType& GetRef(int iRow, int iCol)
    {
        CheckRange(iRow, iCol);
        return *AllocatorType::Data(iRow + GetLD() * iCol);
    }

    const DataType* GetPtr(int iRow, int iCol) const
    {
        CheckRange(iRow, iCol);
        return AllocatorType::Data(iRow + GetLD() * iCol);
    }

    DataType* GetPtr(int iRow, int iCol)
    {
        CheckRange(iRow, iCol);
        return AllocatorType::Data(iRow + GetLD() * iCol);
    }

    Vector<SharedAllocator<DataType>> GetRow(int iRow, const Range& colRange) = delete;
    //Vector<SharedAllocator<DataType>> GetRow(int iRow, const Range& colRange)
    //{
    //    assert(false);
    //    return GetColumn(colRange, iRow);
    //}

    Vector<SharedAllocator<DataType>> GetColumn(const Range& rowRange, int iCol)
    {
        return Vector<SharedAllocator<typename AllocatorType::MemberType>>(GetPtr(rowRange.First(), iCol), rowRange.Count(RowSize() - 1));
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
    using Allocator2DType = typename StorageType::Allocator2DType;
    using SharedMatrixType = typename StorageType::SharedMatrixType;

    Matrix()
        : StorageType() {}

    // Constructor for Dynamic and Static Allocator
    Matrix(int nRow, int nCol)
        : StorageType(nRow, nCol) {}

    // Constructor for Shared2DAllocator Only
    Matrix(const DataType* data, int nRow, int nCol, int ld)
        : StorageType(const_cast<DataType*>(data), nRow, nCol, ld) {}

    // Copy Constructor
    Matrix(const Matrix& mat)
        : StorageType(mat.RowSize(), mat.ColSize())
    {
        Allocator2DType::Allocate(mat.Size());
        Allocator2DType::InitMatCopy(mat, *this);
    }

    template<typename AllocatorType>
    Matrix(Vector<AllocatorType>& vec)
        : StorageType(vec.Data(), vec.Size(), 1, vec.Size()) {}

    // Destructor
    ~Matrix() { Allocator2DType::DeAllocate(); }

    void Resize(int nRow, int nCol)
    {
        if (MatrixDimension::Resized(nRow, nCol))
        {
            Allocator2DType::DeAllocate();
            Allocator2DType::Allocate(MatrixDimension::Size());
        }
    }

    void Reserve(int nRow, int nCol)
    {
        //if (MatrixDimension::RowSize() < nRow || nCol < nCol)
        if (nRow * nCol > MatrixDimension::Size())
        {
            Resize(nRow, nCol);
        }
    }

    void Swap(Matrix<StorageType>& mat)
    {
        StorageType::Swap(mat);
    }

    const DataType& operator()(int iRow, int iCol) const
    {
        return StorageType::GetRef(iRow, iCol);
    }

    DataType& operator()(int iRow, int iCol)
    {
        return StorageType::GetRef(iRow, iCol);
    }

    SharedMatrixType operator()(const Range& rowRange, const Range& colRange)
    {
        int nRow = rowRange.Count(MatrixDimension::RowSize() - 1);
        int nCol = colRange.Count(MatrixDimension::ColSize() - 1);

        return SharedMatrixType(const_cast<DataType*>(StorageType::GetPtr(rowRange.First(), colRange.First())), 
            nRow, nCol, StorageType::GetLD());
    }

    Vector<SharedAllocator<DataType>> operator()(const Range& rowRange, int iCol)
    {
        return StorageType::GetColumn(rowRange, iCol);
    }

    Vector<SharedAllocator<DataType>> operator()(int iRow, const Range& colRange)
    {
        return StorageType::GetRow(iRow, colRange);
    }

    bool IsIdentity() const
    {
        int nRow = MatrixDimension::RowSize();
        int nCol = MatrixDimension::ColSize();

        if (nRow != nCol)
        {
            return false;
        }

        for (int i = 0; i < nRow; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (!IsZero(StorageType::GetRef(i, j))
                {
                    return false;
                }
            }

            if (!IsEqual(StorageType::GetRef(i, j), 1.0))
            {
                return false;
            }

            for (int j = i + 1; j < n; j++)
            {
                if (!IsZero(StorageType::GetRef(i, j)))
                {
                    return false;
                }
            }
        }

        return true;
    }

    bool operator==(const DataType& value) const
    {
        int nRow = MatrixDimension::RowSize();
        int nCol = MatrixDimension::ColSize();

        for (int i = 0; i < nCol; i++)
        {
            for (int j = 0; j < nRow; j++)
            {
                if (!IsZero(StorageType::GetRef(j, i) - value))
                {
                    return false;
                }
            }
        }

        return true;
    }
};

template<typename StorageType>
template<typename StorageType1>
void Matrix<StorageType>::CopyMat(const Matrix<StorageType1>& other)
{
    // Valid for FMatrixStorage only
    if (other.GetLD() == StorageType::GetLD() && StorageType::GetLD() == StorageType::RowSize())
    {
        int size = StorageType::RowSize() * StorageType::ColSize();
        const DataType* pSrc = other.Data();
        DataType* pDest = StorageType::Data();
        
        for (int i = 0; i < size; i++)
        {
            pDest[i] = pSrc[i];
        }
    }
    else 
    {
        const DataType* pSrc = other.Data();
        DataType* pDest = StorageType::Data();        
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