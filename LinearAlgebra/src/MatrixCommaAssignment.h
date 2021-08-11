#pragma once


template<typename StorageType>
class Matrix;

// For assigning values of a matrix with comma(,)
// Column Major Assignment Only -> Need to be modified for Row Major Assignment!!!
//
// --- Example ---
// 
// 
// ---------------


template<typename StorageType>
class MatrixCommaAssignment
{
private:
    using DataType = typename StorageType::DataType;

    Matrix<StorageType>* m_pMat;
    int m_IP;
    DataType m_Value;

    void GetCrd(int& iRow, int& iCol)
    {
        iRow = m_IP % m_pMat->RowSize();
        iCol = m_IP / m_pMat->RowSize();
    }

public:
    MatrixCommaAssignment(Matrix<StorageType>* pMat, int ip, const DataType& value)
        : m_pMat(pMat), m_IP(ip), m_Value(value) {}

    ~MatrixCommaAssignment()
    {
        if (m_pMat && m_IP < m_pMat->Size())
        {
            int iRow, iCol;
            GetCrd(iRow, iCol);

            // Valid for FMatrixStorage only
            DataType* pd = m_pMat->Data();
            int ld = m_pMat->GetLD();
            int nRow = m_pMat->RowSize();
            int nCol = m_pMat->ColSize();

            pd += iCol * ld;

            for (int i = iRow; i < nRow; i++)
            {
                pd[i] = m_Value;
            }

            for (int j = iCol + 1; j < nCol; j++)
            {
                pd += ld;
                for (int i = 0; i < nRow; i++)
                {
                    pd[i] = m_Value;
                }
            }
        }
    }

    operator DataType() { return m_Value; }

    MatrixCommaAssignment operator,(const DataType& value)
    {
        Matrix<StorageType>* pMat = m_pMat;
        if (m_IP < m_pMat->Size())
        {
            int iRow, iCol;
            GetCrd(iRow, iCol);

            (*m_pMat)(iRow, iCol) = value;
            m_pMat = nullptr;
        }

        return MatrixCommaAssignment(pMat, m_IP + 1, value);
    }
};