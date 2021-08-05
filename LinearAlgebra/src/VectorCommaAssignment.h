#pragma once


template<typename AllocatorType>
class Vector;

template<typename AllocatorType>
class VectorCommaAssignment
{
private:
    //typedef typename AllocatorType::MemberType DataType;
    using DataType = typename AllocatorType::MemberType;

    Vector<AllocatorType>* m_pVec;
    int m_IP;
    DataType m_Value;

public:
    VectorCommaAssignment(Vector<AllocatorType>* pVec, int ip, const DataType& value)
        : m_pVec(pVec), m_IP(ip), m_Value(value) {}

    ~VectorCommaAssignment()
    {
        if (m_pVec)
        {
            for (int i = m_IP; i < m_pVec->Size(); i++)
            {
                (*m_pVec)(i) = m_Value;
            }
        }
    }

    operator DataType() { return m_Value; }

    VectorCommaAssignment operator,(DataType value)
    {
        Vector<AllocatorType>* pVec = m_pVec;
        if (m_IP < m_pVec->Size())
        {
            (*m_pVec)(m_IP) = value;
            m_pVec = nullptr;
        }

        return VectorCommaAssignment(pVec, m_IP + 1, value);
    }
};
