#pragma once

#ifdef _DEBUG
  #ifndef _CHECK_RANGE
    #define _CHECK_RANGE
  #endif
#endif


template <typename DataType>
class DynamicAllocator
{
private:
    DataType* m_Data;

protected:
    DynamicAllocator()
        : m_Data(nullptr) {}

    DynamicAllocator(int size) 
    { 
        Allocate(size); 
    }

    //DynamicAllocator(DataType* data, int ID) 
    //{ 
    //    assert(false); 
    //}
    DynamicAllocator(DataType* data, int ID) = delete;

    void Swap(DynamicAllocator& other)
    {
        std::swap(m_Data, other.m_Data);
    }

    void Allocate(int size) 
    { 
        m_Data = new DataType[size]; 
    }

    void DeAllocate()
    {
        if (m_Data)
        {
            delete m_Data;
            m_Data = nullptr;
        }
    }

    // Why is this function needed?
    //int GetID(int ID) const { return ID; }
    
    template<typename MatrixType>
    void InitMatCopy(const MatrixType& src, MatrixType& dest)
    {
        dest = src;
    }

    enum { DEFAULT_SIZE = 0, DEFAULT_ROWS = 0, DEFAULT_COLUMNS = 0 };

public:
    //typedef DataType MemberType;
    using MemberType = DataType;

    DataType* Data(int idx = 0) 
    { 
        return m_Data + idx; 
    }

    const DataType* Data(int idx = 0) const 
    {
        return m_Data + idx; 
    }
};


template <typename DataType, int N = 3>
class StaticAllocator
{
private:
    DataType m_Data[N];

protected:
    StaticAllocator() {}

    StaticAllocator(int size) 
    { 
        Allocate(size); 
    }

    // Why is this constructor needed?
    //StaticAllocator(DataType* data, int ID)
    //{
    //    assert(false);
    //}
    StaticAllocator(DataType* data, int ID) = delete;

    void Allocate(int size) 
    { 
        assert(size <= N); 
    }

    void DeAllocate() {}

    // Why is this function needed?
    //int GetID(int ID) const { return ID; }

    template <typename MatrixType>
    void InitMatCopy(const MatrixType& src, MatrixType& dest)
    {
        dest = src;
    }

    enum { DEFAULT_SIZE = N };

public:
    //typedef DataType MemberType;
    using MemberType = DataType;

    DataType* Data(int idx = 0) 
    { 
        return m_Data + idx; 
    }

    const DataType* Data(int idx = 0) const 
    { 
        return m_Data + idx; 
    }
};


template <typename DataType>
class SharedAllocator
{
private:
    DataType* m_Data;

protected:
    /*SharedAllocator() 
    { 
        assert(false); 
    }*/
    SharedAllocator() = delete;

    //SharedAllocator(int size) 
    //{ 
    //    assert(false); 
    //}
    SharedAllocator(int size) = delete;

    SharedAllocator(DataType* data)
        : m_Data(data) {}

    void Swap(SharedAllocator& other)
    { 
        std::swap(m_Data, other.m_Data); 
    }

    //void Allocate(int size) 
    //{ 
    //    assert(false); 
    //}
    void Allocate(int size) = delete;

    void DeAllocate() {}

    template <typename MatrixType>
    void InitMatCopy(const MatrixType& src, MatrixType& dest) 
    {
        m_Data = src.m_Data;
    }

    enum { DEFAULT_SIZE = 0 };

public:
    //typedef DataType MemberType;
    using MemberType = DataType;

    DataType* Data(int idx = 0)
    {
        return m_Data + idx;
    }

    const DataType* Data(int idx = 0) const
    {
        return m_Data + idx;
    }
};