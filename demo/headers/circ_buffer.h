#include <vector>
#include <functional>
#include <utility>

// Ben - yet another circular buffer implementation

template <class Type>
class CircBuffer
{
public:
    CircBuffer(unsigned int maxSize, bool preAllocateMem = true):
        m_maxSize(maxSize)
    {
        if( preAllocateMem )
        {
            m_buffer.reserve( m_maxSize );
        }

        m_fill = 0;
        m_index = -1;
    }

    void add( Type const & item)
    {
        // Not thread-safe
        if ( m_fill < m_maxSize)
        {
            m_buffer.push_back( item );
            ++m_index;
            ++m_fill;
        }
        else
        {
            m_index = (m_index +1) % m_maxSize;
            m_buffer[m_index] = item;
        }
    }

    void add( Type && item)
    {
        // Not thread-safe
        if ( m_fill < m_maxSize)
        {
            ++m_index;
            m_buffer.emplace_back( item );
            ++m_fill;
        }
        else
        {
            m_index = (m_index +1) % m_maxSize;
            m_buffer[m_index] = std::move(item);;
        }
    }

    Type const & front() const
    {
        // Returns the most recent value
        return m_buffer[m_index];
    }

    std::vector<Type> const & buffer() const
    {
        // Get a (const) handle to the buffer, useful to compute standard vector statistics for instance
        return m_buffer;
    }

    void clear()
    {
        m_fill = m_index = 0;
        m_buffer.clear();
    }

    inline bool empty() const
    {
        return !valid();
    }

    Type const & back() const
    {
        // Returns the oldest value
        size_t const index = m_fill == m_maxSize ? (m_index + 1) % m_maxSize : 0;
        return m_buffer[index];
    }

    unsigned int maxSize() const
    {
        return m_maxSize;
    }

    unsigned int size() const
    {
        return m_fill;
    }

    inline bool valid() const
    {
        return m_fill > 0;
    }

    Type const & operator[](int nIndex) const
    {
        unsigned int const realIndex = ( m_index + nIndex) % m_maxSize;
        return m_buffer[realIndex];
    }

private:
    std::vector<Type> m_buffer;
    unsigned int const m_maxSize;
    unsigned int     m_fill;
    unsigned int     m_index;
};
