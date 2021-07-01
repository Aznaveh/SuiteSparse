/** =========================================================================  /
 *  ======================  paru_allocator ==================================  /
 *  ========================================================================= */
/*! @brief   Defining an allocator for C++ data structures so that I have just
 *           one allocator through all the package (paru_alloc)
 *  @author Aznaveh
 */
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>

#include "paru_internal.hpp"

template <class T>
struct paru_allocator
{
    typedef T value_type;

    paru_allocator() = default;
    template <class U>
    constexpr paru_allocator(const paru_allocator<U>&) noexcept {}

        [[nodiscard]] T* allocate(std::size_t n)
    {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_array_new_length();

        if (auto p = static_cast<T*>(paru_alloc(n, sizeof(T))))
        {
            report(p, n);
            return p;
        }

        throw std::bad_alloc();
    }

    void deallocate(T* p, std::size_t n) noexcept
    {
        report(p, n, 0);
        paru_free(n, sizeof(T), p);
    }

private:
    void report(T* p, std::size_t n, bool alloc = true) const
    {
#ifndef NDEBUG
        std::cout << (alloc ? "Alloc: " : "Dealloc: ") << sizeof(T) * n
                  << " bytes at " << std::hex << std::showbase
                  << reinterpret_cast<void*>(p) << std::dec << '\n';
#endif
    }
};

template <class T, class U>
bool operator==(const paru_allocator<T>&, const paru_allocator<U>&)
{
    return true;
}
template <class T, class U>
bool operator!=(const paru_allocator<T>&, const paru_allocator<U>&)
{
    return false;
}

