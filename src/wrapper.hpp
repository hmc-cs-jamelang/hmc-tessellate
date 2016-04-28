#pragma once

#include <utility>

namespace hmc {
    namespace detail {
        struct WrapperDefaultTag {};
    }

    template <typename WrappedType, typename Tag = detail::WrapperDefaultTag>
    struct Wrapper {
        WrappedType value;

        constexpr Wrapper() = default;
        explicit constexpr Wrapper(WrappedType value) : value(value) {}
        explicit constexpr operator WrappedType() const {return value;}

        friend constexpr bool operator==(const Wrapper lhs, const Wrapper rhs)
        {
            return lhs.value == rhs.value;
        }

        friend constexpr bool operator!=(const Wrapper lhs, const Wrapper rhs)
        {
            return lhs.value != rhs.value;
        }

        friend std::ostream& operator<<(std::ostream& stream, const Wrapper& wrapper)
        {
            return stream << wrapper.value;
        }

        friend void swap(Wrapper& a, Wrapper& b)
        {
            using std::swap;
            swap(a.value, b.value);
        }

        friend struct std::hash<Wrapper>;
    };
}

namespace std {
    template <typename WrappedType, typename Tag>
    struct hash<hmc::Wrapper<WrappedType, Tag>> {
        std::size_t operator()(const hmc::Wrapper<WrappedType, Tag>& w) const {
            return std::hash<WrappedType>()(w.value);
        }
    };
}
