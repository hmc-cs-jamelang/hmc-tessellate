#pragma once

namespace wrapper {
    namespace detail {
        struct DefaultTag {};
    }

    template <typename WrappedType, typename Tag = detail::DefaultTag>
    class Wrapper {
    protected:
        WrappedType value;

    public:
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
            return (stream << wrapper.value);
        }
    };
}
