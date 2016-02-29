#pragma once
#include <utility>
#include <assert.h>

namespace verification {
    template <typename F>
    struct ScopeExit {
        F f;
        ScopeExit(F f) : f(std::move(f)) {}
        ScopeExit(ScopeExit&& rhs) : f(std::move(rhs.f)) {}

        ScopeExit() = delete;
        ScopeExit(const ScopeExit&) = delete;
        ScopeExit& operator=(const ScopeExit&) = delete;

        ~ScopeExit() {f();}
    };

    template <typename G>
    ScopeExit<G> makeScopeExit(G g)
    {
        return {g};
    }
}

#define VERIFY_STRING_JOIN(a, b) VERIFY_DO_STRING_JOIN(a, b)
#define VERIFY_DO_STRING_JOIN(a, b) a ## b

#define VERIFICATION_ANONYMOUS_VARIABLE \
VERIFY_STRING_JOIN(verify_scope_exit_, __LINE__)



#define VERIFICATION_INVARIANT(...) \
VERIFICATION(__VA_ARGS__) \
VERIFICATION_EXIT(__VA_ARGS__)

#define VERIFY_INVARIANT(...) \
VERIFY(__VA_ARGS__) \
VERIFY_EXIT(__VA_ARGS__)

#define V_SCOPE_EXIT(...) \
auto VERIFICATION_ANONYMOUS_VARIABLE \
{verification::makeScopeExit([&]{__VA_ARGS__})};




#if !defined(NDEBUG) && !defined(NVERIFY)

    #define VERIFICATION(...) __VA_ARGS__

    #define NO_VERIFICATION(...)

    #define VERIFICATION_EXIT(...) V_SCOPE_EXIT(__VA_ARGS__)

    #define VERIFY(...) VERIFICATION(assert(__VA_ARGS__);)

    #define VERIFY_EXIT(...) VERIFICATION_EXIT(assert(__VA_ARGS__);)

#else

    #define VERIFICATION(...)

    #define NO_VERIFICATION(...) __VA_ARGS__

    #define VERIFICATION_EXIT(...)

    #define VERIFY(...)

    #define VERIFY_EXIT(...)

#endif
