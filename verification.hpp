#pragma once
#include <utility>

namespace Verification {
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
    ScopeExit<G> MakeScopeExit(G g)
    {
        return {g};
    }
}

#define VERIFY_STRING_JOIN(a, b) VERIFY_DO_STRING_JOIN(a, b)
#define VERIFY_DO_STRING_JOIN(a, b) a ## b

#define VERIFICATION_ANONYMOUS_VARIABLE \
VERIFY_STRING_JOIN(verify_scope_exit_, __LINE__)

#if !defined(NDEBUG) && !defined(NVERIFY)
#define VERIFICATION(...) __VA_ARGS__
#define NO_VERIFICATION(...)
#define VERIFICATION_EXIT(...) \
auto VERIFICATION_ANONYMOUS_VARIABLE {Verification::MakeScopeExit([&]{__VA_ARGS__})};
#else
#define VERIFICATION(...)
#define NO_VERIFICATION(...) __VA_ARGS__
#define VERIFICATION_EXIT(...)
#endif

