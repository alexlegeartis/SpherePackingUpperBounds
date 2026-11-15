#include <cstdint>
#include <cstddef>

// Optional sufficient inequality only.  No branch, sign, or coverage state is
// altered here.  The wrapper validates |high|<2^59, 0<=q<=2^20, 0<=T<2^59.
// Both products are strictly below 2^119 and fit in unsigned __int128.
extern "C" int factor_correlated_i64(const std::int64_t* high,
        const std::uint64_t* q0, const std::uint64_t* q1,
        const std::uint64_t* q2, std::size_t n0, std::size_t n1,
        std::size_t n2, std::uint64_t T) {
    if (!n0 || !n1 || !n2 || T >= (std::uint64_t(1)<<59)) return -1;
    using U128 = unsigned __int128;
    std::size_t at=0;
    for (std::size_t i=0; i<n0; ++i) {
        if (q0[i] > (std::uint64_t(1)<<20)) return -1;
        for (std::size_t j=0; j<n1; ++j) {
            if (q1[j] > (std::uint64_t(1)<<20)) return -1;
            const U128 left_factor=U128(T)*q0[i]*q1[j];
            for (std::size_t k=0; k<n2; ++k,++at) {
                if (q2[k] > (std::uint64_t(1)<<20)) return -1;
                if (high[at]<=0) continue;
                if (high[at] >= (std::int64_t(1)<<59)) return -1;
                const U128 lhs=U128(high[at])<<60;
                const U128 rhs=left_factor*q2[k];
                if (lhs>rhs) return 0;
            }
        }
    }
    return 1;
}
