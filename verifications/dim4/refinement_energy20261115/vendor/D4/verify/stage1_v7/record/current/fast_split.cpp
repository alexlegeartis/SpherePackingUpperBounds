#include <cstdint>
#include <cstddef>
#include <cstring>

// Exact midpoint de Casteljau subdivision of integer Bernstein bounds.
// The caller checks |input| < 2^59. Every intermediate is an average
// rounded toward the selected endpoint, so this invariant is preserved.
// The addition has absolute value below 2^60 and cannot overflow int64.
static inline std::int64_t floor_half(std::int64_t value) {
    // C++20 specifies arithmetic signed right shift, rounding toward
    // negative infinity. This is Python's //2 for both signs.
    return value >> 1;
}

extern "C" int split_bernstein_i64(std::int64_t* work,
        std::int64_t* left, std::int64_t* right,
        std::size_t rows, std::size_t columns, int upper) {
    if (rows == 0 || columns == 0 || (upper != 0 && upper != 1)) return 1;
    const std::size_t n = rows - 1;
    std::memcpy(left, work, columns * sizeof(std::int64_t));
    std::memcpy(right + n*columns, work + n*columns, columns * sizeof(std::int64_t));
    for (std::size_t level = 1; level <= n; ++level) {
        const std::size_t length = (rows - level)*columns;
        for (std::size_t k = 0; k < length; ++k)
            work[k] = floor_half(work[k] + work[k+columns] + upper);
        std::memcpy(left + level*columns, work, columns * sizeof(std::int64_t));
        std::memcpy(right + (n-level)*columns, work + (n-level)*columns,
                    columns * sizeof(std::int64_t));
    }
    return 0;
}
