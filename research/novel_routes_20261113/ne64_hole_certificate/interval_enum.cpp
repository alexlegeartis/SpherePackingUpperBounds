// Rigorous outward-rounded enumeration of a congruence class modulo 3.
// Compile: g++ -O3 -std=c++17 -fno-fast-math -ffp-contract=off interval_enum.cpp -o interval_enum
// Usage: interval_enum ldl_bounds.txt gram.txt coset.txt shard_id shard_count
// All shard IDs 0,...,shard_count-1 must report COMPLETE and found=false.
// The parent minimum is an external theorem; this enumerates only a dual coset.
#include <algorithm>
#include <cfenv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

using namespace std;
static constexpr int N = 64;
static constexpr int SPLIT_DEPTH = 47;  // 16 upper coordinates assigned
static constexpr long long COORD_LIMIT = 1000000;
static constexpr double BOUND = 60.;
static constexpr double INF = numeric_limits<double>::infinity();
static_assert(numeric_limits<double>::is_iec559 && numeric_limits<double>::digits == 53,
              "This verifier requires IEEE-754 binary64.");

long long mu_lo[N][N], mu_hi[N][N], Q[N][N], z[N], residue[N], witness[N];
double d_lo[N];
uint64_t nodes = 0, leaves = 0, prefixes = 0, assigned_prefixes = 0;
int shard_id, shard_count;
bool found = false;
long long witness_norm = -1;
chrono::steady_clock::time_point start;

inline double down(double x) { return nextafter(x, -INF); }
inline double up(double x) { return nextafter(x, INF); }
[[noreturn]] void fail(const char* reason) { cerr << "ABORT " << reason << endl; exit(3); }

void walk(int k, double lower_partial) {
    if (found) return;
    ++nodes;
    if (k == SPLIT_DEPTH) {
        uint64_t index = prefixes++;
        if (index % (uint64_t)shard_count != (uint64_t)shard_id) return;
        ++assigned_prefixes;
    }
    if (nodes % 100000000 == 0) {
        cerr << "nodes " << nodes << " elapsed "
             << chrono::duration<double>(chrono::steady_clock::now() - start).count()
             << " depth " << k << endl;
    }
    if (k < 0) {
        ++leaves;
        __int128 norm = 0;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                norm += (__int128)z[i] * Q[i][j] * z[j];
        if (norm <= 60) {
            if (norm < 0) fail("negative norm");
            witness_norm = (long long)norm;
            copy(z, z + N, witness);
            found = true;
        }
        return;
    }

    // Exact integer sums of dyadic endpoints. With |mu endpoint| < 2^53,
    // |z_j| <= 10^6, and at most 64 terms, the numerator has size < 2^79,
    // safely inside signed __int128. The double conversion then one nextafter
    // encloses the exact value; scaling by 2^-40 is exact in this range.
    __int128 sum_lo = 0, sum_hi = 0;
    for (int j = k + 1; j < N; ++j) {
        if (z[j] >= 0) {
            sum_lo += (__int128)mu_lo[k][j] * z[j];
            sum_hi += (__int128)mu_hi[k][j] * z[j];
        } else {
            sum_lo += (__int128)mu_hi[k][j] * z[j];
            sum_hi += (__int128)mu_lo[k][j] * z[j];
        }
    }
    double center_lo = down(ldexp((double)sum_lo, -40));
    double center_hi = up(ldexp((double)sum_hi, -40));
    if (lower_partial > BOUND) return;
    double radius = up(sqrt(up(up(BOUND - lower_partial) / d_lo[k])));
    if (!isfinite(radius) || radius > COORD_LIMIT ||
        fabs(center_lo) > COORD_LIMIT || fabs(center_hi) > COORD_LIMIT)
        fail("unsafe interval radius");
    long long lo = (long long)ceil(down(-center_hi - radius));
    long long hi = (long long)floor(up(-center_lo + radius));
    if (llabs(lo) > COORD_LIMIT || llabs(hi) > COORD_LIMIT)
        fail("unsafe coordinate");
    long long shift = ((residue[k] - lo) % 3 + 3) % 3;
    for (long long x = lo + shift; x <= hi; x += 3) {
        double a = down((double)x + center_lo);
        double b = up((double)x + center_hi);
        double near = a > 0 ? a : (b < 0 ? -b : 0.);
        double term = down(d_lo[k] * down(near * near));
        double next_lower = down(lower_partial + max(0., term));
        if (next_lower <= BOUND) {
            z[k] = x;
            walk(k - 1, next_lower);
            if (found) return;
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 6) return 2;
    if (fesetround(FE_TONEAREST) != 0) fail("rounding mode");
    shard_id = stoi(argv[4]); shard_count = stoi(argv[5]);
    if (shard_count < 1 || shard_id < 0 || shard_id >= shard_count) return 2;
    ifstream ldl(argv[1]), gram(argv[2]), coset(argv[3]);
    for (int i = 0; i < N; ++i) {
        long long lo, hi;
        ldl >> lo >> hi;
        if (!ldl || lo <= 0 || hi < lo || llabs(lo) >= (1LL << 53) || llabs(hi) >= (1LL << 53))
            fail("invalid diagonal enclosure");
        d_lo[i] = ldexp((double)lo, -40);
        for (int j = 0; j < N; ++j) {
            ldl >> lo >> hi;
            if (!ldl || hi < lo || llabs(lo) >= (1LL << 53) || llabs(hi) >= (1LL << 53))
                fail("invalid coefficient enclosure");
            mu_lo[i][j] = lo; mu_hi[i][j] = hi;
        }
    }
    for (auto& row : Q) for (auto& x : row) {
        gram >> x;
        if (!gram || llabs(x) > 1000000) fail("invalid gram input");
    }
    for (auto& x : residue) {
        coset >> x;
        if (!coset || x < 0 || x > 2) fail("invalid residue input");
    }
    start = chrono::steady_clock::now();
    walk(N - 1, 0.);
    double seconds = chrono::duration<double>(chrono::steady_clock::now() - start).count();
    cout << "{\"status\":\"" << (found ? "EXACT_REJECT" : "COMPLETE")
         << "\",\"found\":" << (found ? "true" : "false")
         << ",\"shard_id\":" << shard_id << ",\"shard_count\":" << shard_count
         << ",\"split_depth\":" << SPLIT_DEPTH << ",\"nodes\":" << nodes
         << ",\"leaves\":" << leaves << ",\"prefixes\":" << prefixes
         << ",\"assigned_prefixes\":" << assigned_prefixes
         << ",\"seconds\":" << seconds << ",\"norm\":" << witness_norm
         << ",\"witness\":[";
    if (found) for (int i = 0; i < N; ++i) cout << (i ? "," : "") << witness[i];
    cout << "]}" << endl;
    return found ? 1 : 0;
}
