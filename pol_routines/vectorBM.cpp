#include <vector>
#include <cstdint>
#include <chrono>
#include <algorithm>

using LONG = int64_t;

// Prevents the compiler from getting too aggressive.
static volatile LONG bench_sink = 0;

extern "C" int benchmarkVectorCreations_C(
    int mLen,
    int degM,
    const LONG *M,
    int uLen,
    int degU,
    const LONG *U,
    int repeats,
    long long *time_m_ns,
    long long *time_u_ns,
    long long *time_r_ns,
    long long *time_t_ns
)
{
    using clock = std::chrono::steady_clock;

    if (!M || !U || !time_m_ns || !time_u_ns || !time_r_ns || !time_t_ns) {
        return -1;
    }
    if (degM < 0 || degU < 0 || repeats <= 0) {
        return -1;
    }
    if (mLen <= 0 || uLen <= 0) {
        return -1;
    }
    if (degM >= mLen || degU >= uLen) {
        return -1;
    }

    *time_m_ns = 0;
    *time_u_ns = 0;
    *time_r_ns = 0;
    *time_t_ns = 0;

    const int wsSize = std::max(degM, degU) + 1;

    for (int k = 0; k < repeats; ++k) {
        auto t0 = clock::now();
        std::vector<LONG> m(M, M + (degM + 1));
        auto t1 = clock::now();
        *time_m_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        t0 = clock::now();
        std::vector<LONG> u(U, U + (degU + 1));
        t1 = clock::now();
        *time_u_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        t0 = clock::now();
        std::vector<LONG> rTmp(wsSize, 0);
        t1 = clock::now();
        *time_r_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        t0 = clock::now();
        std::vector<LONG> tTmp(wsSize, 0);
        t1 = clock::now();
        *time_t_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

        // Touch the data so the compiler cannot treat these as useless.
        bench_sink += m[0] + u[0] + rTmp[0] + tTmp[0];
        bench_sink += static_cast<LONG>(m.size() + u.size() + rTmp.size() + tTmp.size());
    }

    return 0;
}