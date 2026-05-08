// nfsim_rng.h — Per-instance Mersenne Twister RNG for thread-safe NFsim
//
// Replaces the global static RNG (MTRand_int32 with static state[624]) with a
// per-instance RNG that can be owned by each System. This enables:
//   1. Thread safety — no shared mutable state between Systems
//   2. Deterministic reproducibility — same seed → same trajectory, always
//   3. Deterministic parallel behavior — independent of thread scheduling
//
// The MT19937 algorithm is identical to the original NFsim implementation.
// Only the storage has changed from static class members to instance members.

#ifndef NFSIM_RNG_H_
#define NFSIM_RNG_H_

#include <cmath>
#include <ctime>

namespace NFcore {

class NfsimRNG {
public:
    NfsimRNG() : p_(n_), init_(false) {
        // Default: seed with 5489 (same as original MTRand default)
        seed(5489UL);
    }

    explicit NfsimRNG(unsigned long s) : p_(n_), init_(false) {
        seed(s);
    }

    // ─── Seeding ─────────────────────────────────────────────────────────
    void seed(unsigned long s) {
        state_[0] = s & 0xFFFFFFFFUL;
        for (int i = 1; i < n_; ++i) {
            state_[i] = 1812433253UL * (state_[i - 1] ^ (state_[i - 1] >> 30)) + i;
            state_[i] &= 0xFFFFFFFFUL;
        }
        p_ = n_;  // force gen_state() on next call
        init_ = true;
        haveNextGaussian_ = false;
    }

    // ─── Core: 32-bit random integer ─────────────────────────────────────
    unsigned long rand_int32() {
        if (p_ == n_) gen_state();
        unsigned long x = state_[p_++];
        x ^= (x >> 11);
        x ^= (x << 7) & 0x9D2C5680UL;
        x ^= (x << 15) & 0xEFC60000UL;
        return x ^ (x >> 18);
    }

    // ─── Uniform double [0, 1) ───────────────────────────────────────────
    double rand_half_open() {
        return static_cast<double>(rand_int32()) * (1.0 / 4294967296.0);
    }

    // ─── Uniform double [0, 1] ───────────────────────────────────────────
    double rand_closed() {
        return static_cast<double>(rand_int32()) * (1.0 / 4294967295.0);
    }

    // ─── Uniform double (0, 1) ───────────────────────────────────────────
    double rand_open() {
        return (static_cast<double>(rand_int32()) + 0.5) * (1.0 / 4294967296.0);
    }

    // ─── NFutil-compatible API (same semantics as NFutil::RANDOM etc.) ───

    // Uniform on (0, max] — used for reaction selection
    double random(double max) {
        return (1.0 - rand_half_open()) * max;
    }

    // Uniform on (0, 1) — used for dt calculation
    double random_open() {
        return rand_open();
    }

    // Uniform on [0, 1] — used for dt calculation (alternative)
    double random_closed() {
        return rand_closed();
    }

    // Uniform integer on [min, max) — used for molecule selection
    int random_int(unsigned long min, unsigned long max) {
        return static_cast<int>(min + static_cast<unsigned long>((max - min) * rand_half_open()));
    }

    // Standard normal (mean=0, var=1) — Box-Muller polar method
    double random_gaussian() {
        if (haveNextGaussian_) {
            haveNextGaussian_ = false;
            return nextGaussian_;
        }
        double v1, v2, s;
        do {
            v1 = 2.0 * rand_open() - 1.0;
            v2 = 2.0 * rand_open() - 1.0;
            s = v1 * v1 + v2 * v2;
        } while (s >= 1.0 || s == 0.0);

        double multiplier = std::sqrt(-2.0 * std::log(s) / s);
        nextGaussian_ = v2 * multiplier;
        haveNextGaussian_ = true;
        return v1 * multiplier;
    }

private:
    static const int n_ = 624;
    static const int m_ = 397;

    unsigned long state_[624];  // per-instance (was static in MTRand_int32)
    int p_;                      // per-instance (was static)
    bool init_;                  // per-instance (was static)

    // Gaussian cache (was file-static in random.cpp)
    bool haveNextGaussian_ = false;
    double nextGaussian_ = 0.0;

    unsigned long twiddle(unsigned long u, unsigned long v) {
        return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1)
            ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
    }

    void gen_state() {
        for (int i = 0; i < (n_ - m_); ++i)
            state_[i] = state_[i + m_] ^ twiddle(state_[i], state_[i + 1]);
        for (int i = n_ - m_; i < (n_ - 1); ++i)
            state_[i] = state_[i + m_ - n_] ^ twiddle(state_[i], state_[i + 1]);
        state_[n_ - 1] = state_[m_ - 1] ^ twiddle(state_[n_ - 1], state_[0]);
        p_ = 0;
    }
};

}  // namespace NFcore

#endif  // NFSIM_RNG_H_
