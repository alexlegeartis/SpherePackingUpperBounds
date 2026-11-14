// Floating discovery oracle only. A node cap never certifies a distance.
#include <algorithm>
#include <cmath>
#include <cstdint>
static int n;
static double Q[64][64], R[64][64], s[64], best;
static long long z[64], bestz[64];
static uint64_t nodes, cap;
static double norm() {
  double v = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      v += (s[i] - z[i]) * Q[i][j] * (s[j] - z[j]);
  return v;
}
static void go(int k, double partial) {
  if (++nodes > cap)
    return;
  if (k < 0) {
    double value = norm();
    if (value < best - 1e-8) {
      best = value;
      std::copy(z, z + n, bestz);
    }
    return;
  }
  double center = s[k];
  for (int j = k + 1; j < n; j++)
    center -= R[k][j] * (z[j] - s[j]) / R[k][k];
  double rad = sqrt(std::max(0., best - partial)) / R[k][k];
  long long lo = ceil(center - rad - 1e-8), hi = floor(center + rad + 1e-8),
            near = llround(center);
  for (long long step = 0; step <= std::max(near - lo, hi - near); step++)
    for (int side = 0; side < (step ? 2 : 1); side++) {
      long long x = near + (side ? -step : step);
      if (x < lo || x > hi)
        continue;
      z[k] = x;
      double t = R[k][k] * (x - center);
      if (partial + t * t < best + 1e-7)
        go(k - 1, partial + t * t);
      if (nodes > cap)
        return;
    }
}
extern "C" void init_cvp(int dim, const double *q, const double *r) {
  n = dim;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      Q[i][j] = q[i * n + j];
      R[i][j] = r[i * n + j];
    }
}
extern "C" uint64_t nearest_cvp(const double *shift, uint64_t limit,
                                double *distance, int64_t *out) {
  std::copy(shift, shift + n, s);
  cap = limit;
  nodes = 0;
  for (int k = n - 1; k >= 0; k--) {
    double center = s[k];
    for (int j = k + 1; j < n; j++)
      center -= R[k][j] * (z[j] - s[j]) / R[k][k];
    z[k] = llround(center);
  }
  best = norm();
  std::copy(z, z + n, bestz);
  go(n - 1, 0);
  *distance = best;
  for (int i = 0; i < n; i++)
    out[i] = bestz[i];
  return nodes;
}
