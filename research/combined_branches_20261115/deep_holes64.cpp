// Discovery only: Ne64 Gram is unscaled; capped distances are upper bounds.
#include <bits/stdc++.h>
using namespace std;
constexpr int N = 64;
long long Q[N][N], z[N], shift[N], besty[N];
double R[N][N], best;
uint64_t nodes, cap;
bool rejected;
int DEN = 3;
long long CUT = 1120;
long long exactnorm(const long long *y) {
  long long s = 0;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      s += y[i] * Q[i][j] * y[j];
  return s;
}
void visit(int k, double partial) {
  if (rejected || ++nodes > cap)
    return;
  if (k < 0) {
    long long y[N];
    for (int i = 0; i < N; i++)
      y[i] = DEN * z[i] + shift[i];
    auto value = exactnorm(y);
    if (value < best * DEN * DEN - 1e-5) {
      best = value / (double)(DEN * DEN);
      copy(y, y + N, besty);
      if (value <= CUT)
        rejected = true;
    }
    return;
  }
  double center = -shift[k] / (double)DEN;
  for (int j = k + 1; j < N; j++)
    center -= R[k][j] * (z[j] + shift[j] / (double)DEN) / R[k][k];
  double rad = sqrt(max(0., best - partial)) / R[k][k];
  long long lo = ceil(center - rad - 1e-8), hi = floor(center + rad + 1e-8),
            near = llround(center);
  for (long long step = 0; step <= max(near - lo, hi - near); step++)
    for (int side = 0; side < (step ? 2 : 1); side++) {
      long long x = near + (side ? -step : step);
      if (x < lo || x > hi)
        continue;
      z[k] = x;
      double u = R[k][k] * (x - center);
      if (partial + u * u < best + 1e-7)
        visit(k - 1, partial + u * u);
      if (rejected || nodes > cap)
        return;
    }
}
int main(int argc, char **argv) {
  if (argc != 7 && argc != 8)
    return 2;
  if (getenv("DEEP_DENOMINATOR"))
    DEN = atoi(getenv("DEEP_DENOMINATOR"));
  if (getenv("DEEP_CUT"))
    CUT = atoll(getenv("DEEP_CUT"));
  if (DEN < 2 || DEN > 20 || CUT <= 0)
    return 5;
  ifstream qf(argv[1]), rf(argv[2]);
  for (auto &row : Q)
    for (auto &x : row)
      qf >> x;
  for (auto &row : R)
    for (auto &x : row)
      rf >> x;
  if (!qf || !rf)
    return 3;
  uint64_t trials = strtoull(argv[4], 0, 10);
  cap = strtoull(argv[5], 0, 10);
  mt19937_64 rng(strtoull(argv[6], 0, 10));
  ofstream out(argv[3]), witness(string(argv[3]) + ".witnesses");
  ifstream supplied;
  if (argc == 8)
    supplied.open(argv[7]);
  uint64_t totalnodes = 0, eligible = 0, processed = 0;
  long long highest = 0;
  auto t = chrono::steady_clock::now();
  for (uint64_t trial = 0; trial < trials; trial++) {
    uint64_t id = trial;
    if (argc == 8) {
      long long previous, nodesold, completeold;
      if (!(supplied >> id >> previous >> nodesold >> completeold))
        break;
      for (auto &x : shift)
        supplied >> x;
      if (!supplied)
        return 4;
    } else
      for (auto &x : shift)
        x = (int)(rng() % DEN) - DEN / 2;
    processed++;
    for (int k = N - 1; k >= 0; k--) {
      double center = -shift[k] / (double)DEN;
      for (int j = k + 1; j < N; j++)
        center -= R[k][j] * (z[j] + shift[j] / (double)DEN) / R[k][k];
      z[k] = llround(center);
    }
    for (int i = 0; i < N; i++)
      besty[i] = DEN * z[i] + shift[i];
    long long value = exactnorm(besty);
    best = value / (double)(DEN * DEN);
    rejected = value <= CUT;
    nodes = 0;
    if (!rejected) {
      visit(N - 1, 0);
      value = exactnorm(besty);
    }
    totalnodes += nodes;
    if (value > CUT) {
      eligible++;
      out << id << ' ' << value << ' ' << nodes << ' ' << (nodes <= cap) << ' ';
      for (auto x : shift)
        out << x << ' ';
      out << '\n';
      out.flush();
    }
    if (argc == 8) {
      witness << id << ' ' << value << ' ' << nodes << ' ' << (nodes <= cap)
              << ' ';
      for (auto x : besty)
        witness << x << ' ';
      witness << '\n';
    }
    if (value > highest) {
      highest = value;
      cout << "best_upper_bound " << value << " physical "
           << value / (1. * DEN * DEN) << " trial " << trial << " nodes "
           << nodes << " uncapped " << (nodes <= cap) << endl;
    }
    if (trial % 10000 == 9999)
      cout << "progress " << trial + 1 << " candidates " << eligible
           << " nodes " << totalnodes << " seconds "
           << chrono::duration<double>(chrono::steady_clock::now() - t).count()
           << endl;
  }
  cout << "done trials " << processed << " denominator " << DEN << " cut "
       << CUT << " candidates " << eligible << " nodes " << totalnodes
       << " seconds "
       << chrono::duration<double>(chrono::steady_clock::now() - t).count()
       << endl;
}
