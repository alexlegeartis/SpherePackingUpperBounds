#include <bits/stdc++.h>
using namespace std;
using U = uint32_t;
struct Word {
  U mask;
  int mult;
};
int main(int argc, char **argv) {
  if (argc != 6)
    return 1;
  ifstream in(argv[1]), catalog(argv[2]);
  int n, k, s6, s7, s9;
  in >> n >> k >> s6 >> s7 >> s9;
  vector<unsigned char> cost(1u << n, 0);
  for (int i = 0; i < s6; ++i) {
    U s;
    in >> s;
    cost[s] = 10;
  }
  for (int i = 0; i < s7; ++i) {
    U s;
    in >> s;
    for (U t = s; t;) {
      U b = t & -t;
      t -= b;
      cost[s ^ b] = max<int>(cost[s ^ b], 1);
    }
  }
  for (int i = 0; i < s9; ++i) {
    U s;
    in >> s;
    cost[s] = 1;
  }
  int count, nc;
  catalog >> count >> nc;
  if (nc != n || !in)
    return 2;
  int steps = atoi(argv[4]);
  mt19937_64 rng(strtoull(argv[5], 0, 10));
  string prefix = argv[3];
  ofstream summary(prefix + ".summary");
  int global = INT_MAX;
  auto start = chrono::steady_clock::now();
  for (int ci = 0; ci < count; ++ci) {
    vector<vector<int>> G(k, vector<int>(n));
    for (auto &r : G)
      for (int &v : r)
        catalog >> v;
    if (!catalog)
      return 3;
    map<U, int> counts;
    int N = 1;
    for (int i = 0; i < k; ++i)
      N *= 3;
    for (int a = 1; a < N; ++a) {
      int q = a, highest = 0;
      vector<int> coeff(k);
      for (int i = 0; i < k; ++i) {
        coeff[i] = q % 3;
        if (coeff[i])
          highest = coeff[i];
        q /= 3;
      }
      if (highest != 1)
        continue;
      U mask = 0;
      for (int j = 0; j < n; ++j) {
        int v = 0;
        for (int i = 0; i < k; ++i)
          v += coeff[i] * G[i][j];
        if (v % 3)
          mask |= 1u << j;
      }
      int w = __builtin_popcount(mask);
      if (w == 6 || w == 9)
        ++counts[mask];
    }
    vector<Word> source;
    for (auto [m, c] : counts)
      source.push_back({m, c});
    vector<int> perm(n), inverse(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng);
    for (int i = 0; i < n; ++i)
      inverse[perm[i]] = i;
    vector<Word> words;
    int energy = 0;
    for (auto w : source) {
      U mask = 0;
      for (int j = 0; j < n; ++j)
        if (w.mask >> j & 1)
          mask |= 1u << perm[j];
      words.push_back({mask, w.mult});
      energy += w.mult * cost[mask];
    }
    int best = energy;
    vector<int> bestinverse = inverse;
    int beststep = 0;
    int tabu[20][20]{};
    for (int step = 1; step <= steps && best; ++step) {
      int delta[20][20]{};
      for (auto w : words) {
        int old = cost[w.mask];
        U zeros = ((1u << n) - 1) ^ w.mask;
        for (U a = w.mask; a;) {
          U i = a & -a;
          a -= i;
          int x = __builtin_ctz(i);
          for (U b = zeros; b;) {
            U j = b & -b;
            b -= j;
            int y = __builtin_ctz(j);
            delta[min(x, y)][max(x, y)] +=
                w.mult * ((int)cost[w.mask ^ i ^ j] - old);
          }
        }
      }
      int dmin = INT_MAX, bi = -1, bj = -1, ties = 0;
      for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
          int d = delta[i][j];
          if (tabu[i][j] > step && energy + d >= best)
            continue;
          if (d < dmin) {
            dmin = d;
            bi = i;
            bj = j;
            ties = 1;
          } else if (d == dmin && rng() % ++ties == 0) {
            bi = i;
            bj = j;
          }
        }
      if (bi < 0)
        break;
      U flip = (1u << bi) | (1u << bj);
      for (auto &w : words)
        if (((w.mask >> bi) ^ (w.mask >> bj)) & 1)
          w.mask ^= flip;
      energy += dmin;
      swap(inverse[bi], inverse[bj]);
      tabu[bi][bj] = step + 5 + rng() % 15;
      if (energy < best) {
        best = energy;
        bestinverse = inverse;
        beststep = step;
      }
    }
    string path = prefix + "_c" + to_string(ci) + ".txt";
    ofstream out(path);
    out << n << ' ' << k << ' ' << best << ' ' << beststep << '\n';
    for (auto &r : G) {
      for (int j = 0; j < n; ++j) {
        if (j)
          out << ' ';
        out << r[bestinverse[j]];
      }
      out << '\n';
    }
    out.close();
    if (best < global) {
      global = best;
      ifstream src(path);
      ofstream dst(prefix + ".best.txt");
      dst << src.rdbuf();
    }
    double seconds =
        chrono::duration<double>(chrono::steady_clock::now() - start).count();
    summary << ci << ' ' << best << ' ' << source.size() << ' ' << seconds
            << '\n';
    summary.flush();
    cout << "class " << ci << " best " << best << " global " << global
         << " words " << source.size() << " seconds " << seconds << endl;
    if (!best)
      break;
  }
}
