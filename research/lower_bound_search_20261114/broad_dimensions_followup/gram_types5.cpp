// Exact signed-permutation class generation for 5x5 min-root Gram candidates.
#include <bits/stdc++.h>
using namespace std;
long long determinant(vector<vector<long long>> A) {
  int n = A.size();
  long long previous = 1;
  for (int k = 0; k < n - 1; k++) {
    if (A[k][k] == 0)
      return 0;
    long long p = A[k][k];
    for (int i = k + 1; i < n; i++)
      for (int j = k + 1; j < n; j++)
        A[i][j] = (A[i][j] * p - A[i][k] * A[k][j]) / previous;
    previous = p;
  }
  return A[n - 1][n - 1];
}
int main(int argc, char **argv) {
  if (argc != 2 && argc != 3)
    return 1;
  const int n = 5, mu = argc == 3 ? atoi(argv[2]) : 4, base = mu + 1,
            total = (int)llround(pow(base, 10));
  if (mu != 4 && mu != 6)
    return 2;
  vector<uint8_t> seen(total);
  vector<pair<int, int>> edges;
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      edges.push_back({i, j});
  array<int, n> perm = {0, 1, 2, 3, 4};
  vector<array<int, n>> perms;
  do {
    perms.push_back(perm);
  } while (next_permutation(perm.begin(), perm.end()));
  ofstream out(argv[1]);
  int classes = 0, kept = 0;
  for (int code = 0; code < total; code++) {
    if (seen[code])
      continue;
    classes++;
    int q = code;
    int K[n][n]{};
    for (int i = 0; i < n; i++)
      K[i][i] = mu;
    for (auto [i, j] : edges) {
      K[i][j] = K[j][i] = q % base - mu / 2;
      q /= base;
    }
    for (auto &p : perms)
      for (int mask = 0; mask < 16; mask++) {
        int s[n] = {1, 1, 1, 1, 1};
        for (int i = 1; i < n; i++)
          if (mask >> (i - 1) & 1)
            s[i] = -1;
        int key = 0, pow = 1;
        for (auto [i, j] : edges) {
          key += (K[p[i]][p[j]] * s[i] * s[j] + mu / 2) * pow;
          pow *= base;
        }
        seen[key] = 1;
      }
    bool okay = true;
    for (int encoded = 1; encoded < 243 && okay; encoded++) {
      int t = encoded;
      int x[n];
      for (int i = 0; i < n; i++) {
        x[i] = t % 3 - 1;
        t /= 3;
      }
      int norm = 0;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          norm += x[i] * K[i][j] * x[j];
      bool nonzero = false;
      for (auto a : x)
        nonzero |= a != 0;
      if (nonzero && norm < mu)
        okay = false;
    }
    // encoded zero is (-1,-1,-1,-1,-1), and must also be checked.
    int allnorm = 0;
    for (auto &r : K)
      for (int a : r)
        allnorm += a;
    if (allnorm < mu)
      okay = false;
    for (int d = 1; d <= n && okay; d++) {
      vector<vector<long long>> A(d, vector<long long>(d));
      for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
          A[i][j] = K[i][j];
      if (determinant(A) <= 0)
        okay = false;
    }
    if (okay) {
      kept++;
      for (auto &r : K)
        for (int a : r)
          out << a << ' ';
      out << '\n';
    }
  }
  cout << "classes " << classes << " relaxed_types " << kept << endl;
}
