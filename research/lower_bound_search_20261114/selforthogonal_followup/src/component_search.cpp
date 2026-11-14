#include <bits/stdc++.h>
using namespace std;
using U = uint32_t;
struct Component {
  array<U, 3> triple;
  U mask;
  array<U, 30> six;
};
vector<unsigned char> cost;
vector<Component> valid;
vector<U> triples;
int n;
array<U, 4> tetratriples(U b) {
  array<U, 4> r;
  U t = b;
  int i = 0;
  while (t) {
    U x = t & -t;
    t -= x;
    r[i++] = b ^ x;
  }
  return r;
}
Component component(U a, U b, U c) {
  Component x;
  x.triple = {a, b, c};
  x.mask = a | b | c;
  x.six[0] = a | b;
  x.six[1] = a | c;
  x.six[2] = b | c;
  int k = 3;
  for (U u = a; u;) {
    U i = u & -u;
    u -= i;
    for (U v = b; v;) {
      U j = v & -v;
      v -= j;
      for (U w = c; w;) {
        U l = w & -w;
        w -= l;
        x.six[k++] = x.mask ^ i ^ j ^ l;
      }
    }
  }
  return x;
}
int score(const Component &c) {
  int s = 4 * cost[c.mask];
  for (int i = 0; i < 30; ++i)
    s += (i < 3 ? 2 : 1) * cost[c.six[i]];
  return s;
}
bool combine17(const Component &c, U a, U b) {
  auto ta = tetratriples(a), tb = tetratriples(b);
  for (U x : ta)
    for (U y : tb)
      if (cost[x | y])
        return false;
  for (U x : c.triple) {
    for (U y : ta)
      if (cost[x | y])
        return false;
    for (U y : tb)
      if (cost[x | y])
        return false;
  }
  for (U x : c.six) {
    for (U y : ta)
      if (cost[x | y])
        return false;
    for (U y : tb)
      if (cost[x | y])
        return false;
  }
  for (U x : c.triple)
    for (U y : ta)
      for (U z : tb)
        if (cost[x | y | z])
          return false;
  return true;
}
bool combine18(const Component &a, const Component &b) {
  for (U x : a.triple)
    for (U y : b.triple)
      if (cost[x | y])
        return false;
  for (U x : a.triple)
    for (U y : b.six)
      if (cost[x | y])
        return false;
  for (U x : b.triple)
    for (U y : a.six)
      if (cost[x | y])
        return false;
  return true;
}
vector<vector<int>> generator(const Component &c) {
  vector<vector<int>> G(4, vector<int>(n));
  for (int j = 0; j < 3; ++j) {
    int i = 0;
    for (int k = 0; k < n; ++k)
      if (c.triple[j] >> k & 1) {
        G[j][k] = 1;
        G[3][k] = (i++ == 0 ? 1 : i == 2 ? 2 : 0);
      }
  }
  return G;
}
void addtetra(vector<vector<int>> &G, U a) {
  vector<int> x(n), y(n);
  int p = 0;
  int u[4] = {1, 1, 1, 0}, v[4] = {0, 1, 2, 1};
  for (int i = 0; i < n; ++i)
    if (a >> i & 1) {
      x[i] = u[p];
      y[i] = v[p++];
    }
  G.push_back(x);
  G.push_back(y);
}
int main(int argc, char **argv) {
  if (argc != 3 && argc != 4)
    return 1;
  bool componentonly = argc == 4 && string(argv[3]) == "component";
  ifstream f(argv[1]);
  int k, s6, s7, s9;
  f >> n >> k >> s6 >> s7 >> s9;
  if (!f || n > 18 || (componentonly ? n < 9 : n < 17))
    return 2;
  cost.assign(1u << n, 0);
  for (int i = 0; i < s6; ++i) {
    U s;
    f >> s;
    cost[s] = 10;
  }
  for (int i = 0; i < s7; ++i) {
    U s;
    f >> s;
    for (U t = s; t;) {
      U b = t & -t;
      t -= b;
      cost[s ^ b] = max<int>(1, cost[s ^ b]);
    }
  }
  for (int i = 0; i < s9; ++i) {
    U s;
    f >> s;
    cost[s] = 1;
  }
  if (!f)
    return 3;
  for (U s = 0; s < (1u << n); ++s)
    if (__builtin_popcount(s) == 3)
      triples.push_back(s);
  uint64_t tested = 0, combinations = 0;
  int least = INT_MAX;
  Component best;
  auto start = chrono::steady_clock::now();
  for (int i = 0; i < (int)triples.size(); ++i)
    for (int j = i + 1; j < (int)triples.size(); ++j) {
      U a = triples[i], b = triples[j];
      if (a & b)
        continue;
      for (int l = j + 1; l < (int)triples.size(); ++l) {
        U c = triples[l];
        if (c & (a | b))
          continue;
        auto x = component(a, b, c);
        int s = score(x);
        ++tested;
        if (s < least) {
          least = s;
          best = x;
        }
        if (!s)
          valid.push_back(x);
      }
    }
  cout << "n " << n << " components_tested " << tested
       << " least_component_score " << least << " valid_components "
       << valid.size() << " seconds "
       << chrono::duration<double>(chrono::steady_clock::now() - start).count()
       << endl;
  vector<vector<int>> answer;
  if (componentonly) {
    if (!valid.empty())
      answer = generator(valid.front());
  } else if (n == 17) {
    for (auto &c : valid) {
      U rem = ((1u << n) - 1) ^ c.mask;
      U first = rem & -rem;
      for (U a = rem; a; a = (a - 1) & rem) {
        if (__builtin_popcount(a) != 4 || !(a & first))
          continue;
        U b = rem ^ a;
        ++combinations;
        if (combine17(c, a, b)) {
          answer = generator(c);
          addtetra(answer, a);
          addtetra(answer, b);
          break;
        }
      }
      if (!answer.empty())
        break;
    }
  } else {
    unordered_map<U, vector<int>> bymask;
    for (int i = 0; i < (int)valid.size(); ++i)
      bymask[valid[i].mask].push_back(i);
    for (int i = 0; i < (int)valid.size() && answer.empty(); ++i) {
      U comp = ((1u << n) - 1) ^ valid[i].mask;
      for (int j : bymask[comp])
        if (j > i) {
          ++combinations;
          if (combine18(valid[i], valid[j])) {
            answer = generator(valid[i]);
            auto G = generator(valid[j]);
            answer.insert(answer.end(), G.begin(), G.end());
            break;
          }
        }
    }
  }
  ofstream out(argv[2]);
  out << (answer.empty() ? "UNSAT" : "SAT") << " n " << n
      << " tested_components " << tested << " minimum_component_score " << least
      << " valid_components " << valid.size() << " combinations "
      << combinations << "\n";
  if (!answer.empty()) {
    for (auto &r : answer) {
      for (int i = 0; i < n; ++i) {
        if (i)
          out << ' ';
        out << r[i];
      }
      out << '\n';
    }
  } else
    out << best.triple[0] << ' ' << best.triple[1] << ' ' << best.triple[2]
        << '\n';
  cout << (answer.empty() ? "UNSAT" : "SAT") << " combinations " << combinations
       << endl;
}
