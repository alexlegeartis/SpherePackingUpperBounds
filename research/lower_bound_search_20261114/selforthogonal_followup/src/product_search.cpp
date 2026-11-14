#include <bits/stdc++.h>
using namespace std;
using U = uint32_t;
int n, target;
vector<U> blocks;
vector<vector<int>> adj;
vector<unsigned char> bad;
vector<int> answer;
U repeat = 0;
uint64_t nodes = 0, edges = 0;
array<U, 4> triples(U a) {
  array<U, 4> r;
  int i = 0;
  for (U b = a; b;) {
    U x = b & -b;
    b -= x;
    r[i++] = a ^ x;
  }
  return r;
}
bool pairgood(U a, U b) {
  if (a & b)
    return false;
  auto x = triples(a), y = triples(b);
  for (U u : x)
    for (U v : y)
      if (bad[u | v])
        return false;
  return true;
}
bool threegood(U a, U b, U c) {
  auto x = triples(a), y = triples(b), z = triples(c);
  for (U u : x)
    for (U v : y)
      for (U w : z)
        if (bad[u | v | w])
          return false;
  return true;
}
bool finish(const vector<int> &chosen, U used) {
  if (n % 4 != 3) {
    answer = chosen;
    return true;
  }
  U r = ((1u << n) - 1) ^ used;
  assert(__builtin_popcount(r) == 3);
  for (int a : chosen)
    for (U t : triples(blocks[a]))
      if (bad[r | t])
        return false;
  for (int i = 0; i < (int)chosen.size(); ++i)
    for (int j = i + 1; j < (int)chosen.size(); ++j)
      for (U u : triples(blocks[chosen[i]]))
        for (U v : triples(blocks[chosen[j]]))
          if (bad[r | u | v])
            return false;
  answer = chosen;
  repeat = r;
  return true;
}
bool go(vector<int> &chosen, const vector<int> &cand, U used) {
  ++nodes;
  if ((int)chosen.size() == target)
    return finish(chosen, used);
  if (chosen.size() + cand.size() < target)
    return false;
  for (int j = 0; j < (int)cand.size(); ++j) {
    if (chosen.size() + cand.size() - j < target)
      break;
    int a = cand[j];
    bool okay = true;
    for (int i = 0; i < (int)chosen.size() && okay; ++i)
      for (int k = i + 1; k < (int)chosen.size(); ++k)
        if (!threegood(blocks[a], blocks[chosen[i]], blocks[chosen[k]])) {
          okay = false;
          break;
        }
    if (!okay)
      continue;
    vector<int> rest;
    set_intersection(cand.begin() + j + 1, cand.end(), adj[a].begin(),
                     adj[a].end(), back_inserter(rest));
    chosen.push_back(a);
    if (go(chosen, rest, used | blocks[a]))
      return true;
    chosen.pop_back();
  }
  return false;
}
int main(int argc, char **argv) {
  if (argc != 3)
    return 1;
  ifstream f(argv[1]);
  int dummy, a, b, c;
  f >> n >> dummy >> a >> b >> c;
  if (!f || n < 15 || n > 20)
    return 2;
  target = n / 4;
  bad.assign(1u << n, 0);
  for (int i = 0; i < a; ++i) {
    U s;
    f >> s;
    bad[s] = 1;
  }
  for (int i = 0; i < b; ++i) {
    U s;
    f >> s;
    for (U t = s; t;) {
      U x = t & -t;
      t -= x;
      bad[s ^ x] = 1;
    }
  }
  for (int i = 0; i < c; ++i) {
    U s;
    f >> s;
    bad[s] = 1;
  }
  if (!f)
    return 3;
  for (U s = 0; s < (1u << n); ++s)
    if (__builtin_popcount(s) == 4)
      blocks.push_back(s);
  adj.resize(blocks.size());
  for (int i = 0; i < (int)blocks.size(); ++i)
    for (int j = i + 1; j < (int)blocks.size(); ++j)
      if (pairgood(blocks[i], blocks[j])) {
        adj[i].push_back(j);
        ++edges;
      }
  vector<int> cand(blocks.size());
  iota(cand.begin(), cand.end(), 0);
  vector<int> chosen;
  bool sat = go(chosen, cand, 0);
  ofstream out(argv[2]);
  out << (sat ? "SAT" : "UNSAT") << " n " << n << " edges " << edges
      << " nodes " << nodes << '\n';
  if (sat) {
    for (int i : answer)
      out << blocks[i] << '\n';
    if (repeat)
      out << "REPEAT " << repeat << '\n';
  }
  return 0;
}
