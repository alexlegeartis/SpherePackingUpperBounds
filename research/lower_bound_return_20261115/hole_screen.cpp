// Discovery only. Nonempty results carry exactly checked upper witnesses.
// A capped or empty floating traversal does not certify a hole.
#include <bits/stdc++.h>
using namespace std;
constexpr int N = 64;
long long Q[N][N], z[N], c[N], witness[N], exact_value;
double mu[N][N], d[N], zd[N], radius2, pruning_slack = -1;
uint64_t nodes, cap;
int residue_modulus;
bool found, capped;

bool walk(int i, double partial) {
    if (++nodes > cap) { capped = true; return true; }
    if (i < 0) {
        __int128 value = 0;
        for (int a = 0; a < N; ++a)
            for (int b = 0; b < N; ++b)
                value += (__int128)z[a] * Q[a][b] * z[b];
        if (value <= (long long)radius2) {
            exact_value = (long long)value;
            copy(z, z + N, witness);
            found = true;
            return true;
        }
        return false;
    }
    double center = 0;
    for (int j = i + 1; j < N; ++j) center -= mu[i][j] * zd[j];
    double fraction = 1.;
    if (pruning_slack >= 0) {
        double j = N-i;
        fraction = min(1., j/N + pruning_slack * sqrt(2*j*(N-j)/(N*N*(N+2.))));
    }
    double remaining = radius2 * fraction + 1e-6 - partial;
    if (remaining < 0) return false;
    long long base = c[i] + residue_modulus * llround((center - c[i]) / residue_modulus);
    double delta = base - center;
    if (d[i] * delta * delta > remaining) return false;
    z[i] = base; zd[i] = (double)base;
    if (walk(i - 1, partial + d[i] * delta * delta)) return true;
    int sign = delta >= 0 ? -1 : 1;
    for (long long step = 1;; ++step) {
        bool any = false;
        for (int side : {1, -1}) {
            long long x = base + side * sign * residue_modulus * step;
            double v = x - center;
            if (d[i] * v * v <= remaining) {
                any = true; z[i] = x; zd[i] = (double)x;
                if (walk(i - 1, partial + d[i] * v * v)) return true;
            }
        }
        if (!any) break;
    }
    return false;
}

int main(int argc, char **argv) {
    if (argc != 7 && argc != 8) return 2;
    ifstream gram(argv[1]), input(argv[2]);
    ofstream output(argv[3]);
    cap = stoull(argv[4]); residue_modulus = stoi(argv[5]); radius2 = stoi(argv[6]);
    if (argc == 8) pruning_slack = stod(argv[7]);
    if (residue_modulus < 2 || residue_modulus > 3 || radius2 <= 0) return 2;
    for (auto &row : Q) for (auto &v : row) gram >> v;
    if (!gram) return 3;
    long double L[N][N] = {}, D[N] = {};
    for (int i = 0; i < N; ++i) {
        L[i][i] = 1;
        for (int j = 0; j < i; ++j) {
            long double v = Q[i][j];
            for (int k = 0; k < j; ++k) v -= L[i][k] * D[k] * L[j][k];
            L[i][j] = v / D[j];
        }
        D[i] = Q[i][i];
        for (int j = 0; j < i; ++j) D[i] -= L[i][j] * L[i][j] * D[j];
        if (D[i] <= 0) return 4;
        d[i] = (double)D[i];
        for (int j = 0; j < i; ++j) mu[j][i] = (double)L[i][j];
    }
    int id, processed = 0, rejected = 0;
    uint64_t total_nodes = 0;
    auto start = chrono::steady_clock::now();
    while (input >> id) {
        for (auto &v : c) input >> v;
        if (!input) return 5;
        nodes = 0; found = capped = false; exact_value = -1;
        walk(N - 1, 0);
        total_nodes += nodes; processed++; rejected += found;
        output << "{\"id\":" << id << ",\"status\":\""
               << (found ? "exact_reject" : capped ? "node_cap" : pruning_slack >= 0 ? "empty_pruned_candidate" : "empty_float_candidate")
               << "\",\"nodes\":" << nodes << ",\"norm\":" << exact_value << ",\"witness\":[";
        if (found) for (int i = 0; i < N; ++i) output << (i ? "," : "") << witness[i];
        output << "]}\n"; output.flush();
        if (processed % 100 == 0)
            cout << "processed " << processed << " rejected " << rejected
                 << " nodes " << total_nodes << " seconds "
                 << chrono::duration<double>(chrono::steady_clock::now() - start).count() << endl;
    }
    cout << "done " << processed << " rejected " << rejected << " nodes " << total_nodes << endl;
}
