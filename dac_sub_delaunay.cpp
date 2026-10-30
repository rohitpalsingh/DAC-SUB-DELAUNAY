// DAC-SUB-DELAUNAY: Divide-and-conquer construction of epsilon-constrained,
// sparsified beta^s-complexes (subspace Delaunay/Alpha approximations) using
// k-way Space-Partitioning trees (spherical k-means).
//
// Based on: R. P. Singh and P. A. Wilsey, "Constructing epsilon-Constrained
// Sparsified beta^s-Complexes using Space Partitioning Trees".
//
//   Algorithm 1  -> buildSPTree()      (k-way SP-tree, spherical k-means)
//   Algorithm 2  -> Mesher::rec()      (post-order incremental mesh generation)
//   Algorithm 3  -> Mesher::finish()   (sort candidates by circum-radius, validate)
//   Algorithm 4  -> Mesher::crossSimplices() (eps-neighbourhood cross-simplex enumeration)
//
// For every simplex dimension d = 1..s the program builds the set of valid
// d-simplices (d+1 vertices) of length <= eps whose beta-exclusion ball is empty:
//   beta >= 1 : ball(CC, CR)                              (beta = 1 : Alpha-like)
//   beta <  1 : ball(CC, CR * (1/beta - sqrt(1/beta^2-1))) (beta -> 0 : Vietoris-Rips)
// where CC/CR are the circum-centre/-radius of the simplex inside its own affine hull.
//
// Build:  g++ -O3 -std=c++17 -o dac_sub_delaunay dac_sub_delaunay.cpp
// Run:    ./dac_sub_delaunay --gen sphere --dim 4 --n 500 --s 3 --eps 0.5 --beta 1 --verify

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using Clock = chrono::steady_clock;
static double secs(Clock::time_point a) { return chrono::duration<double>(Clock::now() - a).count(); }

// ---------------------------------------------------------------- point set
struct Points {
    int n = 0, d = 0;
    vector<double> x;
    const double* operator[](int i) const { return &x[(size_t)i * d]; }
};
static inline double dist2(const double* a, const double* b, int d) {
    double s = 0;
    for (int i = 0; i < d; i++) { double t = a[i] - b[i]; s += t * t; }
    return s;
}

// ------------------------------------------------------- Algorithm 1: SP-tree
struct KWayTreeNode {
    vector<double> centroid;                 // average of subtree points
    double radius = 0;                       // R_encl: centroid -> farthest point
    vector<vector<double>> directionVec;     // unit direction (from centroid) of each child
    KWayTreeNode* parent = nullptr;
    int parentToChildDirection = -1;         // index of this node's direction in parent
    vector<unique_ptr<KWayTreeNode>> children;
    vector<int> data;                        // point ids in this subtree
    bool leaf() const { return children.empty(); }
};

static void sphericalKMeans(const Points& P, const vector<int>& idx, const vector<double>& ctr, int k,
                            mt19937& rng, vector<int>& assign, vector<vector<double>>& centers) {
    int m = (int)idx.size(), d = P.d;
    k = min(k, m);
    vector<double> U((size_t)m * d);
    for (int i = 0; i < m; i++) {
        const double* p = P[idx[i]];
        double nrm = 0;
        for (int j = 0; j < d; j++) { U[(size_t)i * d + j] = p[j] - ctr[j]; nrm += U[(size_t)i * d + j] * U[(size_t)i * d + j]; }
        nrm = sqrt(nrm);
        if (nrm > 1e-300) for (int j = 0; j < d; j++) U[(size_t)i * d + j] /= nrm;
    }
    auto dot = [&](int i, const vector<double>& c) { double s = 0; for (int j = 0; j < d; j++) s += U[(size_t)i * d + j] * c[j]; return s; };
    centers.assign(k, vector<double>(d, 0.0));
    int f = uniform_int_distribution<int>(0, m - 1)(rng);
    for (int j = 0; j < d; j++) centers[0][j] = U[(size_t)f * d + j];
    vector<double> D(m, 1e300);
    for (int c = 1; c < k; c++) {  // k-means++ seeding with cosine distance
        double sum = 0;
        for (int i = 0; i < m; i++) { D[i] = min(D[i], max(0.0, 1.0 - dot(i, centers[c - 1]))); sum += D[i]; }
        int pick;
        if (sum <= 1e-15) pick = uniform_int_distribution<int>(0, m - 1)(rng);
        else {
            double r = uniform_real_distribution<double>(0, sum)(rng), acc = 0;
            pick = m - 1;
            for (int i = 0; i < m; i++) { acc += D[i]; if (acc >= r) { pick = i; break; } }
        }
        for (int j = 0; j < d; j++) centers[c][j] = U[(size_t)pick * d + j];
    }
    assign.assign(m, -1);
    for (int it = 0; it < 30; it++) {
        bool changed = false;
        for (int i = 0; i < m; i++) {
            int bj = 0; double best = -1e300;
            for (int j = 0; j < k; j++) { double dp = dot(i, centers[j]); if (dp > best) { best = dp; bj = j; } }
            if (assign[i] != bj) { assign[i] = bj; changed = true; }
        }
        if (!changed) break;
        vector<vector<double>> sums(k, vector<double>(d, 0.0));
        for (int i = 0; i < m; i++) for (int j = 0; j < d; j++) sums[assign[i]][j] += U[(size_t)i * d + j];
        for (int c = 0; c < k; c++) {
            double nrm = 0; for (int j = 0; j < d; j++) nrm += sums[c][j] * sums[c][j];
            nrm = sqrt(nrm);
            if (nrm > 1e-12) for (int j = 0; j < d; j++) centers[c][j] = sums[c][j] / nrm;
        }
    }
}

static unique_ptr<KWayTreeNode> buildSPTree(const Points& P, vector<int> idx, int maxSize, int k,
                                            KWayTreeNode* parent, int dirId, mt19937& rng) {
    auto node = make_unique<KWayTreeNode>();
    node->parent = parent;
    node->parentToChildDirection = dirId;
    node->centroid.assign(P.d, 0.0);
    for (int i : idx) for (int j = 0; j < P.d; j++) node->centroid[j] += P[i][j];
    for (int j = 0; j < P.d; j++) node->centroid[j] /= (double)idx.size();
    double r2 = 0;
    for (int i : idx) r2 = max(r2, dist2(P[i], node->centroid.data(), P.d));
    node->radius = sqrt(r2);
    node->data = idx;
    if ((int)idx.size() <= maxSize) return node;  // leaf (lines 2-3)

    vector<int> assign; vector<vector<double>> centers;
    sphericalKMeans(P, idx, node->centroid, k, rng, assign, centers);
    vector<vector<int>> groups(centers.size());
    for (size_t i = 0; i < idx.size(); i++) groups[assign[i]].push_back(idx[i]);
    vector<vector<int>> parts; vector<vector<double>> dirs;
    for (size_t g = 0; g < groups.size(); g++) if (!groups[g].empty()) { parts.push_back(groups[g]); dirs.push_back(centers[g]); }
    if (parts.size() < 2) {  // degenerate (duplicates): split by halves to guarantee progress
        parts.clear(); dirs.clear();
        size_t h = idx.size() / 2;
        parts.push_back(vector<int>(idx.begin(), idx.begin() + h));
        parts.push_back(vector<int>(idx.begin() + h, idx.end()));
        dirs.assign(2, vector<double>(P.d, 0.0));
    }
    for (size_t g = 0; g < parts.size(); g++) {
        auto child = buildSPTree(P, parts[g], maxSize, k, node.get(), (int)g, rng);
        node->directionVec.push_back(dirs[g]);
        node->children.push_back(std::move(child));
    }
    return node;
}

// eps-neighbourhood search on the SP-tree; f(id) returns true to stop early.
template <class F>
static bool ballQuery(const KWayTreeNode* nd, const Points& P, const double* q, double r, F&& f) {
    double dc = sqrt(dist2(q, nd->centroid.data(), P.d));
    if (dc > (r + nd->radius) * (1 + 1e-9) + 1e-12) return false;
    if (nd->leaf()) {
        double r2 = r * r;
        for (int i : nd->data) if (dist2(q, P[i], P.d) <= r2) if (f(i)) return true;
        return false;
    }
    for (auto& c : nd->children) if (ballQuery(c.get(), P, q, r, f)) return true;
    return false;
}

// ------------------------------------------------------------ geometry / beta
static double exclusionFactor(double beta) {
    if (beta >= 1.0) return 1.0;
    if (beta <= 0.0) return 0.0;
    return 1.0 / beta - sqrt(1.0 / (beta * beta) - 1.0);
}

// Circum-centre / radius of the m points (m-1 = simplex dimension) in their affine hull.
static bool circumsphere(const Points& P, const int* v, int m, vector<double>& cc, double& cr) {
    int d = P.d, q = m - 1;
    const double* p0 = P[v[0]];
    vector<vector<double>> a(q, vector<double>(d));
    for (int i = 0; i < q; i++) for (int j = 0; j < d; j++) a[i][j] = P[v[i + 1]][j] - p0[j];
    vector<vector<double>> G(q, vector<double>(q + 1));
    double scale = 0;
    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) { double s = 0; for (int t = 0; t < d; t++) s += a[i][t] * a[j][t]; G[i][j] = s; }
        G[i][q] = 0.5 * G[i][i];
        scale = max(scale, G[i][i]);
    }
    if (scale <= 0) return false;
    for (int c = 0; c < q; c++) {  // Gaussian elimination, partial pivoting
        int piv = c;
        for (int r = c + 1; r < q; r++) if (fabs(G[r][c]) > fabs(G[piv][c])) piv = r;
        if (fabs(G[piv][c]) < 1e-12 * scale) return false;  // degenerate simplex
        swap(G[c], G[piv]);
        for (int r = c + 1; r < q; r++) {
            double f = G[r][c] / G[c][c];
            for (int t = c; t <= q; t++) G[r][t] -= f * G[c][t];
        }
    }
    vector<double> coef(q);
    for (int r = q - 1; r >= 0; r--) {
        double s = G[r][q];
        for (int t = r + 1; t < q; t++) s -= G[r][t] * coef[t];
        coef[r] = s / G[r][r];
    }
    cc.assign(d, 0.0);
    double r2 = 0;
    for (int j = 0; j < d; j++) {
        double off = 0;
        for (int i = 0; i < q; i++) off += coef[i] * a[i][j];
        cc[j] = p0[j] + off;
        r2 += off * off;
    }
    cr = sqrt(r2);
    return true;
}

// --------------------------------------------------- clique enumeration helper
static void enumCliques(const vector<vector<int>>& adj, int m, const function<void(const vector<int>&)>& cb) {
    vector<int> cur;
    function<void(const vector<int>&)> rec = [&](const vector<int>& cand) {
        for (size_t a = 0; a < cand.size(); a++) {
            int c = cand[a];
            if ((int)cur.size() + 1 == m) { cur.push_back(c); cb(cur); cur.pop_back(); continue; }
            vector<int> nc;
            set_intersection(cand.begin() + a + 1, cand.end(), adj[c].begin(), adj[c].end(), back_inserter(nc));
            if (cur.size() + 1 + nc.size() < (size_t)m) continue;
            cur.push_back(c); rec(nc); cur.pop_back();
        }
    };
    for (size_t t = 0; t < adj.size(); t++) { cur.assign(1, (int)t); rec(adj[t]); }
}

struct Cand { vector<int> v; double cr; vector<double> cc; };
struct Simplex { vector<int> v; double w; };

// ------------------------------------------------ Algorithms 2, 3, 4: mesher
struct Mesher {
    const Points& P; const KWayTreeNode* root; double eps, beta; int d;  // d = simplex dimension
    size_t candidates = 0;
    vector<int> pos;
    Mesher(const Points& P_, const KWayTreeNode* r, double e, double b, int dim)
        : P(P_), root(r), eps(e), beta(b), d(dim), pos(P_.n, -1) {}

    bool valid(const Cand& c) const {  // beta-criterion via SP-tree range query
        double f = exclusionFactor(beta);
        if (f == 0.0) return true;
        double r = c.cr * f * (1 - 1e-9);
        return !ballQuery(root, P, c.cc.data(), r, [](int) { return true; });
    }

    void addCand(const vector<int>& verts, vector<Cand>& out) {
        Cand c; c.v = verts;
        if (!circumsphere(P, verts.data(), (int)verts.size(), c.cc, c.cr)) return;
        out.push_back(std::move(c));
    }

    // Algorithm 3 (lines 4-8): sort by circum-radius, keep validated simplices.
    void finish(vector<Cand>& cands, vector<Simplex>& out) {
        candidates += cands.size();
        sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.cr < b.cr; });
        for (auto& c : cands) if (valid(c)) out.push_back({c.v, c.cr});
    }

    // LOCALMESH: brute-force enumeration inside a leaf partition.
    void localMesh(const KWayTreeNode* nd, vector<Simplex>& out) {
        vector<int> pool = nd->data;
        sort(pool.begin(), pool.end());
        int L = (int)pool.size();
        if (L < d + 1) return;
        vector<vector<int>> adj(L);
        double e2 = eps * eps;
        for (int a = 0; a < L; a++) for (int b = a + 1; b < L; b++)
            if (dist2(P[pool[a]], P[pool[b]], P.d) <= e2) adj[a].push_back(b);
        vector<Cand> cands;
        enumCliques(adj, d + 1, [&](const vector<int>& loc) {
            vector<int> v(loc.size());
            for (size_t i = 0; i < loc.size(); i++) v[i] = pool[loc[i]];
            addCand(v, cands);
        });
        finish(cands, out);
    }

    // Algorithm 4: cross-partition simplices at an internal node, restricted to the
    // eps-neighbourhood: only points having an eps-neighbour in a different child
    // ("exposed" points) can be vertices of a cross simplex.
    void crossSimplices(const KWayTreeNode* nd, vector<Cand>& cands) {
        int nc = (int)nd->children.size();
        vector<int> lab_of; // label by child, stored through a map from point id
        vector<int> E;
        for (int ci = 0; ci < nc; ci++)
            for (int p : nd->children[ci]->data) {
                bool exposed = false;
                for (int cj = 0; cj < nc && !exposed; cj++) {
                    if (cj == ci) continue;
                    exposed = ballQuery(nd->children[cj].get(), P, P[p], eps, [](int) { return true; });
                }
                if (exposed) E.push_back(p);
            }
        if ((int)E.size() < d + 1) return;
        sort(E.begin(), E.end());
        for (size_t t = 0; t < E.size(); t++) pos[E[t]] = (int)t;
        vector<int> lab(P.n, -1);
        for (int ci = 0; ci < nc; ci++) for (int p : nd->children[ci]->data) lab[p] = ci;
        vector<vector<int>> adj(E.size());
        for (size_t t = 0; t < E.size(); t++) {
            for (int cj = 0; cj < nc; cj++)
                ballQuery(nd->children[cj].get(), P, P[E[t]], eps, [&](int j) {
                    int pj = pos[j];
                    if (pj > (int)t) adj[t].push_back(pj);
                    return false;
                });
            sort(adj[t].begin(), adj[t].end());
        }
        enumCliques(adj, d + 1, [&](const vector<int>& loc) {
            int l0 = lab[E[loc[0]]]; bool cross = false;
            for (int x : loc) if (lab[E[x]] != l0) { cross = true; break; }
            if (!cross) return;
            vector<int> v(loc.size());
            for (size_t i = 0; i < loc.size(); i++) v[i] = E[loc[i]];
            addCand(v, cands);
        });
        for (int e : E) pos[e] = -1;
    }

    // Algorithm 2: post-order traversal; local meshes at leaves, merge on the way up.
    void rec(const KWayTreeNode* nd, vector<Simplex>& out) {
        if (nd->leaf()) { localMesh(nd, out); return; }
        for (auto& c : nd->children) rec(c.get(), out);
        vector<Cand> cands;
        crossSimplices(nd, cands);
        finish(cands, out);
    }
};

// -------------------------------------------------------- brute-force baseline
static set<vector<int>> bruteForce(const Points& P, double eps, double beta, int d) {
    int n = P.n; double e2 = eps * eps, f = exclusionFactor(beta);
    vector<vector<int>> adj(n);
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (dist2(P[a], P[b], P.d) <= e2) adj[a].push_back(b);
    set<vector<int>> res;
    enumCliques(adj, d + 1, [&](const vector<int>& v) {
        vector<double> cc; double cr;
        if (!circumsphere(P, v.data(), d + 1, cc, cr)) return;
        bool ok = true;
        if (f != 0.0) {
            double r = cr * f * (1 - 1e-9), r2 = r * r;
            for (int i = 0; i < n && ok; i++) if (dist2(cc.data(), P[i], P.d) <= r2) ok = false;
        }
        if (ok) res.insert(v);
    });
    return res;
}

// ------------------------------------------------------------ data generators
static Points genSphere(int n, int D, double noise, mt19937& rng) {  // S^{D-1} in R^D
    Points P; P.n = n; P.d = D; P.x.resize((size_t)n * D);
    normal_distribution<double> g(0, 1);
    for (int i = 0; i < n; i++) {
        double nr = 0; vector<double> v(D);
        for (int j = 0; j < D; j++) { v[j] = g(rng); nr += v[j] * v[j]; }
        nr = sqrt(nr);
        for (int j = 0; j < D; j++) P.x[(size_t)i * D + j] = v[j] / nr + noise * g(rng);
    }
    return P;
}
static Points genTorus(int n, double noise, mt19937& rng) {
    Points P; P.n = n; P.d = 3; P.x.resize((size_t)n * 3);
    uniform_real_distribution<double> U(0, 2 * M_PI); normal_distribution<double> g(0, 1);
    for (int i = 0; i < n; i++) {
        double u = U(rng), v = U(rng), R = 2, r = 1;
        P.x[3 * i] = (R + r * cos(v)) * cos(u) + noise * g(rng);
        P.x[3 * i + 1] = (R + r * cos(v)) * sin(u) + noise * g(rng);
        P.x[3 * i + 2] = r * sin(v) + noise * g(rng);
    }
    return P;
}
static Points genRandom(int n, int D, mt19937& rng) {
    Points P; P.n = n; P.d = D; P.x.resize((size_t)n * D);
    uniform_real_distribution<double> U(0, 1);
    for (auto& v : P.x) v = U(rng);
    return P;
}
static bool loadCSV(const string& fn, Points& P) {
    ifstream in(fn);
    if (!in) return false;
    string line; vector<double> all; int d = -1;
    while (getline(in, line)) {
        for (char& c : line) if (c == ',' || c == ';' || c == '\t') c = ' ';
        istringstream ss(line); vector<double> row; string tok; bool bad = false;
        while (ss >> tok) { char* e; double v = strtod(tok.c_str(), &e); if (*e) { bad = true; break; } row.push_back(v); }
        if (bad || row.empty()) continue;  // header / blank / comment
        if (d < 0) d = (int)row.size();
        if ((int)row.size() != d) continue;
        all.insert(all.end(), row.begin(), row.end());
    }
    if (d <= 0) return false;
    P.d = d; P.n = (int)(all.size() / d); P.x = std::move(all);
    return P.n > 0;
}

// ------------------------------------------------------------------------ main
static void usage() {
    puts("Usage: dac_sub_delaunay [options]\n"
         "  --input FILE      CSV/whitespace point cloud (rows = points)\n"
         "  --gen TYPE        sphere | torus | random   (default sphere)\n"
         "  --n N --dim D     synthetic size / ambient dimension (sphere: S^{D-1} in R^D)\n"
         "  --noise X         Gaussian noise for synthetic data (default 0.02)\n"
         "  --s S             max simplex dimension, builds dims 1..S (default 2)\n"
         "  --eps E           epsilon_max connectivity bound (default 0.6)\n"
         "  --beta B          sparsification factor, 0 (VR) .. 1 (Alpha-like) (default 1)\n"
         "  --k K --leaf M    SP-tree fan-out and leaf size (default 4, 32)\n"
         "  --seed S          RNG seed\n"
         "  --verify          compare against brute-force construction\n"
         "  --out FILE        write simplices as: dim weight v0 v1 ...");
}

int main(int argc, char** argv) {
    map<string, string> a; set<string> flags;
    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s == "--verify" || s == "--help") flags.insert(s);
        else if (s.rfind("--", 0) == 0 && i + 1 < argc) a[s] = argv[++i];
        else { usage(); return 1; }
    }
    if (flags.count("--help")) { usage(); return 0; }
    auto geti = [&](const string& k, int def) { return a.count(k) ? atoi(a[k].c_str()) : def; };
    auto getd = [&](const string& k, double def) { return a.count(k) ? atof(a[k].c_str()) : def; };
    mt19937 rng(geti("--seed", 1));
    int n = geti("--n", 300), D = geti("--dim", 3), s = geti("--s", 2), k = max(2, geti("--k", 4)), leaf = max(2, geti("--leaf", 32));
    double eps = getd("--eps", 0.6), beta = getd("--beta", 1.0), noise = getd("--noise", 0.02);

    Points P;
    if (a.count("--input")) {
        if (!loadCSV(a["--input"], P)) { fprintf(stderr, "cannot read %s\n", a["--input"].c_str()); return 1; }
    } else {
        string g = a.count("--gen") ? a["--gen"] : "sphere";
        if (g == "sphere") P = genSphere(n, D, noise, rng);
        else if (g == "torus") P = genTorus(n, noise, rng);
        else if (g == "random") P = genRandom(n, D, rng);
        else { usage(); return 1; }
    }
    if (s < 1) s = 1;
    printf("points n=%d ambient dim=%d | s=%d eps=%.4g beta=%.4g k=%d leaf=%d\n", P.n, P.d, s, eps, beta, k, leaf);

    auto t0 = Clock::now();
    vector<int> all(P.n); iota(all.begin(), all.end(), 0);
    auto root = buildSPTree(P, all, leaf, k, nullptr, -1, rng);
    printf("SP-tree built in %.3f s\n", secs(t0));

    ofstream out; if (a.count("--out")) out.open(a["--out"]);
    vector<Simplex> edges;
    bool allOK = true;
    for (int d = 1; d <= s; d++) {
        auto t1 = Clock::now();
        Mesher M(P, root.get(), eps, beta, d);
        vector<Simplex> mesh;
        M.rec(root.get(), mesh);
        double tm = secs(t1);
        printf("dim %d: %zu valid simplices (%zu candidates)  [%.3f s]\n", d, mesh.size(), M.candidates, tm);
        if (d == 1) edges = mesh;
        if (out.is_open()) for (auto& sx : mesh) { out << d << ' ' << sx.w; for (int v : sx.v) out << ' ' << v; out << '\n'; }
        if (flags.count("--verify")) {
            auto t2 = Clock::now();
            auto ref = bruteForce(P, eps, beta, d);
            set<vector<int>> got; for (auto& sx : mesh) got.insert(sx.v);
            bool ok = (got == ref);
            allOK &= ok;
            printf("   verify vs brute force: %s (brute force %zu simplices, %.3f s)\n", ok ? "PASS" : "FAIL", ref.size(), secs(t2));
        }
    }

    // H0 summary from the beta^1 mesh (edge length = 2 * circum-radius), union-find.
    vector<int> uf(P.n); iota(uf.begin(), uf.end(), 0);
    function<int(int)> find = [&](int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
    vector<pair<double, pair<int, int>>> E;
    for (auto& e : edges) E.push_back({2 * e.w, {e.v[0], e.v[1]}});
    sort(E.begin(), E.end());
    vector<double> deaths; int comps = P.n;
    for (auto& e : E) { int x = find(e.second.first), y = find(e.second.second); if (x != y) { uf[x] = y; comps--; deaths.push_back(e.first); } }
    sort(deaths.rbegin(), deaths.rend());
    printf("H0: %d connected component(s) at eps=%.4g; largest finite death values:", comps, eps);
    for (size_t i = 0; i < min<size_t>(5, deaths.size()); i++) printf(" %.4f", deaths[i]);
    printf("\ntotal time %.3f s\n", secs(t0));
    if (flags.count("--verify")) { printf("VERIFY: %s\n", allOK ? "ALL PASS" : "FAILED"); return allOK ? 0 : 2; }
    return 0;
}
