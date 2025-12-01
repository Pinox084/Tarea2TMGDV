#include <bits/stdc++.h>
using namespace std;

/*
 * Implementación del sketch de cuantiles MRL
 * para la Tarea 2.
 *
 * Uso:
 * ./mrl <n> <epsilon> <archivo>
 *
 * Luego de construir el sketch, entra en modo interactivo:
 * rank x
 * select r
 * quantile phi
 * exit
 */

class MRLSketch {
public:
    using Value = long long;

    MRLSketch(long long n, double epsilon)
        : n_(n), eps_(epsilon)
    {
        if (n_ <= 0) {
            throw invalid_argument("n debe ser positivo");
        }
        if (eps_ <= 0.0 || eps_ >= 1.0) {
            throw invalid_argument("epsilon debe ser en (0,1)");
        }

        // k = eps^{-1} * ceil(log2(eps*n)) + 1 (logs base 2)
        double t = eps_ * static_cast<double>(n_);
        if (t <= 0.0) {
            throw runtime_error("epsilon * n <= 0, parametros invalidos");
        }

        int logTerm = static_cast<int>(ceil(log2(t)));
        if (logTerm < 1) logTerm = 1;

        double k_real = (1.0 / eps_) * static_cast<double>(logTerm) + 1.0;
        k_ = static_cast<int>(ceil(k_real));

        // k par
        if (k_ % 2 != 0) k_ += 1;
        if (k_ <= 0) k_ = 2;
        if (k_ > n_) {
            k_ = static_cast<int>(n_);
            if (k_ % 2 != 0) k_ -= 1;
            if (k_ <= 0) k_ = 2;
        }

        double ratio = static_cast<double>(n_) / static_cast<double>(k_);
        if (ratio < 1.0) ratio = 1.0;
        L_ = static_cast<int>(ceil(log2(ratio)));

        compactors_.assign(L_ + 1, vector<Value>());
        for (int j = 0; j <= L_; ++j) {
            compactors_[j].reserve(k_);
        }

        cerr << "MRL inicializado con n = " << n_
             << ", epsilon = " << eps_
             << ", k = " << k_
             << ", L = " << L_ << "\n";
    }

    void insert(Value x) {
        compactors_[0].push_back(x);
        int j = 0;
        while (j <= L_ && (int)compactors_[j].size() >= k_) {
            compact_level(j);
            ++j;
        }
    }

    long long rank_estimate(Value x) const {
        long long ans = 0;
        for (int j = 0; j <= L_; ++j) {
            long long w = (1LL << j);
            for (auto z : compactors_[j]) {
                if (z < x) { 
                    ans += w;
                }
            }
        }
        return ans;
    }

    Value select_estimate(long long r) const {
        if (r <= 0) r = 1;
        if (r > n_) r = n_;

        struct Pair { Value z; long long w; };
        vector<Pair> B;
        B.reserve(total_size());

        for (int j = 0; j <= L_; ++j) {
            long long w = (1LL << j);
            for (auto z : compactors_[j]) {
                B.push_back({z, w});
            }
        }

        if (B.empty()) {
            throw runtime_error("Sketch vacio: no se puede hacer select()");
        }

        sort(B.begin(), B.end(), [](const Pair& a, const Pair& b) {
            if (a.z != b.z) return a.z < b.z;
            return a.w < b.w;
        });

        long long acum = 0;
        for (auto& p : B) {
            acum += p.w;
            if (acum >= r) return p.z;
        }
        return B.back().z;
    }

    Value quantile_estimate(double phi) const {
        if (phi < 0.0) phi = 0.0;
        if (phi > 1.0) phi = 1.0;
        long long r = (long long)floor(phi * (double)n_);
        if (r < 1) r = 1;
        if (r > n_) r = n_;
        return select_estimate(r);
    }

private:
    long long n_;
    double eps_;
    int k_;
    int L_;
    vector<vector<Value>> compactors_;

    size_t total_size() const {
        size_t s = 0;
        for (auto& v : compactors_) s += v.size();
        return s;
    }

    void compact_level(int j) {
        if (j > L_) return;
        auto& Aj = compactors_[j];
        sort(Aj.begin(), Aj.end());

        vector<Value> promoted;
        promoted.reserve(k_ / 2);
        // x1, x3, x5,... en notación 1-based => indices 0,2,4,... en 0-based
        for (int idx = 0; idx + 1 <= k_; idx += 2) {
            promoted.push_back(Aj[idx]);
        }
        Aj.clear();

        if (j + 1 <= L_) {
            auto& Aj1 = compactors_[j + 1];
            for (auto v : promoted) {
                Aj1.push_back(v);
            }
        }
    }
};

int main(int argc, char* argv[]) {

    if (argc != 3) {
        cerr << "Uso: ./analisis <archivo> <epsilon>\n";
        return 1;
    }

    string filename = argv[1];
    double eps = atof(argv[2]);

    // Cargar datos
    ifstream in(filename);
    if (!in) {
        cerr << "No se pudo abrir el archivo.\n";
        return 1;
    }

    vector<long long> data;
    long long x;
    while (in >> x) data.push_back(x);

    long long n = data.size();
    if (n == 0) {
        cerr << "Archivo vacío.\n";
        return 1;
    }

    // Construir el sketch MRL
    MRLSketch sketch(n, eps);
    for (long long v : data) {
        sketch.insert(v);
    }

    // Orden real del flujo
    vector<long long> sorted = data;
    sort(sorted.begin(), sorted.end());

    auto rank_real = [&](long long v) {
        return (long long)(upper_bound(sorted.begin(), sorted.end(), v) - sorted.begin());
    };

   
    //          EXPERIMENTO DE RANK
    
    int NUM_QUERIES = 10000;  
    int exitos = 0, fallos = 0;

    mt19937_64 rng(12345);
    uniform_int_distribution<int> dist(0, n - 1);

    for (int i = 0; i < NUM_QUERIES; i++) {
        long long v = data[dist(rng)];
        long long r_real = rank_real(v);
        long long r_est  = sketch.rank_estimate(v);
        long long err = llabs(r_est - r_real);

        if (err <= eps * n) exitos++;
        else fallos++;
    }

    double tasa_fallo = 100.0 * fallos / NUM_QUERIES;

    
    //EXPERIMENTO DE CUANTILES

    vector<double> quantiles = {0.25, 0.50, 0.75};
    vector<long long> realQ;
    vector<long long> estQ;

    for (double phi : quantiles) {
        long long r = floor(phi * n);
        if (r < 1) r = 1;
        if (r > n) r = n;

        long long q_real = sorted[r - 1];
        long long q_est  = sketch.quantile_estimate(phi);

        realQ.push_back(q_real);
        estQ.push_back(q_est);
    }

    
    //               RESULTADOS

    cout << "\n=====================================================\n";
    cout << " ANALISIS EXPERIMENTAL DEL SKETCH MRL\n";
    cout << " Archivo: " << filename << "\n";
    cout << " n = " << n << "\n";
    cout << " epsilon = " << eps << "\n";
    cout << "=====================================================\n\n";

    cout << "--- Rank ---\n";
    cout << "Consultas realizadas: " << NUM_QUERIES << "\n";
    cout << "Exitos: " << exitos << "\n";
    cout << "Fallos: " << fallos << "\n";
    cout << "Tasa de fallo = " << tasa_fallo << "%\n\n";

    cout << "--- Cuantiles ---\n";
    for (int i = 0; i < quantiles.size(); i++) {
        cout << "phi=" << quantiles[i]
             << "   real=" << realQ[i]
             << "   estimado=" << estQ[i] << "\n";
    }

    cout << "\nListo.\n";
    return 0;
}
