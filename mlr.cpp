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
    // dejo sync_with_stdio en true para que se comporte "normal"
    // y sea más predecible con cout/cin.
    ios::sync_with_stdio(true);
    cin.tie(nullptr);

    if (argc != 4) {
        cerr << "Uso: " << argv[0] << " <n> <epsilon> <archivo>\n";
        cerr << "Ejemplo: " << argv[0] << " 15808577 0.05 chicago2015.txt\n";
        return 1;
    }

    long long n;
    double eps;
    try {
        n = stoll(argv[1]);
        eps = stod(argv[2]);
    } catch (...) {
        cerr << "Error al parsear n o epsilon.\n";
        return 1;
    }

    string filename = argv[3];
    ifstream in(filename);
    if (!in) {
        cerr << "No se pudo abrir el archivo: " << filename << "\n";
        return 1;
    }

    try {
        MRLSketch sketch(n, eps);

        long long count = 0;
        MRLSketch::Value x;
        while (count < n && (in >> x)) {
            sketch.insert(x);
            ++count;
        }

        if (count < n) {
            cerr << "Advertencia: el archivo tenia menos de n elementos ("
                 << count << " < " << n << ").\n";
        }

        cerr << "Construccion del sketch completada. Elementos leidos: "
             << count << "\n";

        cout << "Sketch MRL listo. Ingrese consultas:\n";
        cout << "  rank x        -> estimacion de rank(x)\n";
        cout << "  select r      -> estimacion de select(r)\n";
        cout << "  quantile phi  -> estimacion de quantile(phi)\n";
        cout << "  exit          -> salir\n";
        cout << "Pulse Enter para iniciar..." << endl;

        // por si acaso limpiar stdin (aunque no tocamos cin antes)
        cin.clear();
        string line;
        getline(cin, line); // Consumir el Enter inicial

        while (true) {
            cout << "> " << flush;   // esto asegura que el prompt se vea altiro
            
            if (!getline(cin, line)) break;

            if (line.empty()) continue; // Si es solo enter, vuelve a pedir input

            stringstream ss(line);
            string cmd;
            if (!(ss >> cmd)) continue; // Si la linea solo tenia espacios

            if (cmd == "exit") {
                break;
            } else if (cmd == "rank") {
                MRLSketch::Value v;
                if (!(ss >> v)) {
                    cerr << "Falta argumento para rank\n";
                    break;
                }
                long long r = sketch.rank_estimate(v);
                cout << "rank(" << v << ") ≈ " << r << endl;
            } else if (cmd == "select") {
                long long r;
                if (!(ss >> r)) {
                    cerr << "Falta argumento para select\n";
                    break;
                }
                auto v = sketch.select_estimate(r);
                cout << "select(" << r << ") ≈ " << v << endl;
            } else if (cmd == "quantile") {
                double phi;
                if (!(ss >> phi)) {
                    cerr << "Falta argumento para quantile\n";
                    break;
                }
                auto v = sketch.quantile_estimate(phi);
                cout << "quantile(" << phi << ") ≈ " << v << endl;
            
            } else {
                cout << "Comando no reconocido. Use: rank, select, quantile, exit" << endl;
            }
        }

    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}