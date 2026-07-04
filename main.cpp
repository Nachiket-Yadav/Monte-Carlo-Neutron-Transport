#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

enum class Status { Travelling, Transmitted, Reflected, Absorbed };

class Particle {
    private:
        double x;   // position, in mean free paths from the left face
        double mu;  // direction cosine: +1 straight ahead, -1 straight back
 
    public:
        Status status;
 
        Particle() : x(0.0), mu(1.0), status(Status::Travelling) {} // Initialize to the left edge of the slab, moving right.
 
        void move(double s, double L) {

            x += s * mu;
            if (x > L)      status = Status::Transmitted;
            else if (x < 0) status = Status::Reflected;
            // otherwise still Travelling: the next collision happens right here

        }
 
        void scatter(double u) { 

            mu = 2.0 * u - 1.0; // isotropic: mu uniform on (-1, 1)

        }  

        void absorb() { 

            status = Status::Absorbed; 

        }
};

class Material {

    private:

        double L;  // slab thickness, in mean free paths
        double c;  // scattering probability per collision (Sigma_s / Sigma_t)
 
    public:

        Material(double length, double scatterProb) : L(length), c(scatterProb) {
            if (c < 0.0 || c > 1.0) {
                throw std::invalid_argument("Scattering probability must be between 0 and 1.");
            }
            if (L <= 0.0) {
                throw std::invalid_argument("Material length must be positive.");
            }
        }
 
        double getLength() const { return L; }
        double getScatteringProbability() const { return c; }
        
};

// One clean roll in (0, 1): re-draw the (astronomically rare) exact 0
// so -log(u) can never blow up.
double getRoll(std::mt19937& engine) {

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(engine);
    while (u <= 0.0) {
        u = uniform(engine);
    }
    return u;

}

// Run one particle's history and return end state.
Status simulateOne(const Material& material, std::mt19937& engine) {

    Particle particle;
 
    while (particle.status == Status::Travelling) {
        double s = -std::log(getRoll(engine));       // exponential free path
        particle.move(s, material.getLength());
 
        if (particle.status == Status::Travelling) { // still inside -> a collision happens
            double u = getRoll(engine);
            if (u >= material.getScatteringProbability()) {
                particle.absorb();
            } else {
                particle.scatter(getRoll(engine));   // fresh roll for the new direction
            }
        }
    }
    return particle.status;

}

struct Tally {

    long long transmitted = 0;
    long long reflected   = 0;
    long long absorbed    = 0;
 
    long long total() const { return transmitted + reflected + absorbed; }

};

// Run N independent histories and count the end states. The switch has no
// default on purpose: if Status ever grows a new value, -Wall warns here,
// and a particle that somehow finishes while Travelling throws loudly
// instead of evaporating uncounted.
Tally runBatch(const Material& material, long long N, std::mt19937& engine) {

    Tally tally;
    for (long long i = 0; i < N; ++i) {
        switch (simulateOne(material, engine)) {
            case Status::Transmitted: ++tally.transmitted; break;
            case Status::Reflected:   ++tally.reflected;   break;
            case Status::Absorbed:    ++tally.absorbed;    break;
            case Status::Travelling:
                throw std::logic_error("particle finished while still Travelling");
        }
    }
    if (tally.total() != N) {
        throw std::logic_error("conservation broken: T + R + A != N");
    }
    return tally;

}

void writeCsvHeader(std::ofstream& out) {

    out << "c,L,N,seed,transmitted,reflected,absorbed,"
           "T_frac,R_frac,A_frac,T_err,R_err,A_err\n";

}

// Each row is one full experiment: the settings that produced it, the raw
// counts, the fractions, and their binomial standard errors sqrt(p(1-p)/N)
// -- the error bars for the Python plots.
void writeCsvRow(std::ofstream& out, double c, double L, long long N,
                 unsigned seed, const Tally& tally) {

    // TRA fractions
    double T = static_cast<double>(tally.transmitted) / N;
    double R = static_cast<double>(tally.reflected)   / N;
    double A = static_cast<double>(tally.absorbed)    / N;

    double Terr = std::sqrt(T * (1.0 - T) / N);
    double Rerr = std::sqrt(R * (1.0 - R) / N);
    double Aerr = std::sqrt(A * (1.0 - A) / N);
 
    out << c << ',' << L << ',' << N << ',' << seed << ','
        << tally.transmitted << ',' << tally.reflected << ',' << tally.absorbed << ','
        << T << ',' << R << ',' << A << ','
        << Terr << ',' << Rerr << ',' << Aerr << '\n';

}

int main() {

    const long long N = 1'000'000;  // histories per configuration
    const unsigned seed = 42;       // To make every run reproducible
    std::mt19937 engine(seed); // Mersenne Twister engine, generates random numbers
 
    // Validations. 
    // The c = 0 rows trace out the Beer-Lambert line - the case of the pure absorber.
    // T = e^{-L}; c = 0.5 exercises every branch.
    // c = 1 is the pure scatterer, which never absorbs.
    const std::vector<std::pair<double, double>> configs = {  // {c, L}
        {0.0, 1.0}, {0.0, 2.0}, {0.0, 3.0}, {0.0, 4.0}, {0.0, 5.0},
        {0.5, 5.0},
        {1.0, 5.0},
    };
 
    std::ofstream out("results.csv");
    if (!out) {
        std::cerr << "Could not open results.csv for writing.\n";
        return 1;
    }
    out << std::setprecision(9);
    writeCsvHeader(out);
 
    std::cout << std::setprecision(6);
    bool allPass = true;
 
    for (const auto& [c, L] : configs) {
        Material material(L, c);  // constructor takes (length, scatterProb)
        Tally tally = runBatch(material, N, engine);
        writeCsvRow(out, c, L, N, seed, tally);
        
        // TRA fractions/probabilities
        double T = static_cast<double>(tally.transmitted) / N;
        double R = static_cast<double>(tally.reflected)   / N;
        double A = static_cast<double>(tally.absorbed)    / N;
        double Terr = std::sqrt(T * (1.0 - T) / N);
 
        std::cout << "c = " << c << "  L = " << L
                  << "  |  T = " << T << "  R = " << R << "  A = " << A;
 
        if (c == 0.0) {
            // Pre-registered expectation: T within 3 standard errors of e^{-L},
            // and reflection impossible without a scatter.
            double exact = std::exp(-L);
            bool pass = std::fabs(T - exact) < 3.0 * Terr && tally.reflected == 0;
            std::cout << "  |  exact e^-L = " << exact
                      << (pass ? "  [PASS]" : "  [FAIL]");
            if (!pass) allPass = false;
        } else if (c == 1.0) {
            // Pre-registered expectation: a pure scatterer never absorbs.
            bool pass = (tally.absorbed == 0);
            std::cout << "  |  expect A = 0"
                      << (pass ? "  [PASS]" : "  [FAIL]");
            if (!pass) allPass = false;
        }
        std::cout << '\n';
    }
 
    std::cout << (allPass ? "All checks passed." : "SOME CHECKS FAILED.")
              << "  Wrote results.csv\n";
    return allPass ? 0 : 1;

}