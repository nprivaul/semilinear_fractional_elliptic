#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {
constexpr double PI = 3.14159265358979323846;

// One global RNG (thread_local if you’ll multithread later)
inline std::mt19937_64& rng() {
    static thread_local std::mt19937_64 eng{123456789ULL}; // set fixed seed or std::random_device{}()
    return eng;
}

// Samplers using the shared engine
inline double normal(double mean, double stdev) {
    std::normal_distribution<double> dist(mean, stdev);
    return dist(rng());
}
inline double uniform(double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng());
}
inline double exponential(double lambda) {
    std::exponential_distribution<double> dist(lambda);
    return dist(rng());
}
inline double gamma_sample(double alpha, double beta) {
    std::gamma_distribution<double> dist(alpha, beta);
    return dist(rng());
}

// Hypergeometric 2F1 — simple series (unchanged logic)
double hyperg_2F1(double a, double b, double c, double x) {
    const double TOL = 1e-4; // was 10e-5 (i.e., 1e-4). Keep explicit.
    double term = (a * b * x) / c;
    double value = 1.0 + term;
    int n = 1;

    while (std::abs(term) > TOL) {
        ++a; ++b; ++c; ++n;
        term *= (a * b * x) / (c * n);
        value += term;
    }
    return value;
}

double gamma_pdf(double x, double a) {
    return std::pow(x, a - 1) * std::exp(-x) / std::tgamma(a);
}

} // namespace

// Original globals kept (if you truly need them global)
double k00 = 0.0;
int d = 1;

// ---- Helpers on vectors ----
double norm2(const std::vector<double>& x) {
    double r2 = 0.0;
    for (double xi : x) r2 += xi * xi;
    return std::sqrt(r2);
}

// Brownian step: y = x + sqrt(t) * N(0, I)
std::vector<double> brownian(const std::vector<double>& x, double t) {
    const int n = static_cast<int>(x.size());
    std::vector<double> y(n);
    const double s = std::sqrt(t);
    for (int i = 0; i < n; ++i) y[i] = x[i] + s * normal(0.0, 1.0);
    return y;
}

// --- Your distributions built from standard RNGs ---
double law_uniform(double a, double b) { return uniform(a, b); }
double law_exponential(double lambda) { return exponential(lambda); }
double law_gamma(double alpha, double beta) { return gamma_sample(alpha, beta); }

// Subordinator (kept, but uses std RNG now)
double subordinator(double x, double a, double /*e*/, double t) {
    const double lambda0 = -PI / 2.0;
    const double W = -std::log(uniform(0.0, 1.0));
    const double U = law_uniform(-PI / 2.0, PI / 2.0);
    const double A1 = std::sin((a / 2.0) * (U - lambda0)) / std::pow(std::cos(U), 2.0 / a);
    const double A2 = std::pow(std::cos(U - (a / 2.0) * (U - lambda0)) / W, (2.0 - a) / a);
    return x + std::pow(t, 2.0 / a) * A1 * A2;
}

double law_stable(double x, double a, double /*b*/, double t) {
    const double W = law_exponential(1.0);
    const double U = law_uniform(-PI / 2.0, PI / 2.0);
    const double A1 = std::sin(a * (U - 0.0)) / std::pow(std::cos(U), 1.0 / a);
    const double A2 = std::pow(std::cos(U - a * (U - 0.0)) / W, (1.0 - a) / a);
    double s = x + std::pow(t, 1.0 / a) * A1 * A2;
    if (s == 0.0) std::printf("%.80f\n", s);
    return s;
}

// Problem-specific pieces
double condition_initiale(double a, const std::vector<double>& x) {
    const double r = norm2(x);
    if (r > 1.0) return 0.0;
    const double r2 = r * r;
    return std::pow(1.0 - r2, k00 + a / 2.0);
}

double constante(double a) {
    return -std::pow(2.0, a)
           * std::tgamma((d + a) / 2.0)
           * std::tgamma(k00 + 1 + a / 2.0)
           / std::tgamma(d / 2.0)
           / std::tgamma(k00 + 1);
}

double constantetilde(double a) {
    return -std::pow(2.0, a)
           * std::tgamma((d + a) / 2.0)
           * std::tgamma(k00 + 1 + a / 2.0)
           / std::tgamma(-a / 2.0)
           / std::tgamma(k00 + 1 + (d + a) / 2.0);
}

double sol_funct(double a, const std::vector<double>& x) {
    const double r = norm2(x);
    const double r2 = r * r;
    if (r > 1.0) {
        return -constantetilde(a)
               * hyperg_2F1((d + a) / 2.0, (2 + a) / 2.0, k00 + 1 + (d + a) / 2.0, 1.0 / r2)
               * std::pow(r, -a - d);
    } else {
        return -constante(a)
               * hyperg_2F1((d + a) / 2.0, -k00, d / 2.0, r2);
    }
}

// One MC sample
double solution_basique_elliptique_sample(double mu, double a, double e, const std::vector<double>& x) {
    const double r = norm2(x);
    if (r >= 1.0) return 0.0;

    const double W = law_gamma(mu, 1.0);
    std::vector<double> y = brownian(x, subordinator(0.0, a, e, std::pow(2.0, a / 2.0) * (W / 500.0)));
    if (norm2(y) >= 1.0) return 0.0;

    for (int i = 0; i < 500; ++i) {
        y = brownian(y, subordinator(0.0, a, e, std::pow(2.0, a / 2.0) * (W / 500.0)));
        if (norm2(y) >= 1.0) return 0.0;
    }
    return sol_funct(a, y) / gamma_pdf(W, mu);
}

// MC estimator
double solution_basique_elliptique(double mu, double a, double e, const std::vector<double>& x, int n) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i) acc += solution_basique_elliptique_sample(mu, a, e, x);
    const double avg = acc / static_cast<double>(n);
    std::printf("result=%f\n", avg);
    return avg;
}

int main() {
    // Fixed seed for reproducibility; comment this out to randomize:
    rng().seed(0xC0FFEEULL);

    double alpha = 0.8;
    double mu = 0.5 * (2 - 2 / 1.5);
    double e = 0.0;
    const double constant_factor = std::pow(2.0, -alpha) * std::tgamma(0.5) /
                                   (std::tgamma(1 + alpha / 2.0) * std::tgamma((1 + alpha) / 2.0));

    std::vector<double> x(d);

    std::FILE* mc_estimates = std::fopen("mc_estimates", "w");
    if (!mc_estimates) {
        std::perror("mc_estimates");
        return 1;
    }

    clock_t c1, c2;

    for (int k = 0; k <= 10; ++k) {
        x[0] = -1.0 + k * 0.1;
        for (int j = 1; j < d; ++j) x[j] = 0.0;

        c1 = std::clock();

        std::time_t rawtime;
        std::time(&rawtime);
        std::tm* timeinfo = std::localtime(&rawtime);
        std::printf("Current time: %s", std::asctime(timeinfo));

        std::printf("x=%f\n", x[0]);
        int N = 1000000;
	double est = solution_basique_elliptique(mu, alpha, e, x, N);
        std::printf("y=%f\n", est);

        c2 = std::clock();
        std::printf("x\t\t y\n");

        const double r = norm2(x);
        const double r2 = r * r;
        const double exact_inner = std::pow(1.0 - std::min(r2, 1.0), k00 + alpha / 2.0);
        const double est_scaled = constant_factor * est;
        const double exact_scaled = constant_factor * exact_inner;

        std::printf("%f\t%f\t %f \t Exact=%f\n", x[0], r, est_scaled, exact_scaled);
        std::fprintf(mc_estimates, "%f\t%f\t%f\n", x[0], est_scaled, exact_scaled);
        std::fprintf(mc_estimates, "%f\t%f\t%f\n", -x[0], est_scaled, exact_scaled);

        std::printf("Time taken: %g minutes\n\n", (c2 - c1) / 60.0 / static_cast<double>(CLOCKS_PER_SEC));
    }

    std::fclose(mc_estimates);

    std::FILE* exact_solution = std::fopen("exact_solution", "w");
    if (!exact_solution) {
        std::perror("exact_solution");
        return 1;
    }

    for (int k = 0; k <= 4000; ++k) {
        x[0] = -2.0 + k * 0.001;
        for (int j = 1; j < d; ++j) x[j] = 0.0;

        c1 = std::clock();
        std::time_t rawtime;
        std::time(&rawtime);
        std::tm* timeinfo = std::localtime(&rawtime);
        std::printf("Current time: %s", std::asctime(timeinfo));
        std::printf("x=%f\n", x[0]);
        c2 = std::clock();

        const double r = norm2(x);
        const double r2 = r * r;
        const double exact_inner = std::pow(1.0 - std::min(r2, 1.0), k00 + alpha / 2.0);
        const double exact_scaled = constant_factor * exact_inner;

        std::fprintf(exact_solution, "%f\t%f\n", x[0], exact_scaled);

        std::printf("Time taken: %g minutes\n\n", (c2 - c1) / 60.0 / static_cast<double>(CLOCKS_PER_SEC));
    }

    std::fclose(exact_solution);
    system("gnuplot script.plt");
    return 0;
}
