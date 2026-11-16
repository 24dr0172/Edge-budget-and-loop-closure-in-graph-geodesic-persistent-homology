#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// Exact circular consecutive-pair ranks from unnormalised exponential gaps.
// Each side count stops before its sum reaches the central gap.
void keys(const std::vector<double>& z, std::vector<int>& sym,
          std::vector<int>& mut) {
    const int n = z.size();
    for (int i = 0; i < n; ++i) {
        int left = 0, right = 0;
        double sum = 0;
        for (int j = i - 1, k = 0; k < n - 1; ++k, --j) {
            if (j < 0) j += n;
            sum += z[j];
            if (sum >= z[i]) break;
            ++left;
        }
        sum = 0;
        for (int j = i + 1, k = 0; k < n - 1; ++k, ++j) {
            if (j >= n) j -= n;
            sum += z[j];
            if (sum >= z[i]) break;
            ++right;
        }
        sym[i] = 1 + std::min(left, right);
        mut[i] = 1 + std::max(left, right);
    }
}

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--keys") {
        int n;
        while (std::cin >> n) {
            std::vector<double> z(n); std::vector<int> a(n), b(n);
            for (auto& x : z) std::cin >> x;
            keys(z, a, b);
            for (auto x : a) std::cout << x << ' ';
            for (auto x : b) std::cout << x << ' ';
            std::cout << '\n';
        }
        return 0;
    }
    if (argc != 4) return 2;
    const int n = std::stoi(argv[1]), trials = std::stoi(argv[2]);
    const uint64_t seed = std::stoull(argv[3]);
    std::mt19937_64 rng(seed);
    std::vector<double> z(n); std::vector<int> a(n), b(n);
    int ts = 0, tm = 0, bad = 0;
    double largest_gap = 0;
    for (int rep = 0; rep < trials; ++rep) {
        double sum = 0, biggest = 0;
        for (auto& x : z) {
            // Open (0,1), using 52 random bits and midpoint conversion.
            const double u = (double(rng() >> 12) + 0.5) * 0x1.0p-52;
            x = -std::log(u); sum += x; biggest = std::max(biggest, x);
        }
        largest_gap = std::max(largest_gap, biggest / sum);
        bad += biggest / sum >= 1.0 / 3.0;
        keys(z, a, b);
        const int ma = *std::max_element(a.begin(), a.end());
        const int mb = *std::max_element(b.begin(), b.end());
        ts += std::count(a.begin(), a.end(), ma) > 1;
        tm += std::count(b.begin(), b.end(), mb) > 1;
    }
    std::cout << std::setprecision(17) << n << ',' << trials << ',' << seed
              << ',' << ts << ',' << tm << ',' << bad << ',' << largest_gap << '\n';
}
