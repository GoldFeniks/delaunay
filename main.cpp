#include <chrono>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <functional>
#include "delaunay.hpp"

template<typename T>
void output_vector(std::ofstream& out, const std::vector<T>& values) {
    for (const auto& it : values)
        out << it << ' ';
    out << std::endl;
}

int main() {
    constexpr auto n = 50;
    std::vector<double> x(n), y(n);

    std::default_random_engine eng(static_cast<unsigned long>(time(nullptr)));
    std::uniform_real_distribution<double> dist(-10000, 10000);
    auto random = std::bind(dist, eng);

    for (size_t i = 0; i < n; ++i) {
        x[i] = random();
        y[i] = random();
    }

    std::ofstream out("test.txt");
    output_vector(out, x);
    output_vector(out, y);

    const auto t1 = std::chrono::system_clock::now();
    delaunay_triangulation<double> tri(x, y);
    const auto t2 = std::chrono::system_clock::now();
    std::cout << "Evaluation time(ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << std::endl;
    std::cout << "Number of points   : " << tri.points().size() << std::endl;
    std::cout << "Number of edges    : " << tri.edges().size() << std::endl;
    std::cout << "Number of triangles: " << tri.triangles().size() << std::endl;

    for (const auto& it : tri.edges())
        out << it.a.get().x << ' ' << it.a.get().y << ' '<< it.b.get().x << ' ' << it.b.get().y << std::endl;
}