#pragma once
#include <list>
#include <tuple>
#include <cmath>
#include <vector>
#include <cstddef>
#include <utility>
#include <algorithm>

template<typename T>
class delaunay_triangulation {

private:

    class edge;
    class point;
    class triangle;

    template<typename V>
    class data_value_wrapper;

public:

    template<typename VX, typename VY>
    delaunay_triangulation(const VX& x, const VY& y) : _first_triangle(&_triangles, 0) {
        _points.reserve(x.size());
        for (size_t i = 0; i < x.size(); ++i)
            _points.emplace_back(x[i], y[i]);
        triangulate();
    }

    delaunay_triangulation(const delaunay_triangulation& other) : _first_triangle(&_triangles, 0) {
        *this = other;
    }

    delaunay_triangulation(delaunay_triangulation&& other) : _first_triangle(&_triangles, 0) {
        *this = std::move(other);
    }

    auto& operator=(const delaunay_triangulation& other) {
        _edges = other._edges;
        _points = other._points;
        _triangles = other._triangles;
        _reassign_data();
        return *this;
    }

    auto& operator=(delaunay_triangulation&& other) noexcept {
        _edges = std::move(other._edges);
        _points = std::move(other._points);
        _triangles = std::move(other._triangles);
        _reassign_data();
        return *this;
    }

    const auto& edges() const {
        return _edges;
    }

    const auto& points() const {
        return _points;
    }

    const auto& triangles() const {
        return _triangles;
    }

    const auto find_triangle(const T& x, const T& y) const {
        return find_triangle({x, y});
    }

    const auto find_triangle(const point& p) const {
        const auto b = _first_triangle._data == &_triangles;
        return find_triangle(_first_triangle, p);
    }

    const auto find_triangle(const data_value_wrapper<triangle>& t, const T& x, const T& y) const {
        return _find_triangle(t, x, y);
    }

    const auto find_triangle(const data_value_wrapper<triangle>& t, const point& p) const {
        return _find_triangle(t, p);
    }

private:

    const data_value_wrapper<triangle> _first_triangle;

    std::vector<edge> _edges;
    std::vector<point> _points;
    std::vector<triangle> _triangles;

    void _reassign_data() {
        for (auto& it : _edges) {
            it.a._data = &_points;
            it.b._data = &_points;
            it.t1._data = &_triangles;
            if (it.t2.is_valid())
                it.t2._data = &_triangles;
        }
        for (auto& it : _triangles) {
            it.a._data = &_edges;
            it.b._data = &_edges;
            it.c._data = &_edges;
        }
    }

    static constexpr auto eps = T(1e-10);

    template<typename V>
    class data_value_wrapper {

    public:

        data_value_wrapper() = default;

        data_value_wrapper(std::vector<V>* data, const int64_t index) : _data(data), _index(index) {}

        bool operator==(const data_value_wrapper<V>& other) const {
            return _index == other._index;
        }

        bool operator!=(const data_value_wrapper<V>& other) const {
            return !(*this == other);
        }

        bool is_valid() const {
            return _data != nullptr;
        }

        const auto& get() const {
            return (*_data)[_index];
        }

        auto& get() {
            return (*_data)[_index];
        }

        const auto& index() {
            return _index;
        }

    private:

        int64_t _index = -1;
        std::vector<V>* _data = nullptr;

        friend class delaunay_triangulation;

    };

    struct point {

        T x, y;

        point() = default;
        point(T x, T y) : x(std::move(x)), y(std::move(y)) {}

    };

    struct edge {

        data_value_wrapper<point> a, b;
        data_value_wrapper<triangle> t1, t2;

        edge(data_value_wrapper<point> a, data_value_wrapper<point> b) : a(std::move(a)), b(std::move(b)) {}

        void flip() {
            std::swap(a, b);
        }

    };

    class triangle {

    public:

        data_value_wrapper<edge> a, b, c;

        triangle(data_value_wrapper<edge> a, data_value_wrapper<edge> b, data_value_wrapper<edge> c) :
            a(std::move(a)), b(std::move(b)), c(std::move(c)) {}

        void rearrange(const data_value_wrapper<edge>& e) {
            if (e == b)
                std::swap(a, b);
            else if (e == c)
                std::swap(a, c);
            _flip_edges(e.get().b, b, c);
        }

        auto barycentric_coordinates(const point& p) const {
            const auto& [a, b, c] = get_points();
            const auto v0x = b.x - a.x;
            const auto v0y = b.y - a.y;
            const auto v1x = c.x - a.x;
            const auto v1y = c.y - a.y;
            const auto v2x = p.x - a.x;
            const auto v2y = p.y - a.y;
            const auto d00 = v0x * v0x + v0y * v0y;
            const auto d01 = v0x * v1x + v0y * v1y;
            const auto d11 = v1x * v1x + v1y * v1y;
            const auto d20 = v2x * v0x + v2y * v0y;
            const auto d21 = v2x * v1x + v2y * v1y;
            const auto den = d00 * d11 - d01 * d01;
            const auto v   = (d11 * d20 - d01 * d21) / den;
            const auto w   = (d00 * d21 - d01 * d20) / den;
            return std::make_tuple(T(1) - v - w, v, w);
        }

        auto points() const {
            const auto& p0 = a.get().a;
            const auto& p1 = a.get().b;
            const auto& p2 = b.get().a;
            if (p0 != p2 && p1 != p2)
                return std::tie(p0, p1, p2);
            return std::tie(p0, p1, b.get().b);
        }

        auto get_points() const {
            const auto& [a, b, c] = points();
            return std::tie(a.get(), b.get(), c.get());
        }

        bool contains_point(const point& p) const {
            const auto& [a, b, c] = get_points();
            const auto d1 = _sign(p, a, b);
            const auto d2 = _sign(p, b, c);
            const auto d3 = _sign(p, c, a);

            const auto neg = (d1 < -eps) || (d2 < -eps) || (d3 < -eps);
            const auto pos = (d1 > eps) || (d2 > eps) || (d3 > eps);

            return !(neg && pos);
        }

    private:

        static void _flip_edges(const data_value_wrapper<point>& p, data_value_wrapper<edge>& a, data_value_wrapper<edge>& b) {
            auto& ag = a.get();
            auto& bg = b.get();
            if (p == ag.b)
                ag.flip();
            else if (p == bg.a)
                std::swap(a, b);
            else if (p == bg.b) {
                bg.flip();
                std::swap(a, b);
            }
            if (a.get().b != b.get().a)
                b.get().flip();
        }

        static auto _sign(const point& a, const point& b, const point& c) {
            return (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y);
        }

    };

    void triangulate() {
        const auto p0 = data_value_wrapper(&_points, 0);
        const auto& p0g = p0.get();

        std::vector<std::pair<T, data_value_wrapper<point>>> points;
        points.reserve(_points.size() - 1);

        for (size_t i = 0; i < _points.size() - 1; ++i)
            points.emplace_back(_norm(p0g, _points[i + 1]), data_value_wrapper(&_points, i + 1));
        std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

        auto p1 = points[0].second;
        const auto& p1g = p1.get();
        auto p2 = points[1].second;
        size_t di = 1;
        auto [r, c] = _circle_center(p0g, p1g, p2.get());
        for (size_t i = 2; i < points.size(); ++i) {
            const auto [cr, cc] = _circle_center(p0g, p1g, points[i].second.get());
            if (cr < r) {
                r = cr;
                c = cc;
                p2 = points[i].second;
                di = i;
            }
        }
        points.erase(points.begin() + di);

        if (_is_clockwise(p0g, p1g, p2.get()))
            std::swap(p1, p2);

        for (size_t i = 1; i < points.size(); ++i)
            points[i].first = _norm(c, points[i].second.get());
        std::sort(points.begin() + 1, points.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

        auto e0 = _add_edge(p0, p1);
        auto e1 = _add_edge(p1, p2);
        auto e2 = _add_edge(p2, p0);
        _add_triangle(e0, e1, e2);

        std::list<data_value_wrapper<edge>> hull;
        hull.push_back(e0);
        hull.push_back(e1);
        hull.push_back(e2);

        for (size_t i = 1; i < points.size(); ++i) {
            const auto& p = points[i].second;
            auto& e = hull.front();
            const auto buff = !_is_clockwise(e.get().a.get(), e.get().b.get(), p.get());
            auto a = _find_first(++hull.begin(), p.get(), buff);
            auto b = _find_first(hull.rbegin(), p.get(), buff).base();
            if (!buff) {
                if (b == hull.end()) {
                    auto t = _add_triangle(hull.front(), p);
                    e0 = t.get().a;
                    e1 = t.get().b;
                    _add_triangles(e1, ++hull.begin(), a);
                } else {
                    auto t = _add_triangle(*b, p);
                    e0 = t.get().a;
                    e1 = t.get().b;
                    auto d = b;
                    _add_triangles(e1, ++d, hull.end());
                    _add_triangles(e1, hull.begin(), a);
                    hull.erase(b, hull.end());
                }
                a = hull.erase(hull.begin(), a);
            } else {
                auto t = _add_triangle(*a, p);
                e0 = t.get().a;
                e1 = t.get().b;
                auto d = a;
                _add_triangles(e1, ++d, b);
                a = hull.erase(a, b);
            }
            hull.insert(a, {e0, e1});
        }

        _flip_triangles();
    }

    static auto _norm(const point& a, const point& b) {
        return std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2);
    }

    static auto _circle_center(const point& a, const point& b, const point& c) {
        const auto& [x0, y0] = a;
        const auto& [x1, y1] = b;
        const auto& [x2, y2] = c;
        const auto y01 = y0 - y1;
        const auto y02 = y0 - y2;
        const auto y12 = y1 - y2;
        const auto sx0 = std::pow(x0, 2);
        const auto sx1 = std::pow(x1, 2);
        const auto sx2 = std::pow(x2, 2);
        const auto sy0 = std::pow(y0, 2);
        const auto sy1 = std::pow(y1, 2);
        const auto sy2 = std::pow(y2, 2);
        const auto x = (sx2 * -y01 + sx1 * y02 - (sx0 + y01 * y02) * y12) / (T(2) * (x2 * -y01 + x1 * y02 + x0 * -y12));
        const auto y = (-sx1 * x2 + sx0 * (x2 - x1) + x2 * y01 * (y0 + y1) + x0 * (sx1 - sx2 + sy1 - sy2) + x1 * (sx2 - sy0 + sy2)) /
                       (T(2) * (x2 * y01 + x0 * y12 + x1 * -y02));
        const auto p = point(x, y);
        return std::make_tuple(_norm(a, p), p);
    }

    static auto _is_clockwise(const point& a, const point& b, const point& c) {
        const auto& [x0, y0] = a;
        const auto& [x1, y1] = b;
        const auto& [x2, y2] = c;
        return x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1) < eps;
    }

    auto _add_edge(const data_value_wrapper<point>& a, const data_value_wrapper<point>& b) {
        _edges.emplace_back(a, b);
        return data_value_wrapper(&_edges, _edges.size() - 1);
    }

    auto _add_triangle(data_value_wrapper<edge>& a, data_value_wrapper<edge>& b, data_value_wrapper<edge>& c) {
        _triangles.emplace_back(a, b, c);
        auto res = data_value_wrapper(&_triangles, _triangles.size() - 1);
        a.get().t1 = res;
        b.get().t1 = res;
        c.get().t1 = res;
        return res;
    }

    auto _add_triangle(data_value_wrapper<edge>& c, const data_value_wrapper<point>& p) {
        auto a = _add_edge(c.get().a, p);
        auto b = _add_edge(p, c.get().b);
        c.get().flip();
        _triangles.emplace_back(a, b, c);
        auto res = data_value_wrapper(&_triangles, _triangles.size() - 1);
        a.get().t1 = res;
        b.get().t1 = res;
        c.get().t2 = res;
        return res;
    }

    auto _add_triangle(data_value_wrapper<edge>& a, data_value_wrapper<edge>& c) {
        auto b = _add_edge(a.get().a, c.get().b);
        c.get().flip();
        _triangles.emplace_back(a, b, c);
        auto res = data_value_wrapper(&_triangles, _triangles.size() - 1);
        a.get().t2 = res;
        b.get().t1 = res;
        c.get().t2 = res;
        return res;
    }

    template<typename It>
    static auto _find_first(It it, const point& p, const bool val) {
        while (true) {
            const auto& e = it->get();
            if (_is_clockwise(e.a.get(), e.b.get(), p) == val)
                break;
            ++it;
        }
        return it;
    }

    template<typename It>
    void _add_triangles(data_value_wrapper<edge>& e, It it, const It end) {
        while (it != end) {
            auto t = _add_triangle(e, *it);
            e = t.get().b;
            ++it;
        }
    }

    auto _flip_triangles() {
        std::vector<data_value_wrapper<edge>> edges;
        edges.reserve(_edges.size());
        for (size_t i = 0; i < _edges.size(); ++i)
            if (_edges[i].t2.is_valid())
                edges.emplace_back(data_value_wrapper(&_edges, i));
        bool do_flip = true;
        while (do_flip) {
            do_flip = false;
            for (auto& it : edges)
                do_flip |= _flip_triangles_over_edge(it);
        }
    }

    static bool _flip_triangles_over_edge(data_value_wrapper<edge>& e) {
        auto& eg = e.get();
        auto& t1 = eg.t1.get();
        auto& t2 = eg.t2.get();
        t1.rearrange(e);
        t2.rearrange(e);
        if (_should_flip(t1.b.get().b, eg.b, eg.a, t2.b.get().b)) {
            e.get().a = t2.b.get().b;
            e.get().b = t1.b.get().b;
            std::swap(t1.b, t2.c);
            _reassign_triangle(t1.b, eg.t2, eg.t1);
            _reassign_triangle(t2.c, eg.t1, eg.t2);
            return true;
        }
        return false;
    }

    static void _reassign_triangle(data_value_wrapper<edge>& e, data_value_wrapper<triangle>& t1, data_value_wrapper<triangle>& t2) {
        if (e.get().t1 == t1)
            e.get().t1 = t2;
        else
            e.get().t2 = t2;
    }

    static bool _should_flip(const data_value_wrapper<point>& a, const data_value_wrapper<point>& b,
                             const data_value_wrapper<point>& c, const data_value_wrapper<point>& d) {
        const auto& [xa, ya] = a.get();
        const auto& [xb, yb] = b.get();
        const auto& [xc, yc] = c.get();
        const auto& [xd, yd] = d.get();

        const auto xba = xb - xa;
        const auto yba = yb - ya;
        const auto xca = xc - xa;
        const auto yca = yc - ya;
        const auto xbd = xb - xd;
        const auto ybd = yb - yd;
        const auto xcd = xc - xd;
        const auto ycd = yc - yd;

        const auto cosa = xba * xca + yba * yca;
        const auto cosb = xbd * xcd + ybd * ycd;

        if (cosa < -eps && cosb < -eps)
            return true;

        if (cosa > eps && cosb > eps)
            return false;

        const auto sina = std::abs(xba * yca - yba * xca);
        const auto sinb = std::abs(xbd * ycd - ybd * xcd);

        if (cosa * sinb + sina * cosb < -eps)
            return true;
        return false;
    }

    static const auto _find_triangle(const data_value_wrapper<triangle>& wt, const point& p) {
        const auto& t = wt.get();
        if (t.contains_point(p))
            return wt;
        const auto& rp = t.a.get().a;
        point cp;
        if (rp == t.b.get().a || rp == t.b.get().b)
            cp = _half_way_point(rp.get(), t.c.get());
        else
            cp = _half_way_point(rp.get(), t.b.get());
        if (_does_intersect(cp, p, t.a.get()))
            return _find_triangle(t.a.get(), wt, p);
        if (_does_intersect(cp, p, t.b.get()))
            return _find_triangle(t.b.get(), wt, p);
        return _find_triangle(t.c.get(), wt, p);
    }

    static const auto _find_triangle(const edge& e, const data_value_wrapper<triangle>& wt, const point& p) {
        if (e.t2.is_valid() && wt == e.t1)
            return _find_triangle(e.t2, p);
        return wt;
    }

    static auto _half_way_point(const point& p, const edge& e) {
        const auto& a = e.a.get();
        const auto& b = e.b.get();
        const auto  x = (b.x + a.x) / T(2);
        const auto  y = (b.y + a.y) / T(2);
        return point((x + p.x) / T(2), (y + p.y) / T(2));
    }

    static bool _does_intersect(const point& a, const point& b, const edge& ed) {
        const auto& [ax, ay] = a;
        const auto& [bx, by] = b;
        const auto& [cx, cy] = ed.a.get();
        const auto& [dx, dy] = ed.b.get();
        const auto e = _det(bx - ax, by - ay, dx - cx, dy - cy);
        if (std::abs(e) < eps)
            return false;
        const auto f = _det(dy - cy, dx - cx, cy - ay, cx - ax) / e;
        const auto g = _det(by - ay, bx - ax, cy - ay, cx - ax) / e;
        return f > -eps && f < T(1) + eps && g > -eps && g < T(1) + eps;
    }

    static auto _det(const T& a, const T& b, const T& c, const T& d) {
        return a * d - b * c;
    }

};
