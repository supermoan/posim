#ifndef VECTORND_H
#define VECTORND_H

// https://gist.github.com/migimunz/969315

template <typename T>
struct Vector2d {
    T x, y;
    
    // constructors
    Vector2d() : x(0), y(0) {}
    
    Vector2d(T x, T y) {
        this->x = x;
        this->y = y;
    }
    
    Vector2d(const Vector2d<T> &v) {
        this->x = v.x;
        this->y = v.y;
    }
    Vector2d(const std::pair<T, T>& other) {
        this->x = other.first;
        this->y = other.second;
    }
    Vector2d(const int& other) {
        this->x = other;
        this->y = other;
    }
    // operator overloads
    Vector2d& operator=(const Vector2d<T>& v) {
        this->x = v.x;
        this->y = v.y;
        return *this;
    }
    T& operator[](int i) const { return i == 0 ? x : y; }
    T& operator()(int i) const { return (*this)[i]; }
    Vector2d<T> operator()() const { return normalize(); }
    
    bool operator==(const Vector2d<T>& v) const { return x == v.x && y == v.y; }
    bool operator!=(const Vector2d<T>& v) const { return !(*this == v); }
    bool operator==(const int& other) const { return (x == other && y == other); }
    bool operator!=(const int& other) const { return (x != other && y != other); }
    Vector2d& operator+=(T const& s) { x += s; y += s; return *this; }
    Vector2d& operator-=(T const& s) { x -= s; y -= s; return *this; }
    Vector2d& operator*=(T const& s) { x *= s; y *= s; return *this; }
    Vector2d& operator/=(T const& s) { x /= s; y /= s; return *this; }
    Vector2d& operator+=(const Vector2d& v) { x += v.x; y += v.y; return *this; }
    Vector2d& operator-=(const Vector2d& v) { x -= v.x; y -= v.y; return *this; }
    Vector2d& operator*=(const Vector2d& v) { x *= v.x; y *= v.y; return *this; }
    Vector2d& operator/=(const Vector2d& v) { x /= v.x; y /= v.y; return *this; }
    Vector2d operator+(const Vector2d<T>& v) const {
        return Vector2d<T>(x+v.x, y+v.y);
    }
    
    Vector2d operator-(const Vector2d<T>& v) const {
        return Vector2d<T>(x-v.x, y-v.y);
    }
    
    //scaling
    Vector2d operator+(T v) const { return Vector2d<T>(x+v, y+v); }
    Vector2d operator-(T v) const { return Vector2d<T>(x-v, y-v); }
    Vector2d operator*(T v) const { return Vector2d<T>(x*v, y*v); }
    Vector2d operator/(T v) const { return Vector2d<T>(x/v, y/v); }
    
    //component-wise scaling
    Vector2d operator*(const Vector2d<T> &v) const {
        return Vector2d<T>(x*v.x, y*v.y);
    }
    
    //component-wise scaling
    Vector2d operator/(const Vector2d<T> &v) const {
        return Vector2d<T>(x/v.x, y/v.y);
    }
    
    Vector2d operator-() const {
        return (Vector2d(-x, -y));
    }
    
    Vector2d operator+() const {
        return *this;
    }
    
    T length2() const {
        return x*x + y*y;
    }
    
    T length() const {
        return (T)sqrt(length2());
    }
    
    Vector2d normalize() const {
        T len(length());
        if (len < 0.000001) {
            return *this / 0.000001;
        }
        return *this / len;
    }
    
    Vector2d cw() const {
        return Vector2d(-y, +x);
    }
    
    Vector2d ccw() const {
        return Vector2d(+y, -x);
    }
    
    Vector2d rotate(T angle) const {
        angle = -angle; // rotate clockwise
        T s = sin(angle), c = cos(angle);
        return Vector2d(x * c - y * s, x * s + y * c);
    }
    
    static T dot(const Vector2d<T> &a, const Vector2d<T> &b) {
        return a.x*b.x + a.y*b.y;
    }
    static T det(const Vector2d<T> &a, const Vector2d<T> &b) {
        return a.x * b.y + a.y * b.x;
    }
    static T cross(const Vector2d<T> &a, const Vector2d<T> &b) {
        return a.x*b.y-a.y*b.x;
    }
    T angle(const Vector2d<T> &b) {
        return atan2(cross(*this, b), dot(*this, b));
    }
    float distanceFrom(const Vector2d<T> &other) {
        return sqrt(distanceFromSquared(other));
    }
    float distanceFromSquared(const Vector2d<T> &other) {
        float dx = x - other.x;
        float dy = y - other.y;
        return dx * dx + dy * dy;
    }
};

typedef Vector2d<int> Vector2di;
typedef Vector2d<float> Vector2df;
typedef Vector2d<double> Vector2dd;

#endif
