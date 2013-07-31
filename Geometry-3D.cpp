#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <queue>
#include <stack>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cassert>
#include <algorithm>
using namespace std;

#define T double
#define CP const Point &
#define VP vector<Point>
#define CP3 const Point3 &
#define VP3 vector<Point3>
#define VI vector<int>
#define VVI vector<VI>


const double eps = 1e-8, pi = atan2(0, -1);
int sgn(double a) {
    return a < -eps ? -1 : a < eps ? 0 : 1;
}

struct Point {
    T x, y;
    Point (T xx = 0, T yy = 0) : x(xx), y(yy){}
    Point operator + (CP b) const {
        return Point(x + b.x, y + b.y);
    }
    Point operator - (CP b) const {
        return Point(x - b.x, y - b.y);
    }
    Point operator * (const T & k) const {
        return Point(k * x, k * y);
    }
    bool operator == (CP b) const {
        return fabs(x - b.x) < eps && fabs(y - b.y) < eps;
    }
    bool operator != (CP b) const {
        return !((*this) == b);
    }

    double abs() const {
        return sqrt(x * x + y * y);
    }
    T abs2() const {
        return x * x + y * y;
    }
    T dot(CP b) const {
        return x * b.x + y * b.y;
    }
    T det(CP b) const {
        return x * b.y - y * b.x;
    }
    Point normalized() const {
        //Be careful with zero vectors
        return Point(x / abs(), y / abs());
    }
    Point rotate(double theta) const {
        return Point(x * cos(theta) - y * sin(theta),
            x * sin(theta) + y * cos(theta));
    }
};
typedef Point Vector;

struct Point3 {
    T x, y, z;
    Point3 (T xx = 0, T yy = 0, T zz = 0) : x(xx), y(yy), z(zz){}
    Point3 operator + (CP3 b) const {
        return Point3(x + b.x, y + b.y, z + b.z);
    }
    Point3 operator - (CP3 b) const {
        return Point3(x - b.x, y - b.y, z - b.z);
    }
    Point3 operator * (const T & k) const {
        return Point3(k * x, k * y, k * z);
    }
    bool operator == (CP3 b) const {
        return fabs(x - b.x) < eps && fabs(y - b.y) < eps && fabs(z - b.z) < eps;
    }
    bool operator != (CP3 b) const {
        return !((*this) == b);
    }

    double abs() const {
        return sqrt(x * x + y * y + z * z);
    }
    T abs2() const {
        return x * x + y * y + z * z;
    }
    T dot(CP3 b) const {
        return x * b.x + y * b.y + z * b.z;
    }
    Point3 det(CP3 b) const {
        return Point3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
    Point3 normalized() const {
        //Be careful with zero vectors
        return Point3(x / abs(), y / abs(), z / abs());
    }
    Point3 rotate(Point3 axis, double theta) const {
        Point3 xaxis, yaxis, zaxis, ori;
        zaxis = axis.normalized();
        ori = zaxis * ((*this).dot(zaxis));
        yaxis = ((*this) - ori).normalized();
        xaxis = yaxis.det(zaxis);
        Point proj = Point(((*this) - ori).dot(xaxis),
            ((*this) - ori).dot(yaxis));
        proj = proj.rotate(theta);
        return xaxis * proj.x + yaxis* proj.y + ori;
    }
};
typedef Point3 Vector3;

namespace ConvexHull {
    struct cmp {
        bool operator () (CP p1, CP p2) const {
            if (sgn(p1.x - p2.x) == 0) {
                return sgn(p1.y - p2.y) < 1;
            } else {
                return sgn(p1.x - p2.x) < 1;
            }
        }
    };
    VP solve(VP ps) {
        int n = ps.size(), k = 0;
        if (n <= 1) {
            return ps;
        }
        sort(ps.begin(), ps.end(), cmp());

        VP qs(n * 2);
        for (int i = 0; i < n; qs[k++] = ps[i++]) {
            while(k > 1 && (qs[k - 1] - qs[k - 2]).det(ps[i] - qs[k - 1]) < eps) {
                k--;
            }
        }
        for (int i = n - 2, t = k; i >= 0; qs[k++] = ps[i--]) {
            while (k > t && (qs[k - 1] - qs[k - 2]).det(ps[i] - qs[k - 1]) < eps) {
                k--;
            }
        }

        qs.resize(k - 1);
        return qs;
    }
}

namespace ConvexHull3 {
    //Before unique() there should be sort()
    struct face {
        int index[3];
        int & operator [] (const int n) {
            return index[n];
        }
    };
    double randeps() {
        return (rand() / (double)RAND_MAX - 0.5) * eps;
    }
    vector<face> solve(VP3 ps) {
        int n = ps.size();
        for (int i = 0; i < n; i++) {
            ps[i] = ps[i] + Point3(randeps(), randeps(), randeps());
        }
        VVI vs(n, VI(n));
        vector<face> crt;
        crt.push_back((face){0, 1, 2});
        crt.push_back((face){2, 1, 0});
        for (int i = 3; i < n; i++) {
            vector<face> next;
            for (int j = 0; j < (int)crt.size(); j++) {
                face & t = crt[j];
                int v = (ps[t[1]] - ps[t[0]]).det(ps[t[2]] - ps[t[0]]).dot(
                    ps[i] - ps[t[0]]) < 0 ? -1 : 1;
                if (v < 0) next.push_back(t);
                for (int k = 0; k < 3; k++) {
                    if (vs[t[(k + 1) % 3]][t[k]] == 0) {
                        vs[t[k]][t[(k + 1) % 3]] = v;
                    } else {
                        if (vs[t[(k + 1) % 3]][t[k]] != v) {
                            if (v > 0) next.push_back((face){t[k], t[(k + 1) % 3], i});
                            else next.push_back((face){t[(k + 1) % 3], t[k], i});
                        }
                        vs[t[(k + 1) % 3]][t[k]] = 0;
                    }
                }
            }
            crt = next;
        }
        return crt;
    }
}

VP3 vs;
bool cmp(CP3 p1, CP3 p2) {
    //In order to sort vs
    if (sgn(p1.x - p2.x) != 0) return p1.x < p2.x;
    if (sgn(p1.y - p2.y) != 0) return p1.y < p2.y;
    if (sgn(p1.z - p2.z) != 0) return p1.z < p2.z;
};

int main(void) {
    srand(time(NULL));
    int n;
    while (scanf("%d", &n) != EOF && n) {
        double peak = 0, area = -1;
        vs.resize(n);
        for (int i = 0; i < n; i++) {
            scanf("%lf%lf%lf", &vs[i].x, &vs[i].y, &vs[i].z);
        }
        sort(vs.begin(), vs.end(), cmp);
        n = unique(vs.begin(), vs.end()) - vs.begin();
        vs.resize(n);
        if (n <= 3) {
            printf("%.3lf %.3lf\n", 0, 0);
            continue;
        }
        vector<ConvexHull3::face> faces = ConvexHull3::solve(vs);
        int fn = faces.size();
        for (int i = 0; i < fn; i++) {
            double tmppeak = 0, tmparea = -1;
            ConvexHull3::face & f = faces[i];
            Point3 & v0 = vs[f[0]], & v1 = vs[f[1]], & v2 = vs[f[2]];
            Vector3 e1 = v1 - v0, e2 = v2 - v0;
            Vector3 vert = e1.det(e2);

            for (int l = 0; l < n; l++) {
                if (l != f[0] && l != f[1] && l != f[2]) {
                    if (abs((vs[l] - v0).dot(vert) / vert.abs()) > tmppeak) {
                        tmppeak = abs((vs[l] - v0).dot(vert) / vert.abs());
                    }
                }
            }
            if (sgn(tmppeak - peak) >= 0) {
                VP conv;
                Vector3 ex = e1.normalized(),
                    ey = e1.det(vert).normalized();
                for (int l = 0; l < n; l++) {
                    conv.push_back(Point((vs[l] - v0).dot(ex), (vs[l] - v0).dot(ey)));
                }
                double tmps = 0;
                conv = ConvexHull::solve(conv);
                int sz = conv.size();
                for (int l = 0; l < sz; l++) {
                    tmps += (conv[l].det(conv[(l + 1) % sz])) / 2;
                }
                if (tmparea < 0 || tmps < tmparea) {
                    tmparea = tmps;
                }
                if (sgn(tmppeak - peak) > 0) {
                    area = tmparea;
                } else {
                    if (area < 0 || tmparea < area) {
                        area = tmparea;
                    }
                }
                peak = tmppeak;
            }
        }
        printf("%.3lf %.3lf\n", peak + eps, area + eps);
    }
    return 0;
}
