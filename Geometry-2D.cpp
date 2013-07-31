#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <algorithm>
using namespace std;

#define T double
#define CP const Point &
#define VP vector<Point>
#define VI vector<int>
#define VVI vector<VI>

const double eps = 1e-8, pi = atan2(0, -1);
int sgn(double a){return a < -eps ? -1 : a < eps ? 0 : 1;}

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

double disLP(CP p1, CP p2, CP q) {
    //点到直线的距离
    return fabs((p2 - p1).det(q - p1)) / (p2 - p1).abs();
}

double disSP(CP p1, CP p2, CP q) {
    //点到线段的距离
    if ((p2 - p1).dot(q - p1) < eps)
        return (q - p1).abs();
    if ((p1 - p2).dot(q - p2) < eps)
        return (q - p2).abs();
    return disLP(p1, p2, q);
}

bool crsSS(CP p1, CP p2, CP q1, CP q2) {
    //二线段交叉判定
    //快速排斥试验
    if (max(p1.x, p2.x) + eps < min(q1.x, q2.x)) return false;
    if (max(q1.x, q2.x) + eps < min(p1.x, p2.x)) return false;
    if (max(p1.y, p2.y) + eps < min(q1.y, q2.y)) return false;
    if (max(q1.y, q2.y) + eps < min(p1.y, p2.y)) return false;
    //跨立试验
    return sgn((p2 - p1).det(q1 - p1)) * sgn((p2 - p1).det(q2 - p1)) < eps
        && sgn((q2 - q1).det(p1 - q1)) * sgn((q2 - q1).det(p2 - q1)) < eps;
}

VP isLL(CP p1, CP p2, CP q1, CP q2){
    //二直线交点，无则返回空vector
    //返回空时有共线与相离的区别，用(p2 - p1).det(q1 - p1) == 0判断
    VP ret;
    T d = (q2 - q1).det(p2 - p1);
    if (fabs(d) < eps) return ret;
    ret.push_back(p1 + (p2 - p1) * ((q2 - q1).det(q1 - p1) / d));
    return ret;
}

VP isCL(CP c, double r, CP p1, CP p2) {
    //返回值按到p1的距离从小到大排列
    double x = (p1 - c).dot(p2 - p1);
    double y = (p2 - p1).abs2();
    double d = x * x - y * ((p1 - c).abs2() - r * r);
    if (d < -eps) return VP(0);
    if (d < 0) d = 0;
    Point q1 = p1 - (p2 - p1) * (x / y);
    Point q2 = (p2 - p1) * (sqrt(d) / y);
    VP ret;
    ret.push_back(q1 - q2);
    ret.push_back(q1 + q2);
    return ret;
}

double rad(Point p1, Point p2) {
    //有向极角差
    return atan2(p1.det(p2), p1.dot(p2));
}

int contains(VP ps, CP q) {
    int n = ps.size();
    int res = -1;
    for (int i = 0; i < n; i++) {
        Point a = ps[i] - q, b = ps[(i + 1) % n] - q;
        if (a.y > b.y) {
            Point t = a; a = b; b = t;
        }
        if (a.y < eps && b.y > eps && a.det(b) > eps) {
            res = -res;
        }
        if (abs(a.det(b)) < eps && a.dot(b) < eps) {
            return 0;
        }
    }
    return res;
}

double areaCT(Point p1, Point p2, double r) {
    VP qs = isCL(0, r, p1, p2);
    if (qs.size() == 0) return r * r * rad(p1, p2) / 2;
    bool b1 = p1.abs() > r + eps, b2 = p2.abs() > r + eps;
    if (b1 && b2) {
        if ((p1 - qs[0]).dot(p2 - qs[0]) < eps &&
            (p1 - qs[1]).dot(p2 - qs[1]) < eps) {
            return (r * r * (rad(p1, p2) - rad(qs[0], qs[1])) +
                qs[0].det(qs[1])) / 2;
        } else {
            return r * r * rad(p1, p2) / 2;
        }
    } else if (b1) {
        return (r * r * rad(p1, qs[0]) + qs[0].det(p2)) / 2;
    } else if (b2) {
        return (r * r * rad(qs[1], p2) + p1.det(qs[1])) / 2;
    } else {
        return p1.det(p2) / 2;
    }
}

double areaCT(Point p1, Point p2, Point center, double r) {
    return areaCT(p1 - center, p2 - center, r);
}

double areaCP(VP ps, Point center, double r) {
    int n = ps.size();
    double ans = 0;
    for (int i = 0; i < n; i++) {
        ans += areaCT(ps[i], ps[(i + 1) % n], center, r);
    }
    return ans;
}

VP convexCut(VP ps, CP p1, CP p2) {
    //返回凸多边形被有向直线p1p2切割后左半部分， 可改为在线半平面交
    int n = ps.size();
    VP ret(0);
    for (int i = 0; i < n; i++) {
        int d1 = sgn((p2 - p1).det(ps[i] - p1));
        int d2 = sgn((p2 - p1).det(ps[(i + 1) % n] - p1));
        if (d1 >= 0) ret.push_back(ps[i]);
        if (d1 * d2 < 0) ret.push_back(isLL(p1, p2, ps[i], ps[(i + 1) % n])[0]);
    }
    return ret;
}

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

struct Line {
    //用来做半平面交。。。。。。
    Point s, e;//起点和终点，要用构造函数初始化，不要直接读
    double k; //极角

    Line(Point ss = Point(0, 0), Point ee = Point(0, 0)) : s(ss), e(ee) {
        k = atan2(ee.y - ss.y, ee.x - ss.x);
    }
};
namespace HalfPlaneIntersection {
    // 排序增量法求半平面交，保证结果为有限凸多边形
    const int MAXL = 20010;
    Line que[MAXL];

    VP isLL(Line l1, Line l2) {
        return isLL(l1.s, l1.e, l2.s, l2.e);
    }

    struct cmp {
        bool operator () (Line a, Line b) const {
            if (fabs(a.k - b.k) > eps) return a.k < b.k;
            return ((a.s - b.s).det(b.e - b.s)) < 0;
        }
    };

    bool kequal(Line a, Line b) {
        return fabs(a.k - b.k) < eps;
    }

    VP solve(vector<Line> line) {
        sort(line.begin(), line.end(), cmp());
        int n = unique(line.begin(), line.end(), kequal) - line.begin();
        assert(n > 2);

        int head = 0, tail = 1;
        que[0] = line[0]; que[1] = line[1];
        VP ret;
        for (int i = 2; i < n; i++) {
            if (fabs((que[tail].e - que[tail].s).det(que[tail - 1].e - que[tail - 1].s)) < eps ||
                fabs((que[head].e - que[head].s).det(que[head + 1].e - que[head + 1].s)) < eps) {
                return ret;
            }
            while (head < tail && ((isLL(que[tail], que[tail - 1])[0] - line[i].s)
                .det(line[i].e - line[i].s)) > eps) tail--;
            while (head < tail && ((isLL(que[head], que[head + 1])[0] - line[i].s)
                .det(line[i].e - line[i].s)) > eps) head++;
            que[++tail] = line[i];
        }
        while (head < tail && ((isLL(que[tail], que[tail - 1])[0] - que[head].s)
            .det(que[head].e - que[head].s)) > eps) tail--;
        while (head < tail && ((isLL(que[head], que[head + 1])[0] - que[tail].s)
            .det(que[tail].e - que[tail].s)) > eps) head++;

        if (tail <= head + 1)
            return ret;
        for (int i = head; i < tail; i++) {
            ret.push_back(isLL(que[i], que[i + 1])[0]);
        }
        if (head < tail + 1)
            ret.push_back(isLL(que[head], que[tail])[0]);
        return ret;
    }
}

double convexDiameter(VP ps) {
    int n = ps.size();
    int is = 0, js = 0;
    for (int i = 1; i < n; i++) {
        if (ps[i].x > ps[is].x) is = i;
        if (ps[i].x < ps[js].x) js = i;
    }
    double maxd = (ps[is] - ps[js]).abs();
    int i = is, j = js;
    do {
        if ((ps[(i + 1) % n] - ps[i]).det(ps[(j + 1) % n] - ps[j]) >= 0) {
            j = (j + 1) % n;
        } else {
            i = (i + 1) % n;
        }
        maxd = max(maxd, (ps[i] - ps[j]).abs());
    } while (i != is || j != js);
    return maxd;
}

int main(void) {
    return 0;
}
