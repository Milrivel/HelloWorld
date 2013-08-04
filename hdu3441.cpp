#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
using namespace std;

#define rep(i, n) for (int i = 0; i < (int)n; i++)
#define clr(a, x) memset(a, x, sizeof(a))

const int MOD = 1000000007;
long long powmod(long long a, long long n) {
    a %= MOD;
    long long ret = 1;
    while (n) {
        if (n & 1) {
            ret = ret * a % MOD;
        }
        a = a * a % MOD;
        n >>= 1;
    }
    return ret;
}
long long inv(long long n) {
    return powmod(n, MOD - 2);
}
void get_fac(int n, int fac[], int cnt[], int & tot) {
    int tmp = (int)((double)sqrt(n) + 1);
    tot = 0;
    clr(fac, 0);
    clr(cnt, 0);
    for (int i = 2; i <= tmp; i++) {
        if (n % i == 0) {
            fac[tot] = i;
            cnt[tot] = 0;
            while (n % i == 0) {
                cnt[tot]++;
                n /= i;
            }
            tot++;
        }
    }
    if (n != 1) {
        fac[tot] = n;
        cnt[tot] = 1;
        tot++;
    }
}

int fac1[55], cnt1[55], tot1;
int fac2[55], cnt2[55], tot2;
int fac3[55], cnt3[55], tot3;
int kfac[55], kcnt[55], ktot;
long long squarecalc(long long n, long long c) {
    long long ret = 0;
    ret += powmod(c, n * n);
    ret += 2 * powmod(c, (n * n + 3) / 4) + powmod(c, (n * n + 1) / 2);
    ret = ret % MOD * inv(4) % MOD;
    return ret;
}
long long euler(long long n) {
    rep(it, ktot) {
        if (n % kfac[it] == 0) {
            n -= n / kfac[it];
        }
    }
    return n % MOD;
}
void ddfs(int cur, long long & ret, long long tmp, long long k, long long m) {
    if (cur == ktot) {
        ret += powmod(m, k / tmp) % MOD * euler(tmp) % MOD;
        ret %= MOD;
        return;
    }
    for (int i = 0; i <= kcnt[cur]; i++) {
        ddfs(cur + 1, ret, tmp, k, m);
        tmp *= kfac[cur];
    }
}
long long calc(long long n, int a, int c) {
    long long k = ((long long)a * a - 1) / (n * n);
    long long m = squarecalc(n, c);
    long long tmpk = k;
    ktot = 0;
    clr(kfac, 0);
    clr(kcnt, 0);
    rep(it, tot3) {
        if (tmpk % fac3[it] == 0) {
            kfac[ktot] = fac3[it];
            while (tmpk % fac3[it] == 0) {
                kcnt[ktot]++;
                tmpk /= fac3[it];
            }
            ktot++;
        }
    }
    long long ret = 0;
    ddfs(0, ret, 1, k, m);
    ret = ret * inv(k) % MOD * c % MOD;
    return ret;
}
void dfs(int cur, long long & ans, long long tmp, int a, int c) {
    if (cur == tot3) {
        ans = (ans + calc(tmp, a, c)) % MOD;
        return;
    }
    for (int i = 0; i <= cnt3[cur]; i += 2) {
        dfs(cur + 1, ans, tmp, a, c);
        tmp *= fac3[cur];
    }
}
int solve(int a, int c) {
    if (a == 1) return c;
    get_fac(a - 1, fac1, cnt1, tot1);
    get_fac(a + 1, fac2, cnt2, tot2);
    clr(fac3, 0);
    clr(cnt3, 0);
    tot3 = 0;
    rep(it, tot1) {
        fac3[tot3++] = fac1[it];
    }
    rep(it, tot2) {
        fac3[tot3++] = fac2[it];
    }
    sort(fac3, fac3 + tot3);
    tot3 = unique(fac3, fac3 + tot3) - fac3;
    rep(it, tot3) {
        rep(it1, tot1) {
            if (fac3[it] == fac1[it1]) {
                cnt3[it] += cnt1[it1];
            }
        }
        rep(it2, tot2) {
            if (fac3[it] == fac2[it2]) {
                cnt3[it] += cnt2[it2];
            }
        }
    }
    long long ans = 0;
    dfs(0, ans, 1, a, c);
    return (int)((ans % MOD + MOD) % MOD);
}

int main(void) {
    int t;
    scanf("%d", &t);
    rep(cases, t) {
        int a, c;
        scanf("%d%d", &a, &c);
        printf("Case %d: %d\n", cases + 1, solve(a, c));
    }
    return 0;
}
