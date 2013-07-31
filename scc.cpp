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

const int MAXN = 10000 + 10;
vector<int> oriG[MAXN], revG[MAXN], revDfsNode;
int visited[MAXN], sccIndex[MAXN];

void addEdge(int u, int v) {
    oriG[u].push_back(v);
    revG[v].push_back(u);
}

void oriDfs(int node) {
    visited[node] = 1;
    int n = oriG[node].size();
    for (int i = 0; i < n; i++) {
        if (!visited[oriG[node][i]]) {
            oriDfs(oriG[node][i]);
        }
    }
    revDfsNode.push_back(node);
}

void revDfs(int node, int cnt) {
    visited[node] = 1;
    sccIndex[node] = cnt;
    int n = revG[node].size();
    for (int i = 0; i < n; i++) {
        if (!visited[revG[node][i]]) {
            revDfs(revG[node][i], cnt);
        }
    }
}

int scc(int n) {
    memset(visited, 0, sizeof(visited));
    revDfsNode.clear();
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            oriDfs(i);
        }
    }

    memset(visited, 0, sizeof(visited));
    int cnt = 0;
    for (int i = (int)revDfsNode.size() - 1; i >= 0; i--) {
        if (!visited[revDfsNode[i]]) {
            revDfs(revDfsNode[i], cnt++);
        }
    }
    return cnt;
}

int main(void) {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        int m;
        scanf("%d", &m);
        for (int j = 0; j < m; j++) {
            int v;
            scanf("%d", &v);
            addEdge(i, v);
        }
    }
    int nSCC = scc(n);
    for (int i = 0; i < n; i++) {
        printf("%d : %d\n", i, sccIndex[i]);
    }
    return 0;
}
