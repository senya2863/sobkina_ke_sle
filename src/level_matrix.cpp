#include "level_matrix.h"
#include <queue>
#include <vector>

std::vector<int> computeLevelMatrix(const SparseMatrix& A) {
    int n = A.n;
    std::vector<int> level(n, -1);

    std::queue<int> q;
    q.push(0);
    level[0] = 0;

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        for (const auto& [col, val] : A.rows[current]) {
            if (col < n && level[col] == -1) {
                level[col] = level[current] + 1;
                q.push(col);
            }
        }
    }

    return level;
}