#include <bits/stdc++.h>
#include "GeneralizedMaximumFlow.hpp"

using namespace std;


int main(void) {
    int V, E, F;
    cin >> V >> E >> F;

    GeneralizedMaximumFlow gmf(V);
    for (int i = 0; i < E; ++i) {
        int u, v;
        double capacity, gain;
        cin >> u >> v >> capacity >> gain;
        gmf.add_edge(u, v, capacity, gain);
    }
    unsigned int source = 0;
    unsigned int sink = V - 1;
    gmf.solve(source, sink, F);
    cout << "flow is " << gmf.optimal_cost(sink) << endl;
    gmf.show_flow();

    return 0;
}