#include <bits/stdc++.h>


class GeneralizedMaximumFlow {
    struct Edge {
        const unsigned int from;
        const unsigned int to;
        double flow;             // 流量
        const double cap;        // 容量
        const double gain;       // gain
        const double cost;       // cost
        const unsigned int rev;  // 逆辺のノードid
        const bool is_rev;       // 逆辺かどうか
        Edge(unsigned int from, unsigned int to, double flow, double cap, double gain, unsigned int rev, bool is_rev) : from(from), to(to), flow(flow), gain(gain), cost(-log(gain)), cap(cap), rev(rev), is_rev(is_rev) {
            assert(this->from != this->to);
            assert(0 < this->cap);
            assert(0 < this->gain);
        }
    };

    const unsigned int num_node;          // 頂点数
    unsigned int num_edge = 0;            // 辺数
    std::vector<std::vector<Edge>> graph; // グラフの隣接リスト表現
    const double inf = std::numeric_limits<double>::max() / 3;
    const double epsilon = 1e-10;

public:
    GeneralizedMaximumFlow(unsigned int num_node) : num_node(num_node) {
        graph.resize(num_node);
    }

    // fromからtoへ向かう容量cap、コストcostの辺をグラフに追加する
    void add_edge(const unsigned int from, const unsigned int to, double cap, double gain) {
        assert(gain <= 1);
        graph.at(from).emplace_back(Edge(from, to, 0, cap, gain, graph.at(to).size(), false));
        graph.at(to).emplace_back(Edge(to, from, cap, cap, 1 / gain, graph.at(from).size() - 1, true));
        num_edge++;
    }

    // sinkへの一般化最大流を求める
    void solve(const unsigned int source, const unsigned int sink, double flow) {
        assert(source != sink);
        for (Edge &edge : graph.at(sink)) {
            if (edge.from == sink) {
                assert(edge.is_rev);
            }
        }
        argument_flow(source, sink, flow);
    }

    // 最適解を取得
    double optimal_cost(const unsigned sink) {
        double flow = 0;
        for(int from = 0; from < graph.size(); ++from) {
            for (auto &edge : graph.at(from)) {
                if (not edge.is_rev and edge.to == sink) {
                    flow += edge.flow * edge.gain;
                }
            }
        }
        return flow;
    }

    void show_flow() {
        for(int from = 0; from < graph.size(); ++from) {
            for (auto &edge : graph.at(from)) {
                if (not edge.is_rev) {
                    std::cout << from << "->" << edge.to << "(" << edge.flow << ")" << std::endl;
                }
            }
        }
    }

private:
    void push_flow(Edge &edge, double alpha, const std::vector<double> &labels) {
        edge.flow += labels.at(edge.from) * alpha + epsilon;
        edge.flow = std::max(0.0, std::min(edge.flow, edge.cap));

        Edge &rev_edge = graph.at(edge.to).at(edge.rev);
        rev_edge.flow -= labels.at(rev_edge.from) * alpha + epsilon;
        rev_edge.flow = std::max(0.0, std::min(rev_edge.flow, rev_edge.cap));

        assert(0 <= edge.flow and edge.flow <= edge.cap);
        assert(0 <= rev_edge.flow and rev_edge.flow <= graph.at(edge.to).at(edge.rev).cap);
    }


    void argument_flow(const unsigned int source, const unsigned int sink, double flow) {
        std::vector<double> potential(num_node, 0);

        while (flow > epsilon) {
            // ダイクストラ法を用いて最短距離を計算
            std::vector<int> prev_v(num_node, -1), prev_e(num_node, -1); // 直前の頂点と辺のidx
            std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> que;
            std::vector<double> distance(num_node, inf);      // 最短距離
            distance.at(source) = 0;

            que.push(std::make_pair(0, source));
            while (not que.empty()) {
                double dist = que.top().first;
                int node = que.top().second;
                que.pop();

                if (distance.at(node) < dist) {
                    continue;
                }
                for (int i = 0; i < graph.at(node).size(); ++i) {
                    Edge &e = graph.at(node).at(i);
                    if (abs(e.cap - e.flow) > epsilon and distance.at(e.to) > distance.at(node) + e.cost + potential.at(node) - potential.at(e.to) + epsilon) {
                        assert(e.cap - e.flow > 0);
                        distance.at(e.to) = distance.at(node) + e.cost + potential.at(node) - potential.at(e.to);
                        prev_v.at(e.to) = node;
                        prev_e.at(e.to) = i;
                        que.push(std::make_pair(distance.at(e.to), e.to));
                    }
                }
            }

            // これ以上流せない
            if (distance[sink] == inf) {
                break;
            }

            // ポテンシャルを更新
            for (int v = 0; v < num_node; ++v) {
                potential.at(v) += distance.at(v);
            }

            // alphaとlabelを算出
            double alpha = inf;
            std::vector<double> labels(num_node);
            double sum_of_gain = 0.0;
            labels.at(sink) = 1.0;

            for (int v = sink; v != source; v = prev_v.at(v)) {
                const auto &edge = graph.at(prev_v.at(v)).at(prev_e.at(v));
                sum_of_gain += log(edge.gain);

                const double rest = edge.cap - edge.flow;
                assert(rest > 0);
                const double label = 1.0 / exp(sum_of_gain);
                alpha = std::min(alpha, rest / label);
                labels.at(prev_v.at(v)) = label;
            }

            alpha = std::min(alpha, flow / labels.at(source));

            // sourceからsinkまでflowを流す
            for (int v = sink; v != source; v = prev_v.at(v)) {
                auto &edge = graph.at(prev_v.at(v)).at(prev_e.at(v));
                push_flow(edge, alpha, labels);
            }

            flow -= labels.at(source) * alpha;
        }
    }
};
