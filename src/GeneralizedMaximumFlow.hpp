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
        Edge(unsigned int from, unsigned int to, double flow, double cap, double gain, unsigned int rev, bool is_rev) : from(from), to(to), flow(flow), cap(cap), gain(gain), cost(-log(gain)), rev(rev), is_rev(is_rev) {
            assert(this->from != this->to);
            assert(0 < this->cap);
            assert(0 < this->gain);
        }
    };

    const unsigned int num_node;          // 頂点数
    unsigned int num_edge = 0;            // 辺数
    std::vector<std::vector<Edge>> graph; // グラフの隣接リスト表現
    std::vector<double> excess;
    const double inf = std::numeric_limits<double>::max() / 3;
    const double epsilon = 1e-10;

public:
    GeneralizedMaximumFlow(unsigned int num_node) : num_node(num_node) {
        graph.resize(num_node);
        excess.resize(num_node);
    }

    // fromからtoへ向かう容量cap、コストcostの辺をグラフに追加する
    void add_edge(const unsigned int from, const unsigned int to, double cap, double gain) {
        graph.at(from).emplace_back(Edge(from, to, 0, cap, gain, graph.at(to).size(), false));
        graph.at(to).emplace_back(Edge(to, from, cap * gain, cap * gain, 1 / gain, graph.at(from).size() - 1, true));
        num_edge++;
    }

    void set_node_supply(const unsigned int node, double flow) {
        excess.at(node) += flow;
    }

    // sinkへの一般化最大流を求める
    void solve(const unsigned int sink) {
        for (Edge &edge : graph.at(sink)) {
            if (edge.from == sink) {
                assert(edge.is_rev);
            }
        }
        assert(excess.at(sink) == 0);

        for (unsigned int u = 0; u < num_node; ++u) {
            if (excess.at(u) > 0 and u != sink) {
                canceling_flow_generating_cycle(u);
            }
        }

        bool change = true;
        while (change) {
            change = false;
            for (unsigned int u = 0; u < num_node; ++u) {
                if (excess.at(u) > 0 and u != sink) {
                    change |= argument_flow(u, sink);
                }
            }
        }
    }

    // 最適解を取得
    double optimal_cost(const unsigned sink) {
        return excess.at(sink);
    }

    void show_flow() {
        for(unsigned int from = 0; from < graph.size(); ++from) {
            for (auto &edge : graph.at(from)) {
                if (not edge.is_rev) {
                    std::cout << from << "->" << edge.to << "(" << edge.flow << ")" << std::endl;
                }
            }
        }
    }

private:
    void push_flow(Edge &edge, double alpha, const std::vector<double> &labels) {
        const double flow = labels.at(edge.from) * alpha + epsilon;
        edge.flow += flow;
        edge.flow = std::max(0.0, std::min(edge.flow, edge.cap));

        Edge &rev_edge = graph.at(edge.to).at(edge.rev);
        rev_edge.flow -= flow * edge.gain + epsilon;   // 閉路の場合、実装上の都合で始点と終点が同じラベル値になってしまうので、ラベルを使わないで計算する
        rev_edge.flow = std::max(0.0, std::min(rev_edge.flow, rev_edge.cap));

        assert(0 <= edge.flow and edge.flow <= edge.cap);
        assert(0 <= rev_edge.flow and rev_edge.flow <= graph.at(edge.to).at(edge.rev).cap);
        assert(abs((edge.flow * edge.gain) - (rev_edge.cap - rev_edge.flow)) < epsilon);
    }


    bool argument_flow(const unsigned int source, const unsigned int sink) {
        std::vector<double> potential(num_node, inf);
        potential.at(source) = 0.0;

        // 辺の重みに負があるので、最初の1回はポテンシャルをベルマンフォードで求めておく
        bool update = true;
        while (update) {
            update = false;
            for (unsigned int v = 0; v < num_node; ++v) {
                if (potential.at(v) == inf) {
                    continue;
                }
                for (const auto &e : graph.at(v)) {
                    double new_cost = potential.at(v) - log(e.gain);
                    if (e.cap - e.flow > epsilon and potential.at(e.to) > new_cost + epsilon) {
                        assert(e.cap - e.flow > 0);
                        potential.at(e.to) = new_cost;
                        update = true;
                    }
                }
            }
        }

        bool change = false;
        while (excess.at(source) > epsilon) {
            // ダイクストラ法を用いて最短距離を計算
            std::vector<int> prev_v(num_node, -1), prev_e(num_node, -1); // 直前の頂点と辺のidx
            std::priority_queue<std::pair<double, unsigned int>, std::vector<std::pair<double, unsigned int>>, std::greater<std::pair<double, unsigned int>>> que;
            std::vector<double> distance(num_node, inf);      // 最短距離
            distance.at(source) = 0;

            que.push(std::make_pair(0, source));
            while (not que.empty()) {
                double dist = que.top().first;
                unsigned int node = que.top().second;
                que.pop();

                if (distance.at(node) < dist) {
                    continue;
                }
                for (int i = 0; i < graph.at(node).size(); ++i) {
                    Edge &e = graph.at(node).at(i);
                    if (e.cap - e.flow > epsilon and distance.at(e.to) > distance.at(node) + e.cost + potential.at(node) - potential.at(e.to) + epsilon) {
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
            for (unsigned int v = 0; v < num_node; ++v) {
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

            alpha = std::min(alpha, excess.at(source) / labels.at(source));

            // sourceからsinkまでflowを流す
            for (int v = sink; v != source; v = prev_v.at(v)) {
                auto &edge = graph.at(prev_v.at(v)).at(prev_e.at(v));
                push_flow(edge, alpha, labels);
            }

            excess.at(source) -= labels.at(source) * alpha;
            excess.at(sink) += alpha;
            change = true;
        }
        return change;
    }

    // 増加閉路を消去する
    void canceling_flow_generating_cycle(const unsigned int source) {

        while (true) {
            std::vector<int> prev_v(num_node, -1), prev_e(num_node, -1); // 直前の頂点と辺のidx
            std::vector<double> distance(num_node, inf);
            distance.at(source) = 0;

            int negative_cycle_idx = -1;
            for (int num = 0; num < num_node; ++num) {
                bool change = false;
                for (int u = 0; u < graph.size(); ++u) {
                    for (int i = 0; i < graph.at(u).size(); ++i) {
                        Edge &e = graph.at(u).at(i);

                        if (distance.at(u) == inf) {
                            continue;
                        }

                        double new_dist = distance.at(u) + e.cost + epsilon;
                        if (abs(e.cap - e.flow) > epsilon and distance.at(e.to) > new_dist) {
                            assert(e.cap - e.flow > 0);
                            distance.at(e.to) = new_dist;
                            prev_v.at(e.to) = u;
                            prev_e.at(e.to) = i;
                            change = true;

                            if (num == num_node - 1) {
                                negative_cycle_idx = e.to;
                            }
                        }
                    }
                }
                if (not change) {
                    break;
                }
            }

            // 負閉路がない
            if (negative_cycle_idx == -1) {
                return;
            }

            // 負閉路を取得
            int u = negative_cycle_idx;
            std::vector<unsigned int> cycle;
            std::vector<bool> used(num_node, false);
            while (not used.at(u)) {
                used.at(u) = true;
                cycle.emplace_back(u);
                u = prev_v.at(u);
            }
            negative_cycle_idx = u;

            cycle.clear();
            used.assign(num_node, false);
            while (not used.at(u)) {
                used.at(u) = true;
                cycle.emplace_back(u);
                u = prev_v.at(u);
            }
            assert(u == negative_cycle_idx);

            // alphaとlabelを算出
            double alpha = inf;
            std::vector<double> labels(num_node);
            double sum_of_gain = 0.0;
            labels.at(cycle.at(0)) = 1.0;
            for (unsigned int u : cycle) {
                const auto &edge = graph.at(prev_v.at(u)).at(prev_e.at(u));
                sum_of_gain += log(edge.gain);
                const double rest = edge.cap - edge.flow;
                assert(rest > 0);

                const double label = 1.0 / exp(sum_of_gain);
                labels.at(prev_v.at(u)) = label;
                alpha = std::min(alpha, rest / label);
            }
            assert(alpha != 0);

            for (int u : cycle) {
                Edge &edge = graph.at(prev_v.at(u)).at(prev_e.at(u));
                push_flow(edge, alpha, labels);
            }

            excess.at(u) -= labels.at(source) * alpha;
            excess.at(u) += alpha;
        }
    }
};
