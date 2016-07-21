#include "fast_ig_tools.hpp"

size_t numEdges(const Graph &graph,
                bool undirected) {
    size_t nE = 0;

    for (const auto &edges : graph) {
        nE += edges.size();
    }

    if (undirected) {
        nE /= 2;
    }

    return nE;
}


void write_metis_graph(const Graph &graph,
                       const std::string &filename,
                       bool undirected) {
    std::ofstream out(filename);

    // Count the numder of vertices and the number of edges
    size_t nV = graph.size();
    size_t nE = numEdges(graph, undirected);

    out << nV << " " << nE << " 001\n"; // See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

    for (const auto &edges : graph) {
        for (const auto &edge : edges) {
            out << edge.first + 1 << " " << edge.second << " ";
        }
        out << "\n";
    }
}


void write_metis_graph(const Graph &graph,
                       const std::vector<size_t> &weights,
                       const std::string &filename,
                       bool undirected) {
    assert(graph.size() == weights.size());

    std::ofstream out(filename);

    // Count the numder of vertices and the number of edges
    size_t nV = graph.size();
    size_t nE = numEdges(graph, undirected);

    out << nV << " " << nE << " 011\n"; // See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

    for (size_t i = 0; i < graph.size(); ++i) {
        out << weights[i] << " ";
        for (const auto &edge : graph[i]) {
            out << edge.first + 1 << " " << edge.second << " ";
        }
        out << "\n";
    }
}

bool check_repr_kmers_consistancy(const std::vector<size_t> &answer,
                                  const std::vector<int> &multiplicities,
                                  size_t K, size_t n) {
    if (!std::is_sorted(answer.cbegin(), answer.cend())) return false;

    // Check answer size
    if (answer.size() != n) return false;

    // Check k-mer overlapping
    for (size_t i = 1; i < answer.size(); ++i) {
        if (answer[i] - answer[i - 1] < K) return false;
    }

    // K-mers should belong interval
    for (size_t kmer_i : answer) {
        if (kmer_i >= multiplicities.size()) return false;
    }

    return true;
}

// TODO cover by tests
// and then, refactor it!!!!!!!!111111111111oneoneone
// TODO Rename it to be consistent with the paper
// TODO int -> unsigned or size_t
std::vector<size_t> optimal_coverage(const std::vector<int> &multiplicities,
                                     size_t K, size_t n) {
    assert(n >= 1);
    assert(multiplicities.size() + K - 1 >= n * K);

    const int INF = 1 << 30; // TODO Use exact value

    std::vector<std::vector<int>> imults(n, std::vector<int>(multiplicities.size()));

    // Fill by cummin
    imults[0][0] = multiplicities[0];
    for (size_t i = 1; i < imults[0].size(); ++i) {
        imults[0][i] = std::min(multiplicities[i], imults[0][i - 1]);
    }

    for (size_t j = 1; j < n; ++j) { // n == 1 is useless
        // Kill first K*j elements
        for (size_t i = 0; i < K*j; ++i) {
            imults[j][i] = INF;
        }

        for (size_t i = K*j; i < imults[j].size(); ++i) {
            imults[j][i] = std::min(imults[j][i - 1],
                                    multiplicities[i] + imults[j - 1][i - K]);
        }
    }

    auto ans = imults[n - 1][imults[n - 1].size() - 1];

    VERIFY(ans < INF);

    std::vector<size_t> result(n);
    // Backward reconstruction
    size_t i = imults[n - 1].size() - 1;
    size_t j = n - 1;

    while (j > 0) {
        if (imults[j][i] == multiplicities[i] + imults[j - 1][i - K]) { // Take i-th element
            result[j] = i;
            i -= K;
            j -= 1;
        } else {
            i -= 1;
        }
    }
    assert(j == 0);
    // Find first element
    size_t ii = i;
    while (imults[0][i] != multiplicities[ii]) {
        --ii;
    }
    result[0] = ii;


    // Checking
    int sum = 0;
    for (size_t i : result) {
        sum += multiplicities.at(i);
    }
    assert(ans == sum);

    assert(check_repr_kmers_consistancy(result, multiplicities, K, n));

    return result;
}

// vim: ts=4:sw=4
