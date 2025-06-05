#include <CGAL_MWT/Mutable_index_priority_queue.h>
#include <algorithm>
#include <doctest/doctest.h>
#include <random>

using namespace mipq;

TEST_CASE("[MutableIndexPriorityQueue] test heapsort random keys") {
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 5.0);
    std::vector<std::size_t> sizes{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 100, 1000, 1024};
    for(std::size_t size : sizes) {
        std::vector<double> keys;
        for(std::size_t i = 0; i < size; ++i) {
            keys.push_back(dist(rng));
        }
        std::vector<std::size_t> expected_result;
        expected_result.resize(size);
        std::iota(expected_result.begin(), expected_result.end(), 0);
        std::sort(expected_result.begin(), expected_result.end(),
                  [&](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
        MutableIndexPriorityQueue mipq{[&](std::size_t s) { return keys[s]; }, size};
        std::vector<std::size_t> result;
        while(!mipq.empty()) {
            result.push_back(mipq.pop());
            mipq.verify_invariants();
        }
        CHECK(result == expected_result);
    }
}

TEST_CASE("[MutableIndexPriorityQueue] test random-order key update") {
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 5.0);
    std::vector<std::size_t> sizes{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 100, 1000, 1024};
    for(std::size_t size : sizes) {
        std::vector<double> keys1;
        std::vector<double> keys2;
        for(std::size_t i = 0; i < size; ++i) {
            keys1.push_back(dist(rng));
            keys2.push_back(dist(rng));
        }
        std::vector<std::size_t> expected_result1;
        expected_result1.resize(size);
        std::iota(expected_result1.begin(), expected_result1.end(), 0);
        std::vector<std::size_t> expected_result2 = expected_result1;
        std::sort(expected_result2.begin(), expected_result2.end(),
                  [&](std::size_t a, std::size_t b) { return keys2[a] < keys2[b]; });
        MutableIndexPriorityQueue mipq{[&](std::size_t s) { return keys1[s]; }, size};
        std::vector<std::size_t> random_update_order = expected_result1;
        std::shuffle(random_update_order.begin(), random_update_order.end(), rng);
        for(std::size_t i : random_update_order) {
            keys1[i] = keys2[i];
            mipq.notify_key_changed(i);
            mipq.verify_invariants();
        }
        std::vector<std::size_t> result;
        while(!mipq.empty()) {
            result.push_back(mipq.pop());
            mipq.verify_invariants();
        }
        CHECK(result == expected_result2);
    }
}
