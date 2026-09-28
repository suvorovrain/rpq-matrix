//
// Created by Adrián on 2/5/23.
//

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <iomanip>
#include <chrono>
#include <array>
#include <limits>
#include <stdexcept>
#include <utils.hpp>

#include <bm_baselineGB.hpp>

#include <rpq_solver.hpp>

// #define N 958844164
// #define S 5420 // 1 to 5419
// #define V 296008192 // 1 to...

int main(int argc, char **argv)
{
    GrB_Info info;
    info = GrB_init(GrB_NONBLOCKING);
    if (info != GrB_SUCCESS)
    {
        fprintf(stderr, "Initialization failed!\n");
        GrB_finalize();
        return 1;
    }
    if (argc < 5 || argc > 7)
    {
        std::cerr << "\tUsage: " << argv[0] << " <dataset> <queries> <n_preds> <n_triples> [runs] [warmup_runs]" << std::endl;
        exit(1);
    }

    std::string dataset = argv[1];
    std::string queries = argv[2];
    std::string index = dataset + ".baseline-64";
    uint n_preds = atoi(argv[3]);
    uint n_triples = atoi(argv[4]);
    uint64_t runs = 1;
    uint64_t warmup_runs = 0;
    try {
        for (int argument = 5; argument < argc; ++argument) {
            std::string value = argv[argument];
            if (value.empty() || !std::all_of(value.begin(), value.end(), [](char digit) { return digit >= '0' && digit <= '9'; })) {
                throw std::invalid_argument("invalid run count");
            }
        }
        if (argc >= 6) runs = std::stoull(argv[5]);
        if (argc >= 7) warmup_runs = std::stoull(argv[6]);
        if (runs == 0 || runs > std::numeric_limits<uint64_t>::max() - warmup_runs) throw std::out_of_range("invalid run count");
    } catch (const std::exception &) {
        std::cerr << "runs must be positive and warmup_runs must be non-negative integers" << std::endl;
        GrB_finalize();
        return 1;
    }
    rpq::run_query<bm_baselinegb::wrapper>(dataset, index, queries, n_preds, n_triples, runs, warmup_runs);
    GrB_finalize();
}
