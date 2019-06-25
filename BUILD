cc_library(
    name = "martingale-cs",
    srcs = ["martingale-cs.c"],
    hdrs = ["martingale-cs.h"],
    visibility = ["//visibility:public"],
    deps = [],
)

cc_test(
    name = "martingale-cs_test",
    srcs = ["martingale-cs_test.cc"],
    deps = [
        ":martingale-cs",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "martingale-cs-stat_test",
    srcs = ["martingale-cs-stat_test.cc"],
    size = "enormous",  # each test needs ~15M data points.
    shard_count = 5,
    deps = [
        ":martingale-cs",
        "@com_google_googletest//:gtest_main",
        "@csm//:csm",
    ],
)
