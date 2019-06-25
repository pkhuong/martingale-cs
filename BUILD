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
