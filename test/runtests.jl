# Last-Modified: 2026-05-03

using Test
using PhoXonic
using LinearAlgebra

const TEST_GROUP = isempty(ARGS) ? "all" : lowercase(ARGS[1])
@assert TEST_GROUP in ("all", "core", "ext") "TEST_GROUP must be all/core/ext, got: $TEST_GROUP"

if TEST_GROUP in ("all", "core")
    include("core_tests.jl")
end
if TEST_GROUP in ("all", "ext")
    include("ext_tests.jl")
end
