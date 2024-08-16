#=
julia --project
using Revise; include(joinpath("test", "DataLayouts", "benchmark_copyto.jl"))
=#
using Test
using ClimaCore.DataLayouts
using BenchmarkTools
import ClimaComms
import ClimaCore
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
if ClimaComms.device() isa ClimaComms.CUDADevice
    import CUDA
    device_name = CUDA.name(CUDA.device()) # Move to ClimaComms
else
    device_name = "CPU"
end

include(joinpath(pkgdir(ClimaCore), "benchmarks/scripts/benchmark_utils.jl"))

function benchmarkcopyto!(bm, device, data, val)
    caller = string(nameof(typeof(data)))
    @info "Benchmarking $caller..."
    data_rhs = similar(data)
    fill!(data_rhs, val)
    bc = Base.Broadcast.broadcasted(identity, data_rhs)
    bcp = Base.Broadcast.broadcasted(identity, parent(data_rhs))
    trial = @benchmark ClimaComms.@cuda_sync $device Base.copyto!($data, $bc)
    t_min = minimum(trial.times) * 1e-9 # to seconds
    nreps = length(trial.times)
    n_reads_writes = DataLayouts.ncomponents(data) * 2
    push_info(
        bm;
        kernel_time_s = t_min,
        nreps = nreps,
        caller,
        problem_size = size(data),
        n_reads_writes,
    )
end

@testset "copyto! with Nf = 1" begin
    device = ClimaComms.device()
    device_zeros(args...) = ClimaComms.array_type(device)(zeros(args...))
    FT = Float64
    S = FT
    Nf = 1
    Nv = 63
    Nij = 4
    Nh = 30 * 30 * 6
    Nk = 6
    bm = Benchmark(; float_type = FT, device_name)
    data = DataF{S}(device_zeros(FT, Nf))
    benchmarkcopyto!(bm, device, data, 3)
    @test all(parent(data) .== 3)
    data = IJFH{S, Nij}(device_zeros(FT, Nij, Nij, Nf, Nh))
    benchmarkcopyto!(bm, device, data, 3)
    @test all(parent(data) .== 3)
    data = IFH{S, Nij}(device_zeros(FT, Nij, Nf, Nh))
    benchmarkcopyto!(bm, device, data, 3)
    @test all(parent(data) .== 3)
    # The parent array of IJF and IF datalayouts are MArrays, and can therefore not bm, be passed into CUDA kernels on the RHS.
    # data = IJF{S, Nij}(device_zeros(FT,Nij,Nij,Nf));             benchmarkcopyto!(bm, device, data, 3); @test all(parent(data) .== 3)
    # data = IF{S, Nij}(device_zeros(FT,Nij,Nf));                  benchmarkcopyto!(bm, device, data, 3); @test all(parent(data) .== 3)
    data = VF{S, Nv}(device_zeros(FT, Nv, Nf))
    benchmarkcopyto!(bm, device, data, 3)
    @test all(parent(data) .== 3)
    data = VIJFH{S, Nv, Nij}(device_zeros(FT, Nv, Nij, Nij, Nf, Nh))
    benchmarkcopyto!(bm, device, data, 3)
    @test all(parent(data) .== 3)
    data = VIFH{S, Nv, Nij}(device_zeros(FT, Nv, Nij, Nf, Nh))
    benchmarkcopyto!(bm, device, data, 3)
    @test all(parent(data) .== 3)

    # data = IJKFVH{S}(device_zeros(FT,Nij,Nij,Nk,Nf,Nh)); benchmarkcopyto!(bm, device, data, 3); @test all(parent(data) .== 3) # TODO: test
    # data = IH1JH2{S}(device_zeros(FT,Nij,Nij,Nk,Nf,Nh)); benchmarkcopyto!(bm, device, data, 3); @test all(parent(data) .== 3) # TODO: test
    tabulate_benchmark(bm)
end
