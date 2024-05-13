import Base.CoreLogging

"""
    PreconditionerAlgorithm

Description of how to approximate a `FieldMatrix` or something similar like a
[`LazySchurComplement`](@ref) with a preconditioner `P` for which `P * x = b` is
easy to solve for `x`. If `P` is a diagonal matrix, then `x` can be computed as
`@. inv(P) * b`; otherwise, the `PreconditionerAlgorithm` must specify a
`FieldMatrixSolverAlgorithm` that can be used to solve `P * x = b` for `x`.

# Interface

Every subtype of `PreconditionerAlgorithm` must implement methods for the
following functions:
- [`solver_algorithm`](@ref)
- [`lazy_preconditioner`](@ref)
"""
abstract type PreconditionerAlgorithm end

"""
    solver_algorithm(P_alg)

A `FieldMatrixSolverAlgorithm` that can be used to solve `P * x = b` for `x`,
where `P` is the preconditioner generated by the `PreconditionerAlgorithm`
`P_alg`. If `P_alg` is `nothing` instead of a `PreconditionerAlgorithm`, or if
`P` is a diagonal matrix (and no solver is required to invert it), this returns
`nothing`.
"""
solver_algorithm(::Nothing) = nothing

is_diagonal(P_alg) = isnothing(solver_algorithm(P_alg))

"""
    lazy_preconditioner(P_alg, A)

Constructs a lazy `FieldMatrix` (or a concrete one when possible) that
approximates `A` according to the `PreconditionerAlgorithm` `P_alg`. If `P_alg`
is `nothing` instead of a `PreconditionerAlgorithm`, this returns `one(A)`.
"""
lazy_preconditioner(::Nothing, A::FieldMatrix) = one(A)

"""
    preconditioner_cache(P_alg, A, b)

Allocates the cache required to solve the equation `P * x = b`, where `P` is the
preconditioner generated by the `PreconditionerAlgorithm` `P_alg` for `A`.
"""
function preconditioner_cache(P_alg, A, b)
    is_diagonal(P_alg) && return (;)
    lazy_P = lazy_preconditioner(P_alg, A)
    P = is_lazy(lazy_P) ? Base.Broadcast.materialize(lazy_P) : lazy_P
    P_if_needed = is_lazy(lazy_P) ? (; P) : (;)
    x = similar_to_x(P, b)
    alg = solver_algorithm(P_alg)
    cache = field_matrix_solver_cache(alg, P, b)
    return (; P_if_needed..., b = similar(b), x, cache)
end

"""
    check_preconditioner(P_alg, P_cache, A, b)

Checks that `P` is compatible with `b` in the equation `P * x = b`, where `P` is
the preconditioner generated by the `PreconditionerAlgorithm` `P_alg` for `A`.
If `P_alg` requires a `FieldMatrixSolverAlgorithm` `alg` to solve the equation,
this also calls [`check_field_matrix_solver`](@ref) on `alg`.
"""
function check_preconditioner(P_alg, P_cache, A, b)
    isnothing(P_alg) && return nothing
    lazy_P = lazy_preconditioner(P_alg, A)
    if is_diagonal(P_alg)
        check_block_diagonal_matrix_has_no_missing_blocks(lazy_P, b)
    else
        alg = solver_algorithm(P_alg)
        check_field_matrix_solver(alg, P_cache.cache, lazy_P, b)
    end
end

"""
    lazy_or_concrete_preconditioner(P_alg, P_cache, A)

A wrapper for [`lazy_preconditioner`](@ref) that turns the lazy `FieldMatrix`
`P` into a concrete `FieldMatrix` when the `PreconditionerAlgorithm` `P_alg`
requires a `FieldMatrixSolverAlgorithm` to invert it.
"""
NVTX.@annotate function lazy_or_concrete_preconditioner(P_alg, P_cache, A)
    isnothing(P_alg) && return nothing
    lazy_P = lazy_preconditioner(P_alg, A)
    (is_diagonal(P_alg) || !is_lazy(lazy_P)) && return lazy_P
    @. P_cache.P = lazy_P
    return P_cache.P
end

"""
    apply_preconditioner(P_alg, P_cache, P, lazy_b)

Constructs a lazy `FieldMatrix` (or a concrete one when possible) that
represents the product `@. inv(P) * b`. Here, `lazy_b` is a (possibly lazy)
`FieldVectorView` that represents `b`.
"""
NVTX.@annotate function apply_preconditioner(P_alg, P_cache, P, lazy_b)
    isnothing(P_alg) && return lazy_b
    is_diagonal(P_alg) && return lazy_mul(lazy_inv(P), lazy_b)
    @. P_cache.b = lazy_b
    alg = solver_algorithm(P_alg)
    run_field_matrix_solver!(alg, P_cache.cache, P_cache.x, P, P_cache.b)
    return P_cache.x
end

################################################################################

"""
    MainDiagonalPreconditioner()

A `PreconditionerAlgorithm` that sets `P` to the main diagonal of `A`. This is
also called a "Jacobi" preconditioner.
"""
struct MainDiagonalPreconditioner <: PreconditionerAlgorithm end

solver_algorithm(::MainDiagonalPreconditioner) = nothing
lazy_preconditioner(::MainDiagonalPreconditioner, A::FieldMatrix) =
    lazy_main_diagonal(A)

"""
    BlockDiagonalPreconditioner()

A `PreconditionerAlgorithm` that sets `P` to the block diagonal entries of `A`.
This is also called a "block Jacobi" preconditioner.
"""
struct BlockDiagonalPreconditioner <: PreconditionerAlgorithm end

solver_algorithm(::BlockDiagonalPreconditioner) = BlockDiagonalSolve()
lazy_preconditioner(::BlockDiagonalPreconditioner, A::FieldMatrix) =
    A[matrix_diagonal_keys(keys(A))]

"""
    BlockArrowheadPreconditioner(names₁...; [P_alg₁], [alg₂])

A `PreconditionerAlgorithm` for a 2×2 block matrix:
```math
A = \\begin{bmatrix} A_{11} & A_{12} \\\\ A_{21} & A_{22} \\end{bmatrix}
```
The `FieldName`s in `names₁` correspond to the subscript `₁`, while all other
`FieldName`s correspond to the subscript `₂`. The preconditioner `P` is set to
the following matrix:
```math
P = \\begin{bmatrix} P_{11} & A_{12} \\\\ A_{21} & A_{22} \\end{bmatrix}, \\quad
\\text{where } P_{11} \\text{ is a diagonal matrix}
```
The internal preconditioner `P₁₁` is generated by the `PreconditionerAlgorithm`
`P_alg₁`, which is set to a [`MainDiagonalPreconditioner`](@ref) by default. The
Schur complement of `P₁₁` in `P`, `A₂₂ - A₂₁ * inv(P₁₁) * A₁₂`, is inverted
using the `FieldMatrixSolverAlgorithm` `alg₂`, which is set to a
[`BlockDiagonalSolve`](@ref) by default.
"""
struct BlockArrowheadPreconditioner{
    N <: NTuple{<:Any, FieldName},
    P <: Union{Nothing, PreconditionerAlgorithm},
    A <: FieldMatrixSolverAlgorithm,
} <: PreconditionerAlgorithm
    names₁::N
    P_alg₁::P
    alg₂::A
end
function BlockArrowheadPreconditioner(
    names₁...;
    P_alg₁ = MainDiagonalPreconditioner(),
    alg₂ = BlockDiagonalSolve(),
)
    is_diagonal(P_alg₁) ||
        error("BlockArrowheadPreconditioner requires a preconditioner P_alg₁ \
               that generates a diagonal matrix")
    return BlockArrowheadPreconditioner(names₁, P_alg₁, alg₂)
end

solver_algorithm(P_alg::BlockArrowheadPreconditioner) =
    BlockArrowheadSolve(P_alg.names₁...; P_alg.alg₂)
function lazy_preconditioner(
    P_alg::BlockArrowheadPreconditioner,
    A::FieldMatrix,
)
    A₁₁, A₁₂, A₂₁, A₂₂ = partition_blocks(P_alg.names₁, A)
    lazy_P₁₁ = lazy_preconditioner(P_alg.P_alg₁, A₁₁)
    return lazy_add(lazy_P₁₁, A₁₂, A₂₁, A₂₂)
end

"""
    BlockArrowheadSchurComplementPreconditioner(; [P_alg₁], [alg₂])

A `PreconditionerAlgorithm` that is equivalent to a
[`BlockArrowheadPreconditioner`](@ref), but only applied to the Schur complement
of `A₁₁` in `A`, `A₂₂ - A₂₁ * inv(A₁₁) * A₁₂`, which is represented by a
[`LazySchurComplement`](@ref). Specifically, the preconditioner this generates
is the Schur complement of `P₁₁` in `P`, `A₂₂ - A₂₁ * inv(P₁₁) * A₁₂`, where
`P₁₁` is generated by `P_alg₁`. Unlike the `BlockArrowheadPreconditioner`
constructor, this constructor does not require `names₁` because the block
structure of `A` can be inferred from the `LazySchurComplement`.
"""
struct BlockArrowheadSchurComplementPreconditioner{
    P <: Union{Nothing, PreconditionerAlgorithm},
    A <: FieldMatrixSolverAlgorithm,
} <: PreconditionerAlgorithm
    P_alg₁::P
    alg₂::A
end
function BlockArrowheadSchurComplementPreconditioner(;
    P_alg₁ = MainDiagonalPreconditioner(),
    alg₂ = BlockDiagonalSolve(),
)
    is_diagonal(P_alg₁) ||
        error("BlockArrowheadSchurComplementPreconditioner requires a \
               preconditioner P_alg₁ that generates a diagonal matrix")
    return BlockArrowheadSchurComplementPreconditioner(P_alg₁, alg₂)
end

solver_algorithm(P_alg₂::BlockArrowheadSchurComplementPreconditioner) =
    P_alg₂.alg₂
function lazy_preconditioner(
    P_alg₂::BlockArrowheadSchurComplementPreconditioner,
    A₂₂′::LazySchurComplement,
)
    (; A₁₁, A₁₂, A₂₁, A₂₂) = A₂₂′
    lazy_P₁₁ = lazy_preconditioner(P_alg₂.P_alg₁, A₁₁)
    return lazy_sub(A₂₂, lazy_mul(A₂₁, lazy_inv(lazy_P₁₁), A₁₂))
end

"""
    WeightedPreconditioner(M, unweighted_P_alg)

A `PreconditionerAlgorithm` that sets `P` to `M * P′`, where `M` is a diagonal
`FieldMatrix` and `P′` is the preconditioner generated by `unweighted_P_alg`.
When the entries of `M` are larger than 1, this is called "relaxation" or
"damping"; when the entries are smaller than 1, this is called "extrapolation".
"""
struct WeightedPreconditioner{M <: FieldMatrix, P <: PreconditionerAlgorithm} <:
       PreconditionerAlgorithm
    M::M
    unweighted_P_alg::P
end
function WeightedPreconditioner(M, unweighted_P_alg)
    check_diagonal_matrix(
        M,
        "WeightedPreconditioner cannot use M as a weighting matrix because it",
    )
    P = typeof(unweighted_P_alg)
    return WeightedPreconditioner{typeof(M), P}(M, unweighted_P_alg)
end

solver_algorithm(P_alg::WeightedPreconditioner) =
    solver_algorithm(P_alg.unweighted_P_alg)
lazy_preconditioner(P_alg::WeightedPreconditioner, A) =
    lazy_mul(P_alg.M, lazy_preconditioner(P_alg.unweighted_P_alg, A))

"""
    CustomPreconditioner(M, alg)

A `PreconditionerAlgorithm` that sets `P` to the `FieldMatrix` `M` and inverts
`P` using the `FieldMatrixSolverAlgorithm` `alg`.
"""
struct CustomPreconditioner{
    M <: FieldMatrix,
    A <: Union{Nothing, FieldMatrixSolverAlgorithm},
} <: PreconditionerAlgorithm
    M::M
    alg::A
end
function CustomPreconditioner(M; alg = nothing)
    isnothing(alg) && check_diagonal_matrix(
        M,
        "CustomPreconditioner requires alg to be specified for the matrix M \
         because it",
    )
    return CustomPreconditioner(M, alg)
end

solver_algorithm(P_alg::CustomPreconditioner) = P_alg.alg
lazy_preconditioner(P_alg::CustomPreconditioner, _) = P_alg.M

################################################################################

"""
    StationaryIterativeSolve(; [kwargs...])

A `LazyFieldMatrixSolverAlgorithm` that solves `A * x = b` by setting `x` to
some initial value `x[0]` (usually just the zero vector, ``\\mathbf{0}``) and
then iteratively updating it to
```math
x[n] = x[n - 1] + \\textrm{inv}(P) * (b - A * x[n - 1]).
```
The matrix `P` is called a "left preconditioner" for `A`. In general, this
algorithm converges more quickly when `P` is a close approximation of `A`,
although more complicated approximations often come with a performance penalty.

# Background

Let `x'` denote the value of `x` for which `A * x = b`. Replacing `b` with
`A * x'` in the formula for `x[n]` tells us that
```math
x[n] = x' + (I - \\textrm{inv}(P) * A) * (x[n - 1] - x').
```
In other words, the error on iteration `n`, `x[n] - x'`, can be expressed in
terms of the error on the previous iteration, `x[n - 1] - x'`, as
```math
x[n] - x' = (I - \\textrm{inv}(P) * A) * (x[n - 1] - x').
```
By induction, this means that the error on iteration `n` is
```math
x[n] - x' = (I - \\textrm{inv}(P) * A)^n * (x[0] - x').
```
If we pick some norm ``||\\cdot||``, we find that the norm of the error is
bounded by
```math
||x[n] - x'|| ≤ ||(I - \\textrm{inv}(P) * A)^n|| * ||x[0] - x'||.
```
For any matrix ``M``, the spectral radius of ``M`` is defined as
```math
\\rho(M) = \\max\\{|λ| : λ \\text{ is an eigenvalue of } M\\}.
```
The spectral radius has the property that
```math
||M^n|| \\sim \\rho(M)^n, \\quad \\text{i.e.,} \\quad
\\lim_{n \\to \\infty} \\frac{||M^n||}{\\rho(M)^n} = 1.
```
So, as the value of `n` increases, the norm of the error becomes bounded by
```math
||x[n] - x'|| \\leq \\rho(I - \\textrm{inv}(P) * A)^n * ||x[0] - x'||.
```
This indicates that `x[n]` will converge to `x'` (i.e., that the norm of the
error will converge to 0) when `ρ(I - inv(P) * A) < 1`, and that the convergence
rate is roughly bounded by `ρ(I - inv(P) * A)` for large values of `n`. More
precisely, it can be shown that `x[n]` will converge to `x'` if and only if
`ρ(I - inv(P) * A) < 1`. In practice, though, the convergence eventually stops
due to the limits of floating point precision.

Also, if we assume that `x[n] ≈ x'`, we can use the formula for `x[n]` to
approximate the error on the previous iteration as
```math
x[n - 1] - x' ≈ x[n - 1] - x[n] = \\textrm{inv}(P) * (A * x[n - 1] - b).
```

# Debugging

This algorithm has support for 2 debugging message group names, which can be
passed to the environment variable `JULIA_DEBUG`:
- `error_norm`: prints `||x[n] - x'||₂` on every iteration, approximating the
  error `x[n] - x'` as described above
- `spectral_radius`: prints `ρ(I - inv(P) * A)`, approximating this value with
  the `eigsolve` function from KrylovKit.jl

Because the `eigsolve` function is not compatible with CUDA, debugging the
spectral radius is currently not possible on GPUs.

# Keyword Arguments

There are 4 values that can be included in `kwargs...`:
- `P_alg = nothing`: a `PreconditionerAlgorithm` that specifies how to compute
  `P` and solve `P * x = b` for `x`, or `nothing` if preconditioning
  is not required (in which case `P` is effectively set to `one(A)`)
- `n_iters = 1`: the number of iterations
- `correlated_solves = false`: whether to set `x[0]` to a value of `x` that was
  generated during an earlier call to `field_matrix_solve!`, instead of setting
  it to ``\\mathbf{0}`` (it is always set to ``\\mathbf{0}`` on the first call
  to `field_matrix_solve!`)
- `eigsolve_kwargs = (;)`: keyword arguments for the `eigsolve` function that
  can be used to tune its accuracy and speed (only applicable when debugging the
  spectral radius)
- `debug = nothing`: keyword argument for printing debug information to `@debug`.
  By default, `debug` is set to true if `"error_norm"` or `"spectral_radius"` is in
  `ENV["JULIA_DEBUG"]`, which must be set by users.
"""
struct StationaryIterativeSolve{
    correlated_solves,
    debug,
    P <: Union{Nothing, PreconditionerAlgorithm},
    K <: NamedTuple,
} <: LazyFieldMatrixSolverAlgorithm
    P_alg::P
    n_iters::Int
    eigsolve_kwargs::K
end
function StationaryIterativeSolve(;
    P_alg = nothing,
    n_iters = 1,
    debug = nothing,
    correlated_solves = false,
    eigsolve_kwargs = (;),
)
    B =
        CoreLogging.min_enabled_level(CoreLogging.current_logger()) <=
        CoreLogging.Debug
    _debug = if isnothing(debug)
        B
    else
        B || @warn "`ENV[\"JULIA_DEBUG\"]` must be set to use `debug = true`."
        debug
    end
    # Since Field operations can be much slower than typical Array operations,
    # the default values of krylovdim and maxiter specified in KrylovKit.jl
    # should be replaced with smaller values, and the default value of tol
    # should be replaced with a larger value.
    eigsolve_kwargs′ =
        (; krylovdim = 4, maxiter = 20, tol = 0.01, eigsolve_kwargs...)
    # Make correlated_solves into a type parameter to ensure type-stability.
    params =
        (correlated_solves, _debug, typeof(P_alg), typeof(eigsolve_kwargs′))
    return StationaryIterativeSolve{params...}(P_alg, n_iters, eigsolve_kwargs′)
end
get_debug(::StationaryIterativeSolve{CS, debug}) where {CS, debug} = debug

# Extract correlated_solves as if it were a regular field, not a type parameter.
Base.getproperty(
    alg::StationaryIterativeSolve{correlated_solves},
    name::Symbol,
) where {correlated_solves} =
    name == :correlated_solves ? correlated_solves : getfield(alg, name)

function field_matrix_solver_cache(alg::StationaryIterativeSolve, A, b)
    P_cache = preconditioner_cache(alg.P_alg, A, b)
    # Note: We cannot use similar_to_x here because it doesn't work for some
    # particularly complicated unit tests. For now, we will assume that x is
    # similar to b, rather than just keys(x) == keys(b).
    previous_x_cache =
        alg.correlated_solves ? (; previous_x = zero.(similar(b))) : (;)
    return (; P_cache, previous_x_cache...)
end

check_field_matrix_solver(alg::StationaryIterativeSolve, cache, A, b) =
    check_preconditioner(alg.P_alg, cache.P_cache, A, b)

is_CuArray_type(::Type{T}) where {T} = false
NVTX.@annotate function run_field_matrix_solver!(
    alg::StationaryIterativeSolve,
    cache,
    x,
    A,
    b,
)
    P = lazy_or_concrete_preconditioner(alg.P_alg, cache.P_cache, A)
    using_cuda =
        is_CuArray_type(ClimaComms.array_type(concrete_field_vector(b)))
    if get_debug(alg) && !using_cuda
        @debug begin
            e₀ = concrete_field_vector(b) # Initialize e to any nonzero vector.
            λs, _, info = KrylovKit.eigsolve(e₀, 1; alg.eigsolve_kwargs...) do e
                e_view = field_vector_view(e, keys(b).name_tree)
                lazy_Ae = lazy_mul(A, e_view)
                lazy_Δe = apply_preconditioner(
                    alg.P_alg,
                    cache.P_cache,
                    P,
                    lazy_Ae,
                )
                concrete_field_vector(@. e_view - lazy_Δe) # (I - inv(P) * A) * e
            end
            if info.converged == 0
                (; tol, maxiter) = alg.eigsolve_kwargs
                "Unable to approximate ρ(I - inv(P) * A) to within a tolerance \
                of $(100 * tol) % in $maxiter or fewer iterations"
            else
                "ρ(I - inv(P) * A) ≈ $(abs(λs[1]))"
            end
        end _group = :spectral_radius
    end
    if alg.correlated_solves
        @. x = cache.previous_x
    else
        @. x = zero(x)
    end
    for iter in 1:(alg.n_iters)
        lazy_Δb = lazy_sub(b, lazy_mul(A, x))
        lazy_Δx = apply_preconditioner(alg.P_alg, cache.P_cache, P, lazy_Δb)
        if get_debug(alg)
            @debug begin
                norm_Δx = norm(concrete_field_vector(Base.materialize(lazy_Δx)))
                "||x[$(iter - 1)] - x'||₂ ≈ $norm_Δx"
            end _group = :error_norm
        end
        @. x += lazy_Δx
    end
    if get_debug(alg)
        @debug begin
            lazy_Δb = lazy_sub(b, lazy_mul(A, x))
            lazy_Δx = apply_preconditioner(alg.P_alg, cache.P_cache, P, lazy_Δb)
            norm_Δx = norm(concrete_field_vector(Base.materialize(lazy_Δx)))
            "||x[$(alg.n_iters)] - x'||₂ ≈ $norm_Δx"
        end _group = :error_norm
    end
    if alg.correlated_solves
        @. cache.previous_x = x
    end
end

"""
    ApproximateBlockArrowheadIterativeSolve(names₁...; [P_alg₁], [alg₁], [alg₂], [kwargs...])

Shorthand for constructing a [`SchurComplementReductionSolve`](@ref) whose
`alg₂` is set to a [`StationaryIterativeSolve`](@ref) with a
[`BlockArrowheadSchurComplementPreconditioner`](@ref). The keyword argument
`alg₁` is passed to the constructor for `SchurComplementReductionSolve`, the
keyword arguments `P_alg₁` and `alg₂` are passed to the constructor for
`BlockArrowheadSchurComplementPreconditioner`, and all other keyword arguments
are passed to the constructor for `StationaryIterativeSolve`.

This algorithm is somewhat similar to a `StationaryIterativeSolve` with a
[`BlockArrowheadPreconditioner`](@ref), but it usually converges much more
quickly, i.e., the spectral radius of its iteration matrix (`I - inv(P) * A`)
tends to be significantly smaller. Roughly speaking, this is because it runs the
iterative solver on an equation with fewer variables (the Schur complement
equation, `(A₂₂ - A₂₁ * inv(A₁₁) * A₁₂) * x₂ = b₂′`), which means that, on each
iteration, it accumulates less error due to coupling between variables. However,
even though it converges more quickly, its iterations take longer because they
involve using `alg₁` to invert `A₁₁`. So, when only a few iterations are needed,
a `StationaryIterativeSolve` with a `BlockArrowheadPreconditioner` might be more
performant.

This algorithm is an example of a "segregated" solve, in contrast to the
alternative "coupled" solve. In the context of computational fluid dynamics,
this algorithm can also be viewed as a "SIMPLE" (Semi-Implicit Method for
Pressure-Linked Equations) scheme.  For more information, see Sections 4, 5, and
10 of [Numerical solution of saddle point problems](@cite Benzi2005).
"""
function ApproximateBlockArrowheadIterativeSolve(
    names₁...;
    P_alg₁ = MainDiagonalPreconditioner(),
    alg₁ = BlockDiagonalSolve(),
    alg₂ = BlockDiagonalSolve(),
    kwargs...,
)
    P_alg₂ = BlockArrowheadSchurComplementPreconditioner(; P_alg₁, alg₂)
    outer_alg₂ = StationaryIterativeSolve(; P_alg = P_alg₂, kwargs...)
    return SchurComplementReductionSolve(names₁...; alg₁, alg₂ = outer_alg₂)
end
