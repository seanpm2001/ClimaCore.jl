#=
julia --project
using Revise; include(joinpath("test", "Grids", "convenience_constructors.jl"))
=#
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore: Grids, Topologies, Meshes, DataLayouts
using Test

@testset "Convenience constructors" begin
    grid = Grids.ExtrudedCubedSphereGrid(;
        z_elem = 10,
        z_min = 0,
        z_max = 1,
        radius = 10,
        h_elem = 10,
        n_quad_points = 4,
        horizontal_layout_type = DataLayouts.IJHF,
    )
    @test grid isa Grids.ExtrudedFiniteDifferenceGrid
    @test grid.horizontal_grid isa Grids.SpectralElementGrid2D
    @test Grids.topology(grid.horizontal_grid).mesh isa
          Meshes.EquiangularCubedSphere

    grid = Grids.CubedSphereGrid(; radius = 10, n_quad_points = 4, h_elem = 10)
    @test grid isa Grids.SpectralElementGrid2D
    @test Grids.topology(grid).mesh isa Meshes.EquiangularCubedSphere

    grid = Grids.ColumnGrid(; z_elem = 10, z_min = 0, z_max = 1)
    @test grid isa Grids.FiniteDifferenceGrid

    grid = Grids.Box3DGrid(;
        z_elem = 10,
        x_min = 0,
        x_max = 1,
        y_min = 0,
        y_max = 1,
        z_min = 0,
        z_max = 10,
        periodic_x = false,
        periodic_y = false,
        n_quad_points = 4,
        x_elem = 3,
        y_elem = 4,
    )
    @test grid isa Grids.ExtrudedFiniteDifferenceGrid
    @test grid.horizontal_grid isa Grids.SpectralElementGrid2D
    @test Grids.topology(grid.horizontal_grid).mesh isa Meshes.RectilinearMesh

    grid = Grids.SliceXZGrid(;
        z_elem = 10,
        x_min = 0,
        x_max = 1,
        z_min = 0,
        z_max = 1,
        periodic_x = false,
        n_quad_points = 4,
        x_elem = 4,
    )
    @test grid isa Grids.ExtrudedFiniteDifferenceGrid
    @test grid.horizontal_grid isa Grids.SpectralElementGrid1D
    @test Grids.topology(grid.horizontal_grid).mesh isa Meshes.IntervalMesh

    grid = Grids.RectangleXYGrid(;
        x_min = 0,
        x_max = 1,
        y_min = 0,
        y_max = 1,
        periodic_x = false,
        periodic_y = false,
        n_quad_points = 4,
        x_elem = 3,
        y_elem = 4,
    )
    @test grid isa Grids.SpectralElementGrid2D
    @test Grids.topology(grid).mesh isa Meshes.RectilinearMesh
end
