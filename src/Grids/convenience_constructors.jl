import ..DataLayouts
#=
These convenience constructors accept integer
keyword inputs, so they are dynamically created. You may
want to use a different constructor if you're making the
object in a performance-critical section, and if you know
the type parameters at compile time.

If no convenience constructor exists, then you may need to
create a custom grid using our low-level compose-able API.
=#
check_device_context(context, device) =
    @assert ClimaComms.device(context) == device "The given device and context device do not match."

#####
##### Mesh helpers
#####

"""
    DefaultSliceXMesh(
        ::Type{<:AbstractFloat}; # defaults to Float64
        x_min::Real,
        x_max::Real,
        periodic_x::Bool,
        x_elem::Integer,
    )

A convenience constructor, which builds an `IntervalMesh`.
"""
DefaultSliceXMesh(; kwargs...) = DefaultSliceXMesh(Float64; kwargs...)
function DefaultSliceXMesh(
    ::Type{FT};
    x_min::Real,
    x_max::Real,
    periodic_x::Bool,
    x_elem::Integer,
) where {FT}

    x1boundary = (:east, :west)
    z_boundary_names = (:bottom, :top)
    h_domain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(x_min),
        Geometry.XPoint{FT}(x_max);
        periodic = periodic_x,
        boundary_names = (:east, :west),
    )
    return Meshes.IntervalMesh(h_domain; nelems = x_elem)
end

"""
    DefaultZMesh(
        ::Type{<:AbstractFloat}; # defaults to Float64
        z_min::Real,
        z_max::Real,
        z_elem::Integer,
        stretch::Meshes.StretchingRule = Meshes.Uniform(),
    )

A convenience constructor, which builds an `IntervalMesh`.
"""
DefaultZMesh(; kwargs...) = DefaultZMesh(Float64; kwargs...)
function DefaultZMesh(
    ::Type{FT};
    z_min::Real,
    z_max::Real,
    z_elem::Integer,
    stretch::Meshes.StretchingRule = Meshes.Uniform(),
) where {FT}
    z_boundary_names = (:bottom, :top)
    z_domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(z_min),
        Geometry.ZPoint{FT}(z_max);
        boundary_names = z_boundary_names,
    )
    return Meshes.IntervalMesh(z_domain, stretch; nelems = z_elem)
end

"""
    DefaultRectangleXYMesh(
        ::Type{<:AbstractFloat}; # defaults to Float64
        x_min::Real,
        x_max::Real,
        y_min::Real,
        y_max::Real,
        periodic_x::Bool,
        periodic_y::Bool,
    )

A convenience constructor, which builds a
`RectilinearMesh` with a rectangular domain
composed of interval domains.
"""
DefaultRectangleXYMesh(; kwargs...) = DefaultRectangleXYMesh(Float64; kwargs...)
function DefaultRectangleXYMesh(
    ::Type{FT};
    x_min::Real,
    x_max::Real,
    y_min::Real,
    y_max::Real,
    x_elem::Integer,
    y_elem::Integer,
    periodic_x::Bool,
    periodic_y::Bool,
) where {FT <: AbstractFloat}
    x1boundary = (:east, :west)
    x2boundary = (:south, :north)
    domain = Domains.RectangleDomain(
        Domains.IntervalDomain(
            Geometry.XPoint{FT}(x_min),
            Geometry.XPoint{FT}(x_max);
            periodic = periodic_x,
            boundary_names = x1boundary,
        ),
        Domains.IntervalDomain(
            Geometry.YPoint{FT}(y_min),
            Geometry.YPoint{FT}(y_max);
            periodic = periodic_y,
            boundary_names = x2boundary,
        ),
    )
    return Meshes.RectilinearMesh(domain, x_elem, y_elem)
end

#####
##### Grids
#####

"""
    ExtrudedCubedSphereGrid(
        ::Type{<:AbstractFloat}; # defaults to Float64
        z_elem::Integer,
        z_min::Real,
        z_max::Real,
        radius::Real,
        h_elem::Integer,
        n_quad_points::Integer,
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        stretch::Meshes.StretchingRule = Meshes.Uniform(),
        hypsography::HypsographyAdaption = Flat(),
        global_geometry::Geometry.AbstractGlobalGeometry = Geometry.ShallowSphericalGlobalGeometry(radius),
        quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
        h_mesh = Meshes.EquiangularCubedSphere(Domains.SphereDomain{FT}(radius), h_elem),
        h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(context, h_mesh),
        horizontal_layout_type = DataLayouts.IJFH,
        z_mesh::Meshes.IntervalMesh = DefaultZMesh(FT; z_min, z_max, z_elem, stretch),
        enable_bubble::Bool = false
    )

A convenience constructor, which builds an
`ExtrudedFiniteDifferenceGrid`, given:

 - `FT` the float type (defaults to `Float64`) [`Float32`, `Float64`]
 - `z_elem` the number of z-points
 - `z_min` the domain minimum along the z-direction.
 - `z_max` the domain maximum along the z-direction.
 - `radius` the radius of the cubed sphere
 - `h_elem` the number of elements per side of every panel (6 panels in total)
 - `n_quad_points` the number of quadrature points
 - `device` the [`ClimaComms.device`](@ref)
 - `context` the [`ClimaComms.context`](@ref)
 - `stretch` the mesh [`Meshes.StretchingRule`](@ref) (defaults to `Uniform`) []
 - `hypsography` the hypsography
 - `global_geometry` the global geometry (defaults to `Geometry.CartesianGlobalGeometry`)
 - `quad` the quadrature style (defaults to `Quadratures.GLL{n_quad_points}`)
 - `h_mesh` the horizontal mesh (defaults to `Meshes.EquiangularCubedSphere`)
 - `h_topology` the horizontal topology (defaults to `Topologies.Topology2D`)
 - `horizontal_layout_type` the horizontal DataLayout type (defaults to `DataLayouts.IJFH`)
 - `z_mesh` the z-mesh, defaults to an `Meshes.IntervalMesh` along `z` with given `stretch`
 - `enable_bubble` enables the "bubble correction" for more accurate element areas.

# Example usage

```julia
using ClimaCore: Grids
grid = Grids.ExtrudedCubedSphereGrid(;
    z_elem = 10,
    z_min = 0,
    z_max = 1,
    radius = 10,
    h_elem = 10,
    n_quad_points = 4,
)
```
"""
ExtrudedCubedSphereGrid(; kwargs...) =
    ExtrudedCubedSphereGrid(Float64; kwargs...)

function ExtrudedCubedSphereGrid(
    ::Type{FT};
    z_elem::Integer,
    z_min::Real,
    z_max::Real,
    radius::Real,
    h_elem::Integer,
    n_quad_points::Integer,
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    stretch::Meshes.StretchingRule = Meshes.Uniform(),
    hypsography::HypsographyAdaption = Flat(),
    global_geometry::Geometry.AbstractGlobalGeometry = Geometry.ShallowSphericalGlobalGeometry(
        radius,
    ),
    quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    h_mesh = Meshes.EquiangularCubedSphere(
        Domains.SphereDomain{FT}(radius),
        h_elem,
    ),
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(
        context,
        h_mesh,
    ),
    horizontal_layout_type = DataLayouts.IJFH,
    z_mesh::Meshes.IntervalMesh = DefaultZMesh(
        FT;
        z_min,
        z_max,
        z_elem,
        stretch,
    ),
    enable_bubble::Bool = false,
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)

    z_boundary_names = (:bottom, :top)
    h_grid = Grids.SpectralElementGrid2D(
        h_topology,
        quad;
        horizontal_layout_type,
        enable_bubble,
    )
    z_topology = Topologies.IntervalTopology(context, z_mesh)
    vertical_grid = FiniteDifferenceGrid(z_topology)
    return ExtrudedFiniteDifferenceGrid(
        h_grid,
        vertical_grid,
        hypsography,
        global_geometry,
    )
end

"""
    CubedSphereGrid(
        ::Type{<:AbstractFloat}; # defaults to Float64
        radius::Real,
        h_elem::Integer,
        n_quad_points::Integer,
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
        h_mesh = Meshes.EquiangularCubedSphere(Domains.SphereDomain{FT}(radius), h_elem),
        h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(context, h_mesh),
        horizontal_layout_type = DataLayouts.IJFH,
    )

A convenience constructor, which builds a
`SpectralElementGrid2D` given:

 - `FT` the float type (defaults to `Float64`) [`Float32`, `Float64`]
 - `radius` the radius of the cubed sphere
 - `h_elem` the number of elements per side of every panel (6 panels in total)
 - `n_quad_points` the number of quadrature points
 - `device` the [`ClimaComms.device`](@ref)
 - `context` the [`ClimaComms.context`](@ref)
 - `quad` the quadrature style (defaults to `Quadratures.GLL{n_quad_points}`)
 - `h_mesh` the horizontal mesh (defaults to `Meshes.EquiangularCubedSphere`)
 - `h_topology` the horizontal topology (defaults to `Topologies.Topology2D`)
 - `horizontal_layout_type` the horizontal DataLayout type (defaults to `DataLayouts.IJFH`)

# Example usage

```julia
using ClimaCore: Grids
grid = Grids.CubedSphereGrid(; radius = 10, n_quad_points = 4, h_elem = 10)
```
"""
CubedSphereGrid(; kwargs...) = CubedSphereGrid(Float64; kwargs...)
function CubedSphereGrid(
    ::Type{FT};
    radius::Real,
    h_elem::Integer,
    n_quad_points::Integer,
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    h_mesh = Meshes.EquiangularCubedSphere(
        Domains.SphereDomain{FT}(radius),
        h_elem,
    ),
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(
        context,
        h_mesh,
    ),
    horizontal_layout_type = DataLayouts.IJFH,
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)
    return Grids.SpectralElementGrid2D(h_topology, quad; horizontal_layout_type)
end

"""
    ColumnGrid(
        ::Type{<:AbstractFloat}; # defaults to Float64
        z_elem::Integer,
        z_min::Real,
        z_max::Real,
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        stretch::Meshes.StretchingRule = Meshes.Uniform(),
        z_mesh::Meshes.IntervalMesh = DefaultZMesh(FT; z_min, z_max, z_elem, stretch),
    )

A convenience constructor, which builds a
`FiniteDifferenceGrid`.
"""
ColumnGrid(; kwargs...) = ColumnGrid(Float64; kwargs...)
function ColumnGrid(
    ::Type{FT};
    z_elem::Integer,
    z_min::Real,
    z_max::Real,
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    stretch::Meshes.StretchingRule = Meshes.Uniform(),
    z_mesh::Meshes.IntervalMesh = DefaultZMesh(
        FT;
        z_min,
        z_max,
        z_elem,
        stretch,
    ),
) where {FT}
    check_device_context(context, device)
    z_topology = Topologies.IntervalTopology(context, z_mesh)
    return FiniteDifferenceGrid(z_topology)
end

"""
    Box3DGrid(
        ::Type{<:AbstractFloat}; # defaults to Float64
        z_elem::Integer,
        x_min::Real,
        x_max::Real,
        y_min::Real,
        y_max::Real,
        z_min::Real,
        z_max::Real,
        periodic_x::Bool,
        periodic_y::Bool,
        n_quad_points::Integer,
        x_elem::Integer,
        y_elem::Integer,
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        stretch::Meshes.StretchingRule = Meshes.Uniform(),
        hypsography::HypsographyAdaption = Flat(),
        global_geometry::Geometry.AbstractGlobalGeometry = Geometry.CartesianGlobalGeometry(),
        quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
        horizontal_layout_type = DataLayouts.IJFH,
        [h_topology::Topologies.AbstractDistributedTopology], # optional
        [z_mesh::Meshes.IntervalMesh], # optional
        enable_bubble::Bool = false,
    )

A convenience constructor, which builds a
`ExtrudedFiniteDifferenceGrid` with a
`FiniteDifferenceGrid` vertical grid and a
`SpectralElementGrid2D` horizontal grid, given:

 - `z_elem` the number of z-points
 - `x_min` the domain minimum along the x-direction.
 - `x_max` the domain maximum along the x-direction.
 - `y_min` the domain minimum along the y-direction.
 - `y_max` the domain maximum along the y-direction.
 - `z_min` the domain minimum along the z-direction.
 - `z_max` the domain maximum along the z-direction.
 - `periodic_x` Bool indicating to use periodic domain along x-direction
 - `periodic_y` Bool indicating to use periodic domain along y-direction
 - `n_quad_points` the number of quadrature points
 - `x_elem` the number of x-points
 - `y_elem` the number of y-points
 - `device` the [`ClimaComms.device`](@ref)
 - `context` the [`ClimaComms.context`](@ref)
 - `stretch` the mesh [`Meshes.StretchingRule`](@ref) (defaults to `Uniform`) []
 - `hypsography` the hypsography
 - `global_geometry` the global geometry (defaults to `Geometry.CartesianGlobalGeometry`)
 - `quad` the quadrature style (defaults to `Quadratures.GLL{n_quad_points}`)
 - `h_topology` the horizontal topology (defaults to `Topologies.Topology2D`)
 - `z_mesh` the z-mesh, defaults to an `Meshes.IntervalMesh` along `z` with given `stretch`
 - `enable_bubble` enables the "bubble correction" for more accurate element areas.
 - `horizontal_layout_type` the horizontal DataLayout type (defaults to `DataLayouts.IJFH`)

# Example usage

```julia
using ClimaCore: Grids
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
```
"""
Box3DGrid(; kwargs...) = Box3DGrid(Float64; kwargs...)
function Box3DGrid(
    ::Type{FT};
    z_elem::Integer,
    x_min::Real,
    x_max::Real,
    y_min::Real,
    y_max::Real,
    z_min::Real,
    z_max::Real,
    periodic_x::Bool,
    periodic_y::Bool,
    n_quad_points::Integer,
    x_elem::Integer,
    y_elem::Integer,
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    stretch::Meshes.StretchingRule = Meshes.Uniform(),
    hypsography::HypsographyAdaption = Flat(),
    global_geometry::Geometry.AbstractGlobalGeometry = Geometry.CartesianGlobalGeometry(),
    quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(
        context,
        DefaultRectangleXYMesh(
            FT;
            x_min,
            x_max,
            y_min,
            y_max,
            x_elem,
            y_elem,
            periodic_x,
            periodic_y,
        ),
    ),
    z_mesh::Meshes.IntervalMesh = DefaultZMesh(
        FT;
        z_min,
        z_max,
        z_elem,
        stretch,
    ),
    enable_bubble::Bool = false,
    horizontal_layout_type = DataLayouts.IJFH,
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)
    h_grid = Grids.SpectralElementGrid2D(
        h_topology,
        quad;
        horizontal_layout_type,
        enable_bubble,
    )
    z_topology = Topologies.IntervalTopology(context, z_mesh)
    vertical_grid = FiniteDifferenceGrid(z_topology)
    return ExtrudedFiniteDifferenceGrid(
        h_grid,
        vertical_grid,
        hypsography,
        global_geometry,
    )
end

"""
    SliceXZGrid(
        ::Type{<:AbstractFloat}; # defaults to Float64
        z_elem::Integer,
        x_min::Real,
        x_max::Real,
        z_min::Real,
        z_max::Real,
        periodic_x::Bool,
        n_quad_points::Integer,
        x_elem::Integer,
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        stretch::Meshes.StretchingRule = Meshes.Uniform(),
        hypsography::HypsographyAdaption = Flat(),
        global_geometry::Geometry.AbstractGlobalGeometry = Geometry.CartesianGlobalGeometry(),
        quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    )

A convenience constructor, which builds a
`ExtrudedFiniteDifferenceGrid` with a
`FiniteDifferenceGrid` vertical grid and a
`SpectralElementGrid1D` horizontal grid, given:


 - `FT` the float type (defaults to `Float64`) [`Float32`, `Float64`]
 - `z_elem` the number of z-points
 - `x_min` the domain minimum along the x-direction.
 - `x_max` the domain maximum along the x-direction.
 - `z_min` the domain minimum along the z-direction.
 - `z_max` the domain maximum along the z-direction.
 - `periodic_x` Bool indicating to use periodic domain along x-direction
 - `n_quad_points` the number of quadrature points
 - `x_elem` the number of x-points
 - `device` the [`ClimaComms.device`](@ref)
 - `context` the [`ClimaComms.context`](@ref)
 - `stretch` the mesh [`Meshes.StretchingRule`](@ref) (defaults to `Uniform`) []
 - `hypsography` the hypsography
 - `global_geometry` the global geometry (defaults to `Geometry.CartesianGlobalGeometry`)
 - `quad` the quadrature style (defaults to `Quadratures.GLL{n_quad_points}`)

# Example usage

```julia
using ClimaCore: Grids
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
```
"""
SliceXZGrid(; kwargs...) = SliceXZGrid(Float64; kwargs...)
function SliceXZGrid(
    ::Type{FT};
    z_elem::Integer,
    x_min::Real,
    x_max::Real,
    z_min::Real,
    z_max::Real,
    periodic_x::Bool,
    n_quad_points::Integer,
    x_elem::Integer,
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    stretch::Meshes.StretchingRule = Meshes.Uniform(),
    hypsography::HypsographyAdaption = Flat(),
    global_geometry::Geometry.AbstractGlobalGeometry = Geometry.CartesianGlobalGeometry(),
    quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    horizontal_layout_type = DataLayouts.IFH,
    h_mesh::Meshes.IntervalMesh = DefaultSliceXMesh(
        FT;
        x_min,
        x_max,
        periodic_x,
        x_elem,
    ),
    z_mesh::Meshes.IntervalMesh = DefaultZMesh(
        FT;
        z_min,
        z_max,
        z_elem,
        stretch,
    ),
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)

    h_topology = Topologies.IntervalTopology(context, h_mesh)
    h_grid =
        Grids.SpectralElementGrid1D(h_topology, quad; horizontal_layout_type)
    z_topology = Topologies.IntervalTopology(context, z_mesh)
    vertical_grid = FiniteDifferenceGrid(z_topology)
    return ExtrudedFiniteDifferenceGrid(
        h_grid,
        vertical_grid,
        hypsography,
        global_geometry,
    )
end

"""
    RectangleXYGrid(
        ::Type{<:AbstractFloat}; # defaults to Float64
        x_min::Real,
        x_max::Real,
        y_min::Real,
        y_max::Real,
        periodic_x::Bool,
        periodic_y::Bool,
        n_quad_points::Integer,
        x_elem::Integer, # number of horizontal elements
        y_elem::Integer, # number of horizontal elements
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        hypsography::HypsographyAdaption = Flat(),
        global_geometry::Geometry.AbstractGlobalGeometry = Geometry.CartesianGlobalGeometry(),
        quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    )

A convenience constructor, which builds a
`SpectralElementGrid2D` with a horizontal
`RectilinearMesh` mesh, given:

 - `x_min` the domain minimum along the x-direction.
 - `x_max` the domain maximum along the x-direction.
 - `y_min` the domain minimum along the y-direction.
 - `y_max` the domain maximum along the y-direction.
 - `periodic_x` Bool indicating to use periodic domain along x-direction
 - `periodic_y` Bool indicating to use periodic domain along y-direction
 - `n_quad_points` the number of quadrature points
 - `x_elem` the number of x-points
 - `y_elem` the number of y-points
 - `device` the [`ClimaComms.device`](@ref)
 - `context` the [`ClimaComms.context`](@ref)
 - `hypsography` the hypsography
 - `global_geometry` the global geometry (defaults to `Geometry.CartesianGlobalGeometry`)
 - `quad` the quadrature style (defaults to `Quadratures.GLL{n_quad_points}`)

# Example usage

```julia
using ClimaCore: Grids
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
```
"""
RectangleXYGrid(; kwargs...) = RectangleXYGrid(Float64; kwargs...)
function RectangleXYGrid(
    ::Type{FT};
    x_min::Real,
    x_max::Real,
    y_min::Real,
    y_max::Real,
    periodic_x::Bool,
    periodic_y::Bool,
    n_quad_points::Integer,
    x_elem::Integer, # number of horizontal elements
    y_elem::Integer, # number of horizontal elements
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    hypsography::HypsographyAdaption = Flat(),
    global_geometry::Geometry.AbstractGlobalGeometry = Geometry.CartesianGlobalGeometry(),
    quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    horizontal_layout_type = DataLayouts.IJFH,
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(
        context,
        DefaultRectangleXYMesh(
            FT;
            x_min,
            x_max,
            y_min,
            y_max,
            x_elem,
            y_elem,
            periodic_x,
            periodic_y,
        ),
    ),
    enable_bubble::Bool = false,
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)
    return Grids.SpectralElementGrid2D(
        h_topology,
        quad;
        horizontal_layout_type,
        enable_bubble,
    )
end
