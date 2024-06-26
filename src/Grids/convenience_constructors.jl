import ..DataLayouts
#=
These convenience constructors accept integer
keyword inputs, so they are dynamically created. You may
want to use a different constructor if you're making the
object in a performance-critical section, and if you know
the type parameters at compile time.

If no convenience constructor exists to make the
grid you need, then you may need to use our lower
level compose-able API.
=#
check_device_context(context, device) =
    @assert ClimaComms.device(context) == device "The given device and context device do not match."

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
        h_mesh = Meshes.EquiangularCubedSphere(Domains.SphereDomain(radius), h_elem),
    )

A convenience constructor, which builds a
`ExtrudedFiniteDifferenceGrid`.
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
    horizontal_layout_type = DataLayouts.IJFH,
    h_mesh = Meshes.EquiangularCubedSphere(
        Domains.SphereDomain{FT}(radius),
        h_elem,
    ),
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(context, h_mesh),
    enable_bubble::Bool = false
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)

    z_boundary_names = (:bottom, :top)
    h_grid =
        Grids.SpectralElementGrid2D(h_topology, quad; horizontal_layout_type, enable_bubble)
    z_domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(z_min),
        Geometry.ZPoint{FT}(z_max);
        boundary_names = z_boundary_names,
    )
    z_mesh = Meshes.IntervalMesh(z_domain, stretch; nelems = z_elem)
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
        n_quad_points::Integer,
        h_elem::Integer,
        device::ClimaComms.AbstractDevice = ClimaComms.device(),
        context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
        quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
        h_mesh = Meshes.EquiangularCubedSphere(Domains.SphereDomain(radius), h_elem),
    )

A convenience constructor, which builds a
`SpectralElementGrid2D`.
"""
CubedSphereGrid(; kwargs...) = CubedSphereGrid(Float64; kwargs...)
function CubedSphereGrid(
    ::Type{FT};
    radius::Real,
    n_quad_points::Integer,
    h_elem::Integer,
    device::ClimaComms.AbstractDevice = ClimaComms.device(),
    context::ClimaComms.AbstractCommsContext = ClimaComms.context(device),
    quad::Quadratures.QuadratureStyle = Quadratures.GLL{n_quad_points}(),
    h_mesh = Meshes.EquiangularCubedSphere(
        Domains.SphereDomain{FT}(radius),
        h_elem,
    ),
    horizontal_layout_type = DataLayouts.IJFH,
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(context, h_mesh),
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
) where {FT}
    check_device_context(context, device)
    z_boundary_names = (:bottom, :top)
    z_domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(z_min),
        Geometry.ZPoint{FT}(z_max);
        boundary_names = z_boundary_names,
    )
    z_mesh = Meshes.IntervalMesh(z_domain, stretch; nelems = z_elem)
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
    )

A convenience constructor, which builds a
`ExtrudedFiniteDifferenceGrid` with a
`FiniteDifferenceGrid` vertical grid and a
`SpectralElementGrid2D` horizontal grid.
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
    horizontal_layout_type = DataLayouts.IJFH,
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(context, RectangleXYMesh(FT; x_min, x_max, y_min, y_max, x_elem, y_elem, periodic_x, periodic_y),),
    enable_bubble::Bool = false,
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)
    h_grid =
        Grids.SpectralElementGrid2D(h_topology, quad; horizontal_layout_type, enable_bubble)
    z_domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(z_min),
        Geometry.ZPoint{FT}(z_max);
        boundary_names = (:bottom, :top),
    )
    z_mesh = Meshes.IntervalMesh(z_domain, stretch; nelems = z_elem)
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
`SpectralElementGrid1D` horizontal grid.
 - ``
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
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)

    x1boundary = (:east, :west)
    z_boundary_names = (:bottom, :top)
    h_domain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(x_min),
        Geometry.XPoint{FT}(x_max);
        periodic = periodic_x,
        boundary_names = x1boundary,
    )
    h_mesh = Meshes.IntervalMesh(h_domain; nelems = x_elem)
    h_topology = Topologies.IntervalTopology(context, h_mesh)
    h_grid =
        Grids.SpectralElementGrid1D(h_topology, quad; horizontal_layout_type)
    z_domain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(z_min),
        Geometry.ZPoint{FT}(z_max);
        boundary_names = z_boundary_names,
    )
    z_mesh = Meshes.IntervalMesh(z_domain, stretch; nelems = z_elem)
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
    RectangleXYMesh(
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
RectangleXYMesh(; kwargs...) = RectangleXYMesh(Float64; kwargs...)
function RectangleXYMesh(::Type{FT};
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
`RectilinearMesh` mesh.
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
    h_topology::Topologies.AbstractDistributedTopology = Topologies.Topology2D(context, RectangleXYMesh(FT; x_min, x_max, y_min, y_max, x_elem, y_elem, periodic_x, periodic_y)),
    enable_bubble::Bool = false,
) where {FT}
    @assert horizontal_layout_type <: DataLayouts.AbstractData
    check_device_context(context, device)
    return Grids.SpectralElementGrid2D(h_topology, quad; horizontal_layout_type, enable_bubble)
end
