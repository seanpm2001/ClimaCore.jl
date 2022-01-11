push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

using Test
using StaticArrays, IntervalSets, LinearAlgebra, UnPack

import ClimaCore:
    ClimaCore,
    slab,
    Spaces,
    Domains,
    Meshes,
    Geometry,
    Topologies,
    Spaces,
    Fields,
    Operators,
    DataLayouts

using ClimaCore.Geometry

using Logging: global_logger
using TerminalLoggers: TerminalLogger
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33



global_logger(TerminalLogger())


const R = 6.4e6 # radius
const Ω = 7.2921e-5 # Earth rotation (radians / sec)
const z_top = 3.0e4 # height position of the model top
const grav = 9.8 # gravitational constant
const p_0 = 1e5 # mean sea level pressure

const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const T_tri = 273.16 # triple point temperature
const γ = 1.4 # heat capacity ratio
const cv_d = R_d / (γ - 1)
const cp_d = R_d * γ / (γ - 1)
const T_0 = 300 # isothermal atmospheric temperature
const H = R_d * T_0 / grav # scale height
# P = ρ * R_d * T = ρ * R_d * θ * (P / p_0)^(R_d / C_p) ==>
# (P / p_0)^(1 - R_d / C_p) = R_d / p_0 * ρθ ==>
# P = p_0 * (R_d / p_0)^γ * ρθ^γ
const P_ρθ_factor = p_0 * (R_d / p_0)^γ
# P = ρ * R_d * T = ρ * R_d * (ρe_int / ρ / C_v) = (γ - 1) * ρe_int
const P_ρe_factor = γ - 1




# geopotential
gravitational_potential(z) = grav * z
# Π(ρθ) = cp_d * (R_d * ρθ / p_0)^(R_d / cv_d)
pressure(ρθ) = (ρθ*R_d/p_0)^γ * p_0



const hdiv = Operators.Divergence()
const hwdiv = Operators.Divergence()
const hgrad = Operators.Gradient()
const hwgrad = Operators.Gradient()
const hcurl = Operators.Curl()
const hwcurl = Operators.Curl() # Operator.WeakCurl()
const If2c = Operators.InterpolateF2C()
const Ic2f = Operators.InterpolateC2F(
    bottom = Operators.Extrapolate(),
    top = Operators.Extrapolate(),
)
const vdivf2c = Operators.DivergenceF2C(
    top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
)
const vcurlc2f = Operators.CurlC2F(
    bottom = Operators.SetCurl(Geometry.Contravariant12Vector(0.0, 0.0)),
    top = Operators.SetCurl(Geometry.Contravariant12Vector(0.0, 0.0)),
)
const vgradc2f = Operators.GradientC2F(
    bottom = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
    top = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
)

# initial conditions for density and ρθ
function init_sbr_thermo(z)

    p = p_0 * exp(-z / H)

    ρ = p / (R_d * T_0)

    θ = T_0*(p_0/p)^(R_d/cp_d)

    return (ρ = ρ, ρθ = ρ*θ)
end

function rhs!(dY, Y, p, t)
    @unpack P, Φ, ∇Φ = p

    cρ = Y.Yc.ρ # density on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant12Vector on centers
    cρθ = Y.Yc.ρθ # ρθ on centers

    dρ = dY.Yc.ρ
    dw = dY.w
    duₕ = dY.uₕ
    dρθ = dY.Yc.ρθ

    # # 0) update w at the bottom
    # fw = -g^31 cuₕ/ g^33 ????????

    dρ .= 0 .* cρ
    dw .= 0 .* fw
    duₕ .= 0 .* cuₕ
    dρθ .= 0 .* cρθ

    # hyperdiffusion not needed in SBR

    # 1) Mass conservation

    cw = If2c.(fw)
    cuvw = Geometry.Covariant123Vector.(cuₕ) .+ Geometry.Covariant123Vector.(cw)
    # 1.a) horizontal divergence
    dρ .-= hdiv.(cρ .* (cuvw))
    # 1.b) vertical divergence
    # explicit part
    dρ .-= vdivf2c.(Ic2f.(cρ .* cuₕ))
    # implicit part
    dρ .-= vdivf2c.(Ic2f.(cρ) .* fw)

    # 2) Momentum equation
    # curl term
    # effectively a homogeneous Dirichlet condition on u₁ at the boundary

    cω³ = hcurl.(cuₕ) # Contravariant3Vector
    fω¹² = hcurl.(fw) # Contravariant12Vector
    fω¹² .+= vcurlc2f.(cuₕ) # Contravariant12Vector

    # cross product
    # convert to contravariant
    # these will need to be modified with topography
    fu¹² =
        Geometry.Contravariant12Vector.(
            Geometry.Covariant123Vector.(Ic2f.(cuₕ)),
        ) # Contravariant12Vector in 3D
    fu³ = Geometry.Contravariant3Vector.(Geometry.Covariant123Vector.(fw))
    @. duₕ -= If2c(fω¹² × fu³)
    # Needed for 3D:
    @. duₕ -=
        (f + cω³) ×
        Geometry.Contravariant12Vector(Geometry.Covariant123Vector(cuₕ))
    cp = @. pressure(cρθ)
    @. duₕ -= hgrad(cp) / cρ
    cK = @. (norm(cuvw)^2) / 2
    @. duₕ -= hgrad(cK + Φ)

    @. dw -= fω¹² × fu¹² # Covariant3Vector on faces
    @. dw -= vgradc2f(cp) / Ic2f(cρ)
    @. dw -= (vgradc2f(cK) + ∇Φ)

    # 3) ρθ
    @. dρθ -= hdiv(cuvw * cρθ)
    @. dρθ -= vdivf2c(fw * Ic2f(cρθ))
    @. dρθ -= vdivf2c(Ic2f(cuₕ * cρθ))

    Spaces.weighted_dss!(dY.Yc)
    Spaces.weighted_dss!(dY.uₕ)
    Spaces.weighted_dss!(dY.w)

    return dY
end




function rhs_remainder!(dY, Y, p, t)
    # @info "Remainder part"
    @unpack P, Φ, ∇Φ = p

    cρ = Y.Yc.ρ # density on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant12Vector on centers
    cρθ = Y.Yc.ρθ # ρθ on centers

    dρ = dY.Yc.ρ
    dw = dY.w
    duₕ = dY.uₕ
    dρθ = dY.Yc.ρθ


    # uₕ_phy = Geometry.transform.(Ref(Geometry.UVAxis()), cuₕ)
    # w_phy = Geometry.transform.(Ref(Geometry.WAxis()), fw)
    # @info "maximum vertical velocity is w, u_h", maximum(abs.(w_phy.components.data.:1)), maximum(abs.(uₕ_phy.components.data.:1)), maximum(abs.(uₕ_phy.components.data.:2))


    # # 0) update w at the bottom
    # fw = -g^31 cuₕ/ g^33 ????????



    dρ .= 0 .* cρ
    dw .= 0 .* fw
    duₕ .= 0 .* cuₕ
    dρθ .= 0 .* cρθ

    # hyperdiffusion not needed in SBR

    # 1) Mass conservation
    cw = If2c.(fw)
    cuvw = Geometry.Covariant123Vector.(cuₕ) .+ Geometry.Covariant123Vector.(cw)

    # 1.a) horizontal divergence
    dρ .-= hdiv.(cρ .* (cuvw))
    dρ .-= vdivf2c.(Ic2f.(cρ .* cuₕ))
    
    # 2) Momentum equation

    # curl term
    # effectively a homogeneous Dirichlet condition on u₁ at the boundary
    cω³ = hcurl.(cuₕ) # Contravariant3Vector
    fω¹² = hcurl.(fw) # Contravariant12Vector
    fω¹² .+= vcurlc2f.(cuₕ) # Contravariant12Vector

    # cross product
    # convert to contravariant
    # these will need to be modified with topography
    fu¹² =
        Geometry.Contravariant12Vector.(
            Geometry.Covariant123Vector.(Ic2f.(cuₕ)),
        ) # Contravariant12Vector in 3D
    fu³ = Geometry.Contravariant3Vector.(Geometry.Covariant123Vector.(fw))
    
    @. duₕ -= If2c(fω¹² × fu³)
    # Needed for 3D:
    @. duₕ -=
        (f + cω³) ×
        Geometry.Contravariant12Vector(Geometry.Covariant123Vector(cuₕ))
    cp = @. pressure(cρθ)
    @. duₕ -= hgrad(cp) / cρ
    cK = @. (norm(cuvw)^2) / 2
    @. duₕ -= hgrad(cK + Φ)

    
    @. dw -= fω¹² × fu¹² # Covariant3Vector on faces
    
    @. dw -= vgradc2f(cK)
    


    # 3) ρθ
    @. dρθ -= hdiv(cuvw * cρθ)
    @. dρθ -= vdivf2c(Ic2f(cuₕ * cρθ))


    Spaces.weighted_dss!(dY.Yc)
    Spaces.weighted_dss!(dY.uₕ)
    Spaces.weighted_dss!(dY.w)

    return dY
end



function rhs_implicit!(dY, Y, p, t)
    # @info "Implicit part"
    @unpack P, Φ, ∇Φ = p

    cρ = Y.Yc.ρ # density on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant12Vector on centers
    cρθ = Y.Yc.ρθ # ρθ on centers

    dρ = dY.Yc.ρ
    dw = dY.w
    duₕ = dY.uₕ
    dρθ = dY.Yc.ρθ

    # # 0) update w at the bottom
    # fw = -g^31 cuₕ/ g^33 ????????

    dρ .= 0 .* cρ
    dw .= 0 .* fw
    duₕ .= 0 .* cuₕ
    dρθ .= 0 .* cρθ

    # hyperdiffusion not needed in SBR

    # 1) Mass conservation

    

    # 1.b) vertical divergence
    # we want the total u³ at the boundary to be zero: we can either constrain
    # both to be zero, or allow one to be non-zero and set the other to be its
    # negation

    # TODO implicit
    dρ .-= vdivf2c.(Ic2f.(cρ) .* fw)

    # 2) Momentum equation

    cp = @. pressure(cρθ)
    @. dw -= vgradc2f(cp) / Ic2f(cρ)
    # @. dw -= R_d/cv_d * Ic2f(Π(cρθ)) * vgradc2f(cρθ) / Ic2f(cρ)
    @. dw -= ∇Φ

    # 3) ρθ
    @. dρθ -= vdivf2c(fw * Ic2f(cρθ))

    return dY
end



include("solid_body_rotation_3d_implicit.jl")





