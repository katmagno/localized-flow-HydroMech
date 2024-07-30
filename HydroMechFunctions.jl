# HydroMechFunctions.jl
using Parameters
using ParallelStencil
include("HydroMechDataTypes.jl")

function my_1D_detrend(input)
    x       = collect(1:length(input))
    x_mean  = sum(x)/length(x)
    in_mean = sum(input)/length(input)
    slope   = sum((x .- x_mean).*(input .- in_mean))/sum((x .- x_mean).^2)
    inter   = in_mean - slope*x_mean
    lin_fit = slope.*x .+ inter
    output  = input .- lin_fit
    return output
end

function update_old!(variables::HydroMechVariables)
    @unpack Phi, Phi_o, ∇V, ∇V_o = variables
    @parallel update_old!(Phi_o, ∇V_o, Phi, ∇V)
    @pack! variables = Phi, Phi_o, ∇V, ∇V_o
    return
end

@parallel function update_old!(Phi_o::Data.Array, ∇V_o::Data.Array, Phi::Data.Array, ∇V::Data.Array)
    @all(Phi_o) = @all(Phi)
    @all(∇V_o)  = @all(∇V)
    return
end

function compute_params_∇!(variables::HydroMechVariables, constants::HydroMechConstants)
    @unpack μs, η2μs, R, λPe, k_μf0, ϕ0, nperm, θ_e, θ_k, ρfg, ρsg, ρgBG, dx, dy = constants
    @unpack EtaC, K_muf, Rog, ∇V, ∇qD, Phi, Pf, Pt, Vx, Vy, qDx, qDy  = variables
    @parallel compute_params_∇!(EtaC, K_muf, Rog, ∇V, ∇qD, Phi, Pf, Pt, Vx, Vy, qDx, qDy, μs, η2μs, R, λPe, k_μf0, ϕ0, nperm, θ_e, θ_k, ρfg, ρsg, ρgBG, dx, dy)
    @pack! variables = EtaC, K_muf, Rog, ∇V, ∇qD, Phi, Pf, Pt, Vx, Vy, qDx, qDy
    return
end

@parallel function compute_params_∇!(EtaC::Data.Array, K_muf::Data.Array, Rog::Data.Array, ∇V::Data.Array, ∇qD::Data.Array, Phi::Data.Array, Pf::Data.Array, Pt::Data.Array, Vx::Data.Array, Vy::Data.Array, qDx::Data.Array, qDy::Data.Array, μs::Data.Number, η2μs::Data.Number, R::Data.Number, λPe::Data.Number, k_μf0::Data.Number, ϕ0::Data.Number, nperm::Data.Number, θ_e::Data.Number, θ_k::Data.Number, ρfg::Data.Number, ρsg::Data.Number, ρgBG::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(EtaC)  = (1.0-θ_e)*@all(EtaC)  + θ_e*( μs/@all(Phi)*η2μs*(1.0+0.5*(1.0/R-1.0)*(1.0+tanh((@all(Pf)-@all(Pt))/λPe))) )
    @all(K_muf) = (1.0-θ_k)*@all(K_muf) + θ_k*( k_μf0 * (@all(Phi)/0.02)^nperm)
    #^nperm)
    @all(Rog)   = ρfg*@all(Phi) + ρsg*(1.0-@all(Phi)) - ρgBG
    @all(∇V)    = @d_xa(Vx)/dx  + @d_ya(Vy)/dy
    @all(∇qD)   = @d_xa(qDx)/dx + @d_ya(qDy)/dy
    return
end

function compute_RP!(variables::HydroMechVariables, constants::HydroMechConstants)
    @unpack Pfsc, Pfdmp, min_dxy2, dx, dy = constants
    @unpack dτPf, RPt, RPf, K_muf, ∇V, ∇qD, ∇qD, Pt, Pf, EtaC, Phi = variables
    @parallel compute_RP!(dτPf, RPt, RPf, K_muf, ∇V, ∇qD, Pt, Pf, EtaC, Phi, Pfsc, Pfdmp, min_dxy2, dx, dy)
    @pack! variables = dτPf, RPt, RPf, K_muf, ∇V, ∇qD, ∇qD, Pt, Pf, EtaC, Phi
    return
end

@parallel function compute_RP!(dτPf::Data.Array, RPt::Data.Array, RPf::Data.Array, K_muf::Data.Array, ∇V::Data.Array, ∇qD::Data.Array, Pt::Data.Array, Pf::Data.Array, EtaC::Data.Array, Phi::Data.Array, Pfsc::Data.Number, Pfdmp::Data.Number, min_dxy2::Data.Number, dx::Data.Number, dy::Data.Number)
    @inn(dτPf) = min_dxy2/@maxloc(K_muf)/4.1/Pfsc
    @all(RPt)  =                 - @all(∇V)  - (@all(Pt) - @all(Pf))/(@all(EtaC)*(1.0-@all(Phi)))
    @all(RPf)  = @all(RPf)*Pfdmp - @all(∇qD) + (@all(Pt) - @all(Pf))/(@all(EtaC)*(1.0-@all(Phi)))
    return
end

function compute_P_τ!(variables::HydroMechVariables, constants::HydroMechConstants)
    @unpack μs, β_n, dx, dy, dτPt = constants
    @unpack Pt, Pf, τxx, τyy, σxy, RPt, RPf, dτPf, Vx, Vy, ∇V = variables
    @parallel compute_P_τ!(Pt, Pf, τxx, τyy, σxy, RPt, RPf, dτPf, Vx, Vy, ∇V, dτPt, μs, β_n, dx, dy)
    @pack! variables = Pt, Pf, τxx, τyy, σxy, RPt, RPf, dτPf, Vx, Vy, ∇V
    return
end

@parallel function compute_P_τ!(Pt::Data.Array, Pf::Data.Array, τxx::Data.Array, τyy::Data.Array, σxy::Data.Array, RPt::Data.Array, RPf::Data.Array, dτPf::Data.Array, Vx::Data.Array, Vy::Data.Array, ∇V::Data.Array, dτPt::Data.Number, μs::Data.Number, β_n::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(Pt)  = @all(Pt) +      dτPt *@all(RPt)
    @all(Pf)  = @all(Pf) + @all(dτPf)*@all(RPf)
    @all(τxx) = 2.0*μs*( @d_xa(Vx)/dx - 1.0/3.0*@all(∇V) - β_n*@all(RPt) )
    @all(τyy) = 2.0*μs*( @d_ya(Vy)/dy - 1.0/3.0*@all(∇V) - β_n*@all(RPt) )
    @all(σxy) = 2.0*μs*(0.5*( @d_yi(Vx)/dy + @d_xi(Vy)/dx ))
    return
end

function compute_res!(variables::HydroMechVariables, constants::HydroMechConstants)
    @unpack dampX, dampY, dx, dy = constants
    @unpack Rx, Ry, dVxdτ, dVydτ, τxx, τyy, σxy, Pt, Rog = variables
    @parallel compute_res!(Rx, Ry, dVxdτ, dVydτ, τxx, τyy, σxy, Pt, Rog, dampX, dampY, dx, dy)
    @pack! variables = Rx, Ry, dVxdτ, dVydτ, τxx, τyy, σxy, Pt, Rog
    return
end

@parallel function compute_res!(Rx::Data.Array, Ry::Data.Array, dVxdτ::Data.Array, dVydτ::Data.Array, τxx::Data.Array, τyy::Data.Array, σxy::Data.Array, Pt::Data.Array, Rog::Data.Array, dampX::Data.Number, dampY::Data.Number, dx::Data.Number, dy::Data.Number)
    @all(Rx)    = @d_xi(τxx)/dx + @d_ya(σxy)/dy - @d_xi(Pt)/dx
    @all(Ry)    = @d_yi(τyy)/dy + @d_xa(σxy)/dx - @d_yi(Pt)/dy - @av_yi(Rog)
    @all(dVxdτ) = dampX*@all(dVxdτ) + @all(Rx)
    @all(dVydτ) = dampY*@all(dVydτ) + @all(Ry)
    return
end

function compute_update!(variables::HydroMechVariables, constants::HydroMechConstants)
    @unpack dτV, ρfg, ρgBG, CN, dt, dx, dy = constants
    @unpack Vx, Vy, qDx, qDy, Phi, dVxdτ, dVydτ, K_muf, Pf, Phi_o, ∇V, ∇V_o = variables
    @parallel compute_update!(Vx, Vy, qDx, qDy, Phi, dVxdτ, dVydτ, K_muf, Pf, Phi_o, ∇V, ∇V_o, dτV, ρfg, ρgBG, CN, dt, dx, dy)
    @pack! variables = Vx, Vy, qDx, qDy, Phi, dVxdτ, dVydτ, K_muf, Pf, Phi_o, ∇V, ∇V_o
    return
end

@parallel function compute_update!(Vx::Data.Array, Vy::Data.Array, qDx::Data.Array, qDy::Data.Array, Phi::Data.Array, dVxdτ::Data.Array, dVydτ::Data.Array, K_muf::Data.Array, Pf::Data.Array, Phi_o::Data.Array, ∇V::Data.Array, ∇V_o::Data.Array, dτV::Data.Number, ρfg::Data.Number, ρgBG::Data.Number, CN::Data.Number, dt::Data.Number, dx::Data.Number, dy::Data.Number)
    @inn(Vx)  =  @inn(Vx) + dτV*@all(dVxdτ)
    @inn(Vy)  =  @inn(Vy) + dτV*@all(dVydτ)
    @inn(qDx) = -@av_xi(K_muf)*(@d_xi(Pf)/dx)
    @inn(qDy) = -@av_yi(K_muf)*(@d_yi(Pf)/dy + (ρfg - ρgBG))
    @all(Phi) =  @all(Phi_o) + (1.0-@all(Phi))*(CN*@all(∇V_o) + (1.0-CN)*@all(∇V))*dt
    return
end

@parallel_indices (ix,iy) function bc_x!(A::Data.Array)
    A[1  , iy] = A[2    , iy]
    A[end, iy] = A[end-1, iy]
    return
end

@parallel_indices (ix,iy) function bc_y!(A::Data.Array)
    A[ix, 1  ] = A[ix, 2    ]
    A[ix, end] = A[ix, end-1]
    return
end