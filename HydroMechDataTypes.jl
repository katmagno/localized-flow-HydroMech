# HydroMechDataTypes.jl
using Parameters
using ParallelStencil

@with_kw struct HydroMechConstants
    # Physics - scales
    ns       = 1.0
    ρfg      = 1.0             # fluid rho*g
    k_μf0    = 1.0             # reference permeability
    ηC0      = 1.0             # reference bulk viscosity
    # Physics - non-dimensional parameters
    η2μs     = 10.0            # bulk/shear viscosity ration 10
    R        = 800.0           # Compaction/decompaction strength ratio for bulk rheology 800
    nperm    = 2.7         # Carman-Kozeny exponent
    ϕ0       = 0.02            # reference porosity
    ra       = 2               # radius of initil porosity perturbation
    λ0       = 1.0             # standard deviation of initial porosity perturbation
    t_tot    = 0.5     # total time
    # Physics - dependent scales
    ρsg      = 2.0*ρfg         # solid rho*g 2.0
    lx       = 35              # domain size x
    ly       = ra*lx           # domain size y
    ϕA       = 2*ϕ0          # amplitude of initial porosity perturbation
    λPe      = 0.01            # effective pressure transition zone
    dt       = 1e-5            # physical time-step      ### GOOD
    # Numerics
    CN       = 0.5             # Crank-Nicolson CN=0.5, Backward Euler CN=0.0
    res      = 256             # 256
    nx       = res-1           # numerical grid resolutions; should be a mulitple of 32-1 for optimal GPU perf
    ny       = ra*res-1        # numerical grid resolutions; should be a mulitple of 32-1 for optimal GPU perf
    ε        = 1e-5            # non-linear tolerance       # BAD
    iterMax  = 5e3             # max nonlinear iterations
    nout     = 256             # error checking frequency
    β_n      = 1.0             # numerical compressibility
    Vdmp     = 5.0             # velocity damping for momentum equations
    Pfdmp    = 0.8             # fluid pressure damping for momentum equations
    Vsc      = 2.0             # reduction of PT steps for velocity
    Ptsc     = 2.0             # reduction of PT steps for total pressure    ### BAD
    Pfsc     = 4.0             # reduction of PT steps for fluid pressure
    θ_e      = 9e-1            # relaxation factor for non-linear viscosity
    θ_k      = 1e-1            # relaxation factor for non-linear permeability
    dt_red   = 1e-3            # reduction of physical timestep
    # Derived physics
    μs       = ηC0*ϕ0/η2μs                       # solid shear viscosity ###########################################################################################
    μs0       = ηC0*ϕ0/η2μs 
    λ        = λ0*sqrt(k_μf0*ηC0)                # initial perturbation width
    ρgBG     = ρfg*ϕ0 + ρsg*(1.0-ϕ0)             # Background density
    # Derived numerics
    dx       = lx/(nx-1)                        # grid step in x, y
    dy       = ly/(ny-1)                        # grid step in x, y
    min_dxy2 = min(dx,dy)^2
    dτV      = min_dxy2/μs/(1.0+β_n)/4.1/Vsc     # PT time step for velocity
    dτPt     = 4.1*μs*(1.0+β_n)/max(nx,ny)/Ptsc
    dampX    = 1.0-Vdmp/nx
    dampY    = 1.0-Vdmp/ny
    permeability = 5e-17
end

@with_kw mutable struct HydroMechVariables
    # Array allocations
    Phi_o = 0.0 # Overridden in initialize_variables
    Pt = 0.0 # Overridden in initialize_variables
    Pf = 0.0 # Overridden in initialize_variables
    Rog = 0.0 # Overridden in initialize_variables
    ∇V = 0.0 # Overridden in initialize_variables
    ∇V_o = 0.0 # Overridden in initialize_variables
    ∇qD = 0.0 # Overridden in initialize_variables
    dτPf = 0.0 # Overridden in initialize_variables
    RPt = 0.0 # Overridden in initialize_variables
    RPf = 0.0 # Overridden in initialize_variables
    τxx = 0.0 # Overridden in initialize_variables
    τyy = 0.0 # Overridden in initialize_variables
    σxy = 0.0 # Overridden in initialize_variables
    dVxdτ = 0.0 # Overridden in initialize_variables
    dVydτ = 0.0 # Overridden in initialize_variables
    Rx = 0.0 # Overridden in initialize_variables
    Ry = 0.0 # Overridden in initialize_variables
    Vx = 0.0 # Overridden in initialize_variables
    Vy = 0.0 # Overridden in initialize_variables
    qDx = 0.0 # Overridden in initialize_variables
    # Initial conditions
    qDy = 0.0 # Overridden in initialize_variables
    qDys = 0.0 # Overridden in initialize_variables
    topA = 0.0 # Overridden in initialize_variables
    top = 0.0 # Overridden in initialize_variables
    #Phi_z
    Phi = 0.0 # Overridden in initialize_variables
    #Radc
    EtaC = 0.0 # Overridden in initialize_variables
    K_muf = 0.0 # Overridden in initialize_variables
    ϕ0bc = 0.0 # Overridden in initialize_variables
end

function initialize_variables!(variables::HydroMechVariables, constants::HydroMechConstants)
    @unpack ns, nx, ny, t_tot, ly, dy, ϕA, η2μs, μs, k_μf0, ϕ0, nperm, ρsg, ρfg, μs, μs0 = constants
    @unpack Phi_o, Pt, Pf, Rog, ∇V, ∇V_o, ∇qD, dτPf, RPt, RPf, τxx, τyy, σxy, dVxdτ, dVydτ, Rx, Ry, Vx, Vy, qDx, qDy, qDys, topA, top, Phi, EtaC, K_muf, ϕ0bc = variables
    # Array allocations
    Phi_o    = @zeros(nx  ,ny  )
    Pt       = @zeros(nx  ,ny  )
    Pf       = @zeros(nx  ,ny  )
    Rog      = @zeros(nx  ,ny  )
    ∇V       = @zeros(nx  ,ny  )
    ∇V_o     = @zeros(nx  ,ny  )
    ∇qD      = @zeros(nx  ,ny  )
    dτPf     = @zeros(nx  ,ny  )
    RPt      = @zeros(nx  ,ny  )
    RPf      = @zeros(nx  ,ny  )
    τxx      = @zeros(nx  ,ny  )
    τyy      = @zeros(nx  ,ny  )
    σxy      = @zeros(nx-1,ny-1)
    dVxdτ    = @zeros(nx-1,ny-2)
    dVydτ    = @zeros(nx-2,ny-1)
    Rx       = @zeros(nx-1,ny-2)
    Ry       = @zeros(nx-2,ny-1)
    Vx       = @zeros(nx+1,ny  )
    Vy       = @zeros(nx  ,ny+1)
    qDx      = @zeros(nx+1,ny  )
    # Initial conditions
    qDy      =   zeros(nx  ,ny+1)
    qDys     =   zeros(round(Int, (t_tot/1e-6))  ,nx)
    topA     =   cumsum(rand(nx,1), dims = 1)
    top      =   my_1D_detrend(topA)
    top      =   top .- minimum(top)
    top      =   top./maximum(top)*ly/10*7 # Default is 7 for heterogeneous reservoir-overlying sediment interface
    Phi      =   ϕ0*ones(nx  ,ny)
    μs       =   μs0
    for iy in 1:ny
        for ix in 5:nx-5
            if iy>floor(Int, 0.01*ly/dy) && iy<floor(Int, 0.1*ly/dy + top[ix])
                Phi[ix,iy] = Phi[ix,iy] + ϕA
                Phi[ix,iy] + ϕA
            end
        end
    end
    Phi[[1 end],:] = Phi[[2 end-1],:]
    Phi[:,[1 end]] = Phi[:,[2 end-1]];
    Phi = abs.(Phi)
    EtaC     = μs./Phi.*η2μs
    K_muf    = k_μf0.*(Phi./ϕ0)
    ϕ0bc     = mean.(Phi[:,end])
    qDy[:,[1,end]] .= (ρsg.-ρfg).*(1.0.-ϕ0bc).*k_μf0.*(ϕ0bc./ϕ0).^nperm
    Phi      = Data.Array(Phi)
    EtaC     = Data.Array(EtaC)
    K_muf    = Data.Array(K_muf)
    qDy      = Data.Array(qDy)
    qDys     = Data.Array(qDys)
    @pack! variables = Phi_o, Pt, Pf, Rog, ∇V, ∇V_o, ∇qD, dτPf, RPt, RPf, τxx, τyy, σxy, dVxdτ, dVydτ, Rx, Ry, Vx, Vy, qDx, qDy, qDys, topA, top, Phi, EtaC, K_muf, ϕ0bc
end
