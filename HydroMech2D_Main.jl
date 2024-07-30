# HydroMech2D_Main run simulation and plot visualization.
# The algorithm consists of an outer loop that simulates the progression of time in dt increments and
# an inner loop that runs an approximation algorithm that is optimized using parallelization via the
# ParallelStencil module. 
const USE_GPU = false # Use GPU? If this is set false, then no GPU needs to be available
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
using DelimitedFiles
using DataFrames
using Distributions
using Parameters
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 2)
else
    @init_parallel_stencil(Threads, Float64, 2)
end
using Plots, Printf, Statistics, LinearAlgebra

include("HydroMechFunctions.jl")

##################################################
@views function HydroMech2D()
    constants = HydroMechConstants()
    variables = HydroMechVariables()
    initialize_variables!(variables, constants)
    @unpack dx, dy, ly, dt, ε, iterMax, nout, nx, ny, dt_red, lx, ly, res, t_tot, nperm, μs0, μs = constants
    @unpack dτPf, Vx, Vy, qDx, qDy, Ry, RPf, ∇V, qDys, Phi, Pt, Pf = variables
    # Preparation of visualisation
    ENV["GKSwstype"]="nul";
    if isdir(string(nperm) * "_poro_AdjDepth_viz2D_out_16")==false mkdir(string(nperm) * "_poro_AdjDepth_viz2D_out_16") end; loadpath = "./" * string(nperm) * "_poro_AdjDepth_viz2D_out_16" * "/"; anim = Animation(loadpath,String[])
    println("Animation directory: $(anim.dir)")
    X, Y, Yv = 0:dx:lx, 0:dy:ly, (-dy/2):dy:(ly+dy/2)
    # Time loop
    times = []
    count = 0
    t        = 0.0  # t gives the "time" of the outer (time) loop
    it       = 1  # it gives the number of iterations of the outer (time) loop
    while t< t_tot
        count = count + 1
        update_old!(variables)
        err=2*ε; iter=1; niter=0  # iter/niter give iterations of the inner (err) loop
        while err > ε && iter <= iterMax
            if (iter==11)  global wtime0 = Base.time()  end
            compute_params_∇!(variables, constants)
            compute_RP!(variables, constants)
            @parallel (1:size(dτPf,1), 1:size(dτPf,2))  bc_x!(dτPf)
            @parallel (1:size(dτPf,1), 1:size(dτPf,2))  bc_y!(dτPf)
            compute_P_τ!(variables, constants)
            compute_res!(variables, constants)
            compute_update!(variables, constants)
            @parallel (1:size(Vx,1),  1:size(Vx,2))  bc_y!(Vx)
            @parallel (1:size(Vy,1),  1:size(Vy,2))  bc_x!(Vy)
            @parallel (1:size(qDx,1), 1:size(qDx,2)) bc_y!(qDx)
            @parallel (1:size(qDy,1), 1:size(qDy,2)) bc_x!(qDy)
            if mod(iter,nout)==0
                global norm_Ry, norm_RPf
                norm_Ry = norm(Ry)/length(Ry); norm_RPf = norm(RPf)/length(RPf); err = max(norm_Ry, norm_RPf)
                @printf("iter = %d, err = %1.3e [norm_Ry=%1.3e, norm_RPf=%1.3e] \n", iter, err, norm_Ry, norm_RPf)
            end
            iter+=1; niter+=1
        end
        # Performance
        wtime    = Base.time()-wtime0
        A_eff    = (8*2)/1e9*nx*ny*sizeof(Data.Number)  # Effective main memory access per iteration [GB] (Lower bound of required memory access: Te has to be read and written: 2 whole-array memaccess; Ci has to be read: : 1 whole-array memaccess)
        wtime_it = wtime/(niter-10)                     # Execution time per iteration [s]
        T_eff    = A_eff/wtime_it                       # Effective memory throughput [GB/s]
        @printf("it = %d, time = %1.3e sec (@ T_eff = %1.2f GB/s) \n", it, wtime, round(T_eff, sigdigits=2))
	    # Time
	    dt = dt_red/(1e-10+maximum(abs.(∇V)))
        t  = t + dt
        push!(times,t)
        it+=1
        # Visualisation
        default(size=(800,700))
        if mod(it,5)==0
            #p1 = heatmap(X, Y,  Array(Phi)'  , aspect_ratio=1, xlims=(X[1],X[end]), ylims=(Y[1],Y[end]), c=:viridis, title="Porosity")
            #p2 = heatmap(X, Y,  Array(Pf)', aspect_ratio=1, xlims=(X[1],X[end]), ylims=(Y[1],Y[end]), c=:viridis, xlabel="L_x",  ylabel="L_y", title="Fluid Pressure, time= " * string(round(t, digits = 3)))
            p3 = heatmap(X, Yv, Array(qDy)'  , aspect_ratio=1, xlims=(X[1],X[end]), ylims=(Yv[1],Yv[end]), c=:viridis, xlabel="L_x",  ylabel="L_y", title="Vertical Darcy Flux, time= " * string(round(t, digits = 3)))
            #p4 = heatmap(X, Yv, Array(Vy)'   , aspect_ratio=1, xlims=(X[1],X[end]), ylims=(Yv[1],Yv[end]), c=:viridis, title="vertical velocity")
            display(plot(p3));
            frame(anim);
        end
    end
    save_folder = "Darcy"
    path_name = "/home/users/kmagno/pockmark-degassing/For_Quals/" * save_folder * "/"
    fnm2 = string(nperm)
    fnm3 = string(t_tot)
    fnm4 = string(res)
    fnm_qDys = path_name * "qDys_nperm" * fnm2 * "_ttot" * fnm3 * "_res" * fnm4 * ".csv"
    writedlm(fnm_qDys, Array(qDys), ',')
    fnm_times = path_name * "times_nperm" * fnm2 * "_ttot" * fnm3 * "_res" * fnm4 * ".csv"
    fnm_gif = path_name * "HydroMech2D_nperm" * fnm2 * "_ttot" * fnm3 * "_res" * fnm4 * ".gif"
    gif(anim, fnm_gif, fps = 5)
    return
end

println("hello world")
HydroMech2D()
println("done")
