#added on laptop... 3/30
struct OptControlProb
    locs::Matrix{Float64}
    C::SparseMatrixCSC{Int64}
    Z::SparseMatrixCSC{Int64}
    F::SparseMatrixCSC{Bool, Int64}
    GFlipped::SparseMatrixCSC{Bool, Int64}
    Bmax::Float64
    Q0::Float64
    B0::Float64
    tag::String
end
struct OptControlSolution
    locs::Matrix{Float64}
    C::SparseMatrixCSC{Int64}
    Z::SparseMatrixCSC{Int64}
    F::SparseMatrixCSC{Bool, Int64}
    path::Vector{Int64} #this remains Int until we do OptControl on the trajectory as well...
    gen::Vector{Float64} #needs to be float becuase the gen can be real number from 0 to 1
    timevec::Vector{Float64}
    batt_state::Vector{Float64}
    fuel_state::Vector{Float64}
    noise_restrictions::Vector{Bool}
    tag::String
    time_to_solve::Float64
end

"""
Give an OCP::OptControlProb, path and gen

Output:
out::OptControlSolution
"""
function solve_gen_optimal_control(OCP::OptControlProb, path::Vector{Int64}, gen_MILP::Vector{Bool}, N::Int64; linear = true)
    myModel = Model(Ipopt.Optimizer)
    set_silent(myModel)
    
    # N = length(path) # make u[1] = 0... but need length(path) for indexing wrt path
    tf = 1 #assume path takes 1 second
    dt = (tf)/(N-1)
    
    #pull from OCP
    C,Z,GFlipped = OCP.C, OCP.Z, OCP.GFlipped
    Cvec = C_time_discretized(C, path, N)
    Zvec = C_time_discretized(Z, path, N)
    g_min,g_max = 0, OCP.Q0     #Bounds on the States
    b_min, b_max = 6, 99.9
    
    u_min, u_max = 0, 1     # Bounds on Control

    noiseR_along_path = Float64.([!OCP.GFlipped[path[idx-1],path[idx]] for idx in 2:length(path)])
    noise_intervals = noiseR_to_interval_sets(noiseR_along_path)

    @JuMP.variables(myModel, begin
        # dt ≥ 0, (start = dt_guess) # Delta T for time steps
        g_min  ≤ g[1:N]  ≤ g_max    # fuel level, fuel is monotomically decreasing
        b_min  ≤ b[1:N] ≤ b_max      # battery state of charge
        u_min  ≤ u[1:N] ≤ u_max      # generator setting
    end)

    fix(b[1], b_max; force = true)  
    fix(g[1], OCP.Q0; force = true)  
    fix(u[1], 0;  force = true)  



    mdot_normed(uu,Zz) = (86uu^2 + 8.3*uu + 0.083)*Zz/94.383 #from http://dx.doi.org/10.1051/matecconf/201925206009
    register(myModel, :mdot_normed, 2, mdot_normed, autodiff = true)
        
    obV = get_one_by_Vsp()
    # Pijmax = maximum(nonzeros(Z) .- nonzeros(C))*1.01
    Pijmax = maximum(Zvec .- Cvec)*1.01
    Pnormed(u,Zij,Cij) = 20*(Zij*u - Cij)/Pijmax
    Λ(b,Pij) = obV(b,abs(Pij))/obV(100,0)

    register(myModel, :Λ, 2, Λ, autodiff=true)
    register(myModel, :Pnormed, 3, Pnormed, autodiff=true)
    register(myModel, :sign, 1, sign, autodiff=true)

    for timei in 2:N
        # nodej = path[timei]
        # nodei = path[timei-1]
        if linear
            @constraint(myModel, b[timei] == b[timei-1]+u[timei]*Z[timei-1]-C[timei-1])   #simple linear constraints....
            @constraint(myModel, g[timei] == g[timei-1] - u[timei]*Z[timei-1])
            # @constraint(myModel, b[timei] == b[timei-1]+u[timei]*Z[nodei,nodej]-C[nodei,nodej])   #simple linear constraints....
            # @constraint(myModel, g[timei] == g[timei-1] - u[timei]*Z[nodei,nodej])
        else
            @NLconstraint(myModel, g[timei] == g[timei-1] - u[timei]*Z[timei-1]*mdot_normed(u[timei], Z[timei-1]))
            @NLconstraint(myModel, b[timei] == b[timei-1] +  (Z[timei-1]*u[timei] - C[timei-1])*Λ(b[timei], Pnormed(u[timei], Z[timei-1], C[timei-1]) ) *sign(Pnormed(u[timei], Z[timei-1], C[timei-1])))
            # @constraint(myModel, b[timei] == b[timei-1]+u[timei]*Z[nodei,nodej]-C[nodei,nodej])   #simple linear constraints....
            # @constraint(myModel, g[timei] == g[timei-1] - u[timei]*Z[nodei,nodej])
        end
        time = (timei-1)*dt
        noise_i = get_noise_i(time, noise_intervals)
        @constraint(myModel, u[timei] <= noise_i)
    end
    
    # @constraint(myModel,[timei=2:N], u[timei] <= noiseR_along_path[timei]) #noise restrictions
    # @objective(myModel, Max, b[n]) #maximize final battery
    # @objective(myModel, Min, sum(u[t]*Z[path[t-1],path[t]] for t=2:N)) #minumize fuel use
    @objective(myModel, Max, g[end])
    
    gen_split = u_discretized(gen_MILP, N)
    set_start_value.(u,gen_split) 

    time_to_solve = @elapsed JuMP.optimize!(myModel)    
    # solution_summary(myModel)       # Provide a summary of the solution
    
    # --- Store Solution from Optimal Control Search --- #
    u_star, b_star, g_star = value.(u)[:], value.(b)[:], value.(g)[:]
    t_star = LinRange(0,tf,N) # Store off the Time Vector
    u_star[u_star .<= 0] .= 0
    sol_out = OptControlSolution(OCP.locs, OCP.C, OCP.Z, OCP.F, path, u_star, t_star, b_star, g_star, noiseR_along_path, OCP.tag, time_to_solve)
    return sol_out   
end


"""
Takes MILP instance and maps to optimal control...

Input:
N, k
prob:: string for problem loading
algo:: string for algo type
conn:: string for _4conn or ""
heur:: strig for heuristic used, "" if astar
"""
function MILP_to_opt_ctrl(N::Int64, k::Int64; prob::String = "euc_probs2D", algo::String = "_label", conn::String = "", heur::String = "")
    @load "Solutions\\$(prob)\\$(N)$(conn)_$(k)$(algo)$(heur)" pathL gen
    @load "Problems\\$(prob)\\$(N)$(conn)_$(k)" euc_inst
    C,Z = euc_inst.C, euc_inst.Z
    locs = euc_inst.locs
    GFlipped = euc_inst.GFlipped
    F = euc_inst.F
    B0, Q0, Bmax = euc_inst.B0, euc_inst.Q0, euc_inst.Bmax

    tag = get_tag(N,k, prob=prob, algo=algo, conn=conn, heur=heur)
    OCP =  OptControlProb(locs, C,Z,F, GFlipped, B0, Q0, Bmax, tag)


    return OCP, pathL, gen
end


"""
give N, k, prob, algo, conn, and heur
with_path: 1 for getting Solution\\...
with_path: 0 for just problem taf (20_4conn_1_node_euc)

"""
function get_tag(N::Int64, k::Int64; prob::String = "euc_probs2D", algo::String = "_label", conn::String = "", heur::String = "", with_path::Bool = false)
    if with_path
        return "Solutions\\OCP\\$(N)$(conn)_$(k)$(algo)$(heur)"
    else
        return "$(N)$(conn)_$(k)$(algo)$(heur)"
    end
end


# Functions needed to do perform the interpolations used lated in the verification
function interp1(X, V, Xq)
    # Interp1 is a function which provides a linear interpolation from two variables: X and V. 
    # When providing a variable lookup Xq, the provided output is the interpolated value Vq.
    knots = (X,)
    itp = interpolate(knots, V, Gridded(Linear()))
    itp[Xq]
end

function interp2(X, Y, V, Xq, Yq)
    knots = (X,Y)
    itp = interpolate(knots, V, Gridded(Linear()))
    itp[Xq, Yq]
end
function ezDynamics!(dx,x,p,t)
    #---------------------------------------------------------------------------------
    # This function develops the dynamics for the generator pattern problem (GPP)
    
    # The provided inputs are:
    # dx - derivative of the states a vector of length 2
    # x  - the states in the relative polar system a vector of length 2
    # p  - the control and the associated time vecotor in matrix form
    # t  - the current simulation time used for looking up the required 
    #---------------------------------------------------------------------------------
        # Get the interpolated control
        t_star = p[:,1];                # Pull out the time vector from the optimal control
        u_star = p[:,2];                # Pull out the control value from the optimal control
        u = interp1(t_star, u_star, t)  # interpolate the control for this current instant in time, t
        # Pull the current states
        # X_coord = x[1]
        # Y_coord = x[2]
        # Dynamics
        D_X_coord = vt*cos(u)
        D_Y_coord = vt*sin(u)
        # Write out the derivatives --- it is important to keep track thtat the 
        # derivative of a state assignment is the same in the order of the state assignment
        # that is dX[1] is the derivative of x[1] and similarly for all the states.
        dx[1] = D_X_coord
        dx[2] = D_Y_coord
    end


"""
gets noise restrictions as a vector of intervals, each labeled with 0 or 1 
0: no gen
1: gen allowed
"""
function noiseR_to_interval_sets(noise_along_path)
    intervals = Tuple{Int64, Vector{Float64}}[]
    start_time = 0
    start_val  = noise_along_path[1]
    for i in 1:length(noise_along_path)
        if noise_along_path[i] == start_val 
            continue
        else #if different, then 
            end_time = (i-1)/length(noise_along_path)
            push!(intervals, (start_val,[start_time, end_time]))
            start_time = (i-1)/length(noise_along_path)
            start_val = noise_along_path[i]
            if i == length(noise_along_path) #if last is its own chunk, add manually
                push!(intervals, (start_val,[start_time, 1]))
            end
        end
    end
    return intervals::Vector{Tuple{Int64, Vector{Float64}}}
end


"""
Give a param matrix, and the MILP path.  Give new number of time points.
Outputs vector of param values over path with new time points.
"""
function C_time_discretized(C::SparseMatrixCSC{Int64, Int64}, path::Vector{Int64}, number_timepoints_new::Int64)
    C_new = zeros(Float64, number_timepoints_new-1) #time[1] = 0....
    C_path = [C[path[idx-1],path[idx]] for idx in 2:length(path)]
    
    Δorig = 1/(length(path)-1)
    Δnew = 1/(number_timepoints_new-1)
    time_orig = 0:Δorig:1
    time_new = 0:Δnew:1
    at_prior = 0
    at_prior_orig_time = time_orig[time_orig .<= at_prior][end]
    for i in 2:(length(time_new))-1
        at_now = time_new[i]
        at_now_orig_time = time_orig[time_orig .<= at_now][end]

        if at_now_orig_time == at_prior_orig_time #then we are in same time chunk for all time, so just get weighted chunk (normalized time)
            idx_i = Int(round(at_now_orig_time*(length(path)-1))) +1 #time zero is path[1]
            idx_j = idx_i + 1
            nodei, nodej = path[idx_i], path[idx_j]
            C_new[i-1] = (Δnew/Δorig)*C[nodei, nodej]
        else #else we do a weighted average....
            t1 = at_now - at_now_orig_time
            t2 = at_now_orig_time - at_prior
            
            idx_i = Int(round(at_prior_orig_time*(length(path)-1)))+1
            idx_j = idx_i + 1
            nodei, nodej = path[idx_i], path[idx_j]
            C2bar = C[nodei, nodej]*(t2/Δorig)

            idx_i = Int(round(at_now_orig_time*(length(path)-1))) + 1
            idx_j = idx_i + 1
            nodei, nodej = path[idx_i], path[idx_j]
            C1bar = C[nodei, nodej]*(t1/Δorig)


            
            C_new[i-1] = C1bar + C2bar
        end
        at_prior = at_now
        at_prior_orig_time = time_orig[time_orig .<= at_prior][end]
    end
    # if C_new[end] == 0
    C_new[end] = (Δnew/Δorig)*C[path[end-1], path[end]] 
    # end
    return C_new
end

"""
take u from MILP and split it into a bunch of time points....
"""
function u_discretized(u_MILP::Vector{Bool}, N::Int64)
    #u[1] is at time[2] (inputting raw MILP gen)
    u_new = zeros(N) #first element should remain 0
    Δorig = 1/length(u_MILP)
    Δnew = 1/(N-1)
    Δnew > Δorig && error("need more time points..... does not work when new time steps larger than old (MILP) time steps")
    time_orig = 0:Δorig:1
    time_new = 0:Δnew:1
    at_prior = 0
    at_prior_orig_time = time_orig[time_orig .<= at_prior][end]
    for i in 2:(length(time_new))-1
        at_now = time_new[i]
        at_now_orig_time = time_orig[time_orig .<= at_now][end]
        
        if at_now_orig_time == at_prior_orig_time #then we are in same time chunk for all time, so just get weighted chunk (normalized time)
            idx_i = Int(round(at_now_orig_time*(length(u_MILP))))+1
            u_new[i] = u_MILP[idx_i]
        else #else we do a weighted average....
            t1 = at_now - at_now_orig_time
            t2 = at_now_orig_time - at_prior
            
            idx_j = Int(round(at_prior_orig_time*(length(u_MILP))))+1
            u2bar = u_MILP[idx_j]


            idx_i = Int(round(at_now_orig_time*(length(u_MILP))))+1
            u1bar = u_MILP[idx_i]
            u_new[i] = (u1bar*t1 + u2bar*t2)/(t1+t2)
            println(u_new[i])
        end
        at_prior = at_now
        at_prior_orig_time = time_orig[time_orig .<= at_prior][end]
    end
    u_new[end] = u_MILP[end]
    return u_new::Vector{Float64}
end
u_MILP = rand(Bool, 10)

N = 25

u_new = u_discretized(u_MILP, N)

function get_noise_i(timei::Float64, noise::Vector{Tuple{Int64, Vector{Float64}}})
    for noise_i in noise
        if noise_i[2][1]< timei <= noise_i[2][2]
            return Bool(noise_i[1])
        end
    end
    return 0
end