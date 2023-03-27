#this comment was added on laptop.  Test to see if git pushing between machines is working!
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
function solve_gen_optimal_control(OCP::OptControlProb, path::Vector{Int64}, gen::Vector{Bool}; linear = true)
    myModel = Model(Ipopt.Optimizer)
    set_silent(myModel)
    
    N = length(path) # make u[1] = 0... but need length(path) for indexing wrt path
    tf = 1 #assume path takes 1 second
    dt = (tf)/(N-1)
    
    #pull from OCP
    C,Z,GFlipped = OCP.C, OCP.Z, OCP.GFlipped

    g_min,g_max = 0, OCP.Q0     # Bounds on the States
    b_min, b_max = 6, 99.9
    
    u_min, u_max = 0, 1     # Bounds on Control

    noiseR_along_path = Float64.([!OCP.GFlipped[path[idx-1],path[idx]] for idx in 2:length(path)])
    pushfirst!(noiseR_along_path, 0)

    @JuMP.variables(myModel, begin
        # dt ≥ 0, (start = dt_guess) # Delta T for time steps
        g_min  ≤ g[1:N]  ≤ g_max    # fuel level, fuel is monotomically decreasing
        b_min  ≤ b[1:N] ≤ b_max      # battery state of charge
        u_min  ≤ u[1:N] ≤ u_max      # generator setting
    end)

    fix(b[1], OCP.B0; force = true)  
    fix(g[1], OCP.Q0; force = true)  
    fix(u[1], 0;  force = true)  



    mdot_normed(uu,Zz) = (86uu^2 + 8.3*uu + 0.083)*Zz/94.383 #from http://dx.doi.org/10.1051/matecconf/201925206009
    register(myModel, :mdot_normed, 2, mdot_normed, autodiff = true)
        
    obV = get_one_by_Vsp()
    Pijmax = maximum(nonzeros(Z) .- nonzeros(C))
    Pnormed(u,Zij,Cij) = 20*(Zij*u - Cij)/Pijmax
    Λ(b,Pij) = obV(b,abs(Pij))/obV(100,0)

    register(myModel, :Λ, 2, Λ, autodiff=true)
    register(myModel, :Pnormed, 3, Pnormed, autodiff=true)
    register(myModel, :sign, 1, sign, autodiff=true)

    for timei in 2:N
        nodej = path[timei]
        nodei = path[timei-1]
        
        # println(noiseR_along_path[timei])

        if linear
            @constraint(myModel, b[timei] == b[timei-1]+u[timei]*Z[nodei,nodej]-C[nodei,nodej])   #simple linear constraints....
            @constraint(myModel, g[timei] == g[timei-1] - u[timei]*Z[nodei,nodej])
        else
            @NLconstraint(myModel, g[timei] == g[timei-1] - u[timei]*Z[nodei,nodej]*mdot_normed(u[timei], Z[nodei,nodej]))
            @NLconstraint(myModel, b[timei] == b[timei-1] +  (Z[nodei,nodej]*u[timei] - C[nodei, nodej])*Λ(b[timei], Pnormed(u[timei], Z[nodei,nodej], C[nodei,nodej]) ) *sign(Pnormed(u[timei], Z[nodei,nodej], C[nodei,nodej])))
        end
        
    end
    @constraint(myModel,[i=2:N], u[timei] <= noiseR_along_path[timei]) #noise restrictions
    # @objective(myModel, Max, b[n]) #maximize final battery
    Zvec = [Z[path[i-1],path[i]]  for i = 2:N]
    @objective(myModel, Min, sum(u[t]*Z[path[t-1],path[t]] for t=2:N)) #minumize fuel use

    set_start_value.(u,[0;gen]) 

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


