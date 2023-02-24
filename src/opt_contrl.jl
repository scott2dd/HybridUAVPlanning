
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
    path::Vector{Float64}
    gen::Vector{Float64}
    timevec::Vector{Float64}
    batt_state::Vector{Float64}
    fuel_state::Vector{Float64}
    noise_restrictions::Vector{Bool}
    tag::String
    time_to_solve::Float64
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

function MILP_to_opt_ctrl(N::Int64, k::Int64; prob::String = "euc_probs2D", algo::String = "_label", conn::String = "", heur::String = "")
    @load "Solutions\\$(prob)\\$(N)$(conn)_$(k)$(algo)$(heur)" pathL gen
    @load "Problems\\$(prob)\\$(N)$(conn)_$(k)" euc_inst
    C,Z = euc_inst.C, euc_inst.Z
    locs = euc_inst.locs
    GFlipped = euc_inst.GFlipped
    F = euc_inst.F
    B0, Q0, Bmax = euc_inst.B0, euc_inst.Q0, euc_inst.Bmax

    tag = "$(N)$(conn)_$(k)$(algo)$(heur)"
    OCP =  OptControlProb(locs, C,Z,F, GFlipped, B0, Q0, Bmax, tag)


    return OCP, pathL, gen
end
