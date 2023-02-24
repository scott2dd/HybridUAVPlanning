using JuMP
# using CPLEX


function get_A(Alist)
    N = length(Alist)

    A = zeros(Bool, N,N)

    for i in 1:N
        for j in Alist[i]
            A[i,j] = 1
            A[j,i] = 1
        end
    end

        return A
end

function MILP_hybrid(def::EucGraphInt; tlim = 900)
    S, E, Alist, F, C, Gf, Z = def.S, def.E, def.Alist, def.F, def.C, def.GFlipped, def.Z
    G = float( .! Gf)
    Bstart, Qstart = def.B0, def.Q0
    Bmax = Bstart
    SC = def.StartCost
    Cmin = minimum(nonzeros(C))
    Qstart = 999
    Bstart = Int(round(5*mean(nonzeros(C)))*9999999)
    Bmax = Bstart
    SC = 0
    anchor = def.anchor_list  #if a grid problem, this is zeros(N)
    A = get_A(Alist)    
    N = length(Alist)
    bigM = 99999
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 1, "CPXPARAM_Emphasis_Memory" => 1, "CPXPARAM_WorkMem" => 12000,  "CPXPARAM_MIP_Pool_Capacity" => 0, "CPX_PARAM_TILIM" => tlim)) 
    #, "CPX_PARAM_WORKDIR" => "C:\\Users\\Student\\CPLEX_WDIR", "CPX_PARAM_NODEFILEIND" => 1))
    @variable(m, x[i=1:N,j=1:N], Bin)
    @variable(m, g[i=1:N,j=1:N], Bin)

    @variable(m, Bmax >= b[i=1:N] >= 0)
    @variable(m, Qstart >= q[i=1:N] >= 0)

    @constraint(m, [i=1:N], x[i,i] == 0)
    @constraint(m, [i = 1:N, j = 1:N; A[i,j] == 0], x[i,j] == 0)
    @constraint(m, [i = 1:N, j = 1:N; G[i,j] == 0], g[i,j] == 0)

    @constraint(m, sum(x[S,j] for j=1:N) == 1)
    @constraint(m, sum(x[j,E] for j=1:N) == 1)

    @constraint(m, sum(x[j,S] for j=1:N) == 0)
    @constraint(m, sum(x[E,j] for j=1:N) == 0)

    @constraint(m, [i = setdiff(1:N, [S,E])], sum(x[i,j] for j=1:N)-sum(x[j,i] for j=1:N) == 0)


    @constraint(m, b[S] == Bstart)
    # @constraint(m, [i = 1:N, j = 1:N],  b[j] <=  b[i] - C[i,j] + C[i,j]*F[i,j] + Z[i,j]*g[i,j]
                                                #  - SC*(1 - sum(g[k,i] for k = setdiff(1:N, j))) + bigM*(1-x[i,j]))

    @constraint(m, [i = 1:N, j = 1:N], b[j] <=  b[i] - C[i,j] + C[i,j]*F[i,j] + Z[i,j]*g[i,j] + bigM*(1-x[i,j]))

    @constraint(m, q[S] == Qstart)
    @constraint(m, [i = 1:N, j = 1:N], q[j] <=  q[i] - Z[i,j]*g[i,j] + bigM*(1-x[i,j]))

    @constraint(m, [i=1:N, j = 1:N], g[i,j] <= x[i,j])
    @constraint(m, [i=1:N, j = 1:N], g[i,j] <= G[i,j])


    # @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N) + sum(g[i,j]*(0.1*C[i,j]) for  i=1:N, j=1:N))
    @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N) )

    #add constraints such that xij = 0 for all edges that don't exist
    for i = 1:N
        for j = 1:N
            if A[i,j] == 0
                @constraint(m, x[i,j] == 0)
                @constraint(m, g[i,j] == 0)
            end

        end
    end


    
    optimize!(m)
    xsol = value.(x)
    gsol = value.(g)
    bstate = value.(b)
    qstate = value.(q)

    time = solve_time(m)
    cost = objective_value(m)
    # println(xsol)
    path, gen = IP_to_vec(xsol, gsol, S, E)
    return time, cost, path, gen
end


# @load "Problems\\euc_probs2D\\50_2" euc_inst
# tdp, cost, pathL, gen = MILP_hybrid(euc_inst)











# ## Load Prob
# using JLD2
# println("MILP\n")
# xx, yy, zz = 5,5,5
# k = 5
# @load "Problems\\grid_probs3D\\$(xx)x$(yy)x$(zz)_$(k)" prob
# @load "Problems\\grid_probs3D\\GraphDef_$(xx)x$(yy)x$(zz)" graphdef

# full_def = FullDef3D(prob.S, prob.E, graphdef.Alist, graphdef.A, graphdef.F, graphdef.C, .!prob.GFlipped, graphdef.Z, prob.B0, 500, prob.Bmax, prob.StartCost, prob.anchor_list, prob.Dim)

# # CPLEX Test
# outMILP = []
# try
#     outMILP = MILP_hybrid(full_def, tlim = 15*60)
#     println("Integer Solutions found...")
# catch e
#     println("0 Solutions in model")
# end
