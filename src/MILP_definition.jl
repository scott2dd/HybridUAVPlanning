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

function MILP_edge_notation(def::FullDef3D; tlim = 900)
    S, E, A, Alist, F, C, G, Z = def.S, def.E, def.A, def.Alist, def.F, def.C, def.G, def.Z
    Bstart, Qstart = def.B0, def.Q0
    Bmax = Bstart
    SC = def.StartCost
    anchor = def.anchor_list  #if a grid problem, this is zeros(N)
    Dim = def.Dim    # Dim = [Int(sqrt(size(G,1))), Int(sqrt(size(G,1)))]
    
    N = Int(prod(Dim)) #number of nodes 
    edict = get_edge_dict(Alist) #number of edges...
    nE = length(elist)
    Glist = get_G_list(G, edict)
    
    
    bigM = 99
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 1, "CPXPARAM_Emphasis_Memory" => 1, "CPXPARAM_WorkMem" => 12000,  "CPXPARAM_MIP_Pool_Capacity" => 0, "CPX_PARAM_TILIM" => tlim, "CPX_PARAM_WORKDIR" => "C:\\Users\\Student\\CPLEX_WDIR", "CPX_PARAM_NODEFILEIND" => 1))
    
    @variable(m, e[i=1:nE], Bin)
    @variable(m, g[i=1:nE], Bin)
    @variable(m, u[i=1:N])
    @variable(m, Bmax >= b[i=1:N] >= 0)
    @variable(m, Qstart >= q[i=1:N] >= 0)

    # @constraint(m, [i=1:N], x[i,i] == 0) #don't need this, no edges from i->i exist to begin w/
    # @constraint(m, [i = 1:N, j = 1:N; A[i,j] == 0], x[i,j] == 0)  #nodes without adjacency cant have an edge
    @constraint(m, [i = 1:N, j = 1:N; G[i,j] == 0], g[i,j] == 0)

    
    # 2 constriants orig: must leave start, must arrive at end////
    # @constraint(m, sum(x[S,j] for j=1:N) == 1)
    # @constraint(m, sum(x[j,E] for j=1:N) == 1)

    @constraint(m, sum(e[j] for j in get_edges_on_node(edict, S)) == 1) #exactly one edge connected to S
    @constraint(m, sum(e[j] for j in get_edges_on_node(edict, E)) == 1) #exactly one edge connected to E


    #ORIG Constraints: cant go into start, cant leave end -> we do not need these, as there are no edges in/out, only edges
    # @constraint(m, sum(x[j,S] for j=1:N) == 0)
    # @constraint(m, sum(x[E,j] for j=1:N) == 0)

    #if edge in, edge must be out
    #for every node: 
        #for each edge, edge - sum(all other edges on this node) = 0
    for i = 1:N
        edges = get_edges_on_node(edict, i)
        for edge in edges
            other_edges = setdiff(edges, edge)
            println(other_edges)
            @constraint(m, sum(e[j] for j ∈ other_edges) - e[edge] >= 0)
        end
    end
    #orig: 
    #@constraint(m, [i = setdiff(1:N, [S,E])], sum(x[i,j] for j=1:N)-sum(x[j,i] for j=1:N) == 0)
    


    @constraint(m, b[S] == Bstart)
    @constraint(m, [i = 1:N, j = 1:N],  b[j] <=  b[i] - C[i,j] + C[i,j]*F[i,j] + Z[i,j]*g[i,j] - SC*(1 - sum(g[k,i] for k = 1:N)) + bigM*(1-x[i,j]))

    # @constraint(m, [i = 1:N, j = 1:N], b[j] <=  b[i] - C[i,j] + C[i,j]*F[i,j] + Z[i,j]*g[i,j] + bigM*(1-x[i,j]))

    @constraint(m, q[S] == Qstart)
    @constraint(m, [i = 1:N, j = 1:N], q[j] <=  q[i] - Z[i,j]*g[i,j] + bigM*(1-x[i,j]))

    @constraint(m, [i=1:N, j = 1:N], g[i,j] <= x[i,j])
    @constraint(m, [i=1:N, j = 1:N], g[i,j] <= G[i,j])


    # @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N) + sum(g[i,j]*(0.1*C[i,j]) for  i=1:N, j=1:N))
    @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N) )

    optimize!(m)
    xsol = value.(x)
    gsol = value.(g)
    bstate = value.(b)
    qstate = value.(q)

    time = solve_time(m)
    cost = objective_value(m)
    return cost, xsol, gsol, bstate, qstate, time
end
 
 
 
 function get_edge_dict(Alist)
    list = Set()
    count = 0
    for i in 1:length(Alist)
        L = Alist[i]
        for Ll in L
            a = sort([Ll, i])
            count += 1
            push!(list, a)
        end
    end
    dict = Dict()
    count = 1
    for el in list
        dict[count] = el
        count += 1
    end
    return dict
 end

 function get_G_list(G, edict)
    Gvec = zeros(length(edict))
    for (nE, edge) in edict
        Gvec[nE] = G[CartesianIndex(Tuple(edge))]
        
    end

 end

 function get_edges_on_node(edict, I)
    edgelist = []
    for (e, edge) in edict
        if I ∈ edge
            push!(edgelist, e)
        end
    end
    return edgelist
 end



 function lower_boundLP(def::EucGraph; tlim = 900)
    S, E, Alist, F, C, GFlipped, Z = def.S, def.E, def.Alist, def.F, def.C, def.GFlipped, def.Z
    Bstart, Qstart = def.B0, def.Q0
    Bmax = Bstart
    SC = def.StartCost
    anchor = def.anchor_list  #if a grid problem, this is zeros(N)
 
    G = .!GFlipped

    N = length(Alist)
    bigM = 99
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPXPARAM_Emphasis_Memory" => 1, "CPXPARAM_WorkMem" => 15000,  "CPXPARAM_MIP_Pool_Capacity" => 0, "CPX_PARAM_TILIM" => tlim, "CPX_PARAM_WORKDIR" => "C:\\Users\\Student\\CPLEX_WDIR", "CPX_PARAM_NODEFILEIND" => 1))
    @variable(m, x[i=1:N,j=1:N]) #  , Bin)
    @variable(m, g[i=1:N,j=1:N]) # , Bin)

    @variable(m, Bmax >= b[i=1:N] >= 0)
    @variable(m, Qstart >= q[i=1:N] >= 0)

    @constraint(m, [i=1:N], x[i,i] == 0)
    # @constraint(m, [i = 1:N, j = 1:N; A[i,j] == 0], x[i,j] == 0)
    @constraint(m, [i = 1:N, j = 1:N; G[i,j] == 0], g[i,j] == 0)

    @constraint(m, sum(x[S,j] for j=1:N) == 1)
    @constraint(m, sum(x[j,E] for j=1:N) == 1)

    @constraint(m, sum(x[j,S] for j=1:N) == 0)
    @constraint(m, sum(x[E,j] for j=1:N) == 0)

    @constraint(m, [i = setdiff(1:N, [S,E])], sum(x[i,j] for j=1:N)-sum(x[j,i] for j=1:N) == 0)


    @constraint(m, b[S] == Bstart)
    @constraint(m, [i = 1:N, j = 1:N],  b[j] <=  b[i] - C[i,j] + C[i,j]*F[i,j] + Z[i,j]*g[i,j]
                                                - SC*(1 - sum(g[k,i] for k = setdiff(1:N, j))) + bigM*(1-x[i,j]))

    # @constraint(m, [i = 1:N, j = 1:N], b[j] <=  b[i] - C[i,j] + C[i,j]*F[i,j] + Z[i,j]*g[i,j] + bigM*(1-x[i,j]))

    @constraint(m, q[S] == Qstart)
    @constraint(m, [i = 1:N, j = 1:N], q[j] <=  q[i] - Z[i,j]*g[i,j] + bigM*(1-x[i,j]))

    @constraint(m, [i=1:N, j = 1:N], g[i,j] <= x[i,j])
    @constraint(m, [i=1:N, j = 1:N], g[i,j] <= G[i,j])


    # @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N) + sum(g[i,j]*(0.1*C[i,j]) for  i=1:N, j=1:N))
    @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N) )

    #add constraints such that xij = 0 for all edges that don't exist
    for i = 1:N
        for j = 1:N
            if j ∈ Alist[i]
                @constraint(m, x[i,j] == 0)
                @constraint(m, g[i,j] == 0)
            end

        end
    end


    optimize!(m)
    #  xsol = value.(x)
    #  gsol = value.(g)
    #  bstate = value.(b)
    #  qstate = value.(q)

    time = solve_time(m)
    cost = objective_value(m)
    # path, gen = IP_to_vec(xsol, gsol, S, E)
    return time, cost#, path, gen
 end

 

function feas_check(euc_inst)
    graph =make_graph(euc_inst)
    ww
    #first check if theres a path from start to finish...
    astar_out = a_star(graph, euc_inst.S, euc_inst.E, euc_inst.C)
    # println(astar_out)
    if isempty(astar_out)
        return 0
    end
    #now do LP check...
    # G = euc_inst.G
    # S = euc_inst.S
    # E = euc_inst.E
    # Alist = euc_inst.Alist
    # #now check if theres a path _without_ the noise restricted zones
    # if sum(G[E,:] == 0) && sum(G[E,:] == 0) #if E and S are NOT noise restricted
    #     #for every noise restricted node, remove it from the new graph
    #     for i in 1:length(Alist)
    #         if !sum(G[E,:] > 0) #if the node is not&& X[3] <= Y[3] nois restricted then continue
    #             continue
    #         end
    #         #if here, then i is noise restricted
    #         for j in Alist  #for every node connected to i, remove edges (i,j) and (j,i)
    #             deleteat!(Alist[j], Alist[j] .== i)
    #             deleteat!(Alist[i], Alist[i] .== j)
    #         end
    #     end
    #     #now have altered Alist.... get Alist 
    # end
    # #if E node is in noise restricted... need another plan...
    return 1
end






