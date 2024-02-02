function get_path(label::Label, came_from::Vector{Vector{Tuple{Int64,Int64}}}, start::Int64)
    path = Int64[]
    here = label.node_idx
    cf_idx_here = label.came_from_idx

    here == start && return [start]
    push!(path, here)
    while here != start
        next = came_from[here][cf_idx_here][1]
        cf_idx_next = came_from[here][cf_idx_here][2]

        push!(path, next)
        here = copy(next)
        cf_idx_here = copy(cf_idx_next)
    end
    return reverse(path)
end
function get_gen(label::Label, gen_track::Vector{Vector{Tuple{Int64,Int64}}})
    #get generator pattern from recursive data struct
    genOut = Bool[]
    PL = label.pathlength #path length of tracked label (in # of nodes) ... so if PL=1 then no edges, 
    gt_idx = label.gentrack_idx #index for gen_track
    while PL > 1
        gen_now = gen_track[PL][gt_idx][1]
        push!(genOut, gen_now)
        gt_idx = gen_track[PL][gt_idx][2]
        PL -= 1
    end
    return reverse(genOut)
end


function update_path_and_gen!(new_label::L, came_from::Vector{Vector{Tuple{Int64,Int64}}}, gen_track::Vector{Vector{Tuple{Int64,Int64}}}) where L<:Label
    #correct path...
    pnode = new_label.prior_node_idx
    nextnode = new_label.node_idx
    p_came_from_idx = new_label._hold_came_from_prior
    path_pointer = findall(x -> x == [pnode, p_came_from_idx], came_from[nextnode])
    if isempty(path_pointer) #if no other label has used this same path...
        push!(came_from[nextnode], (pnode, p_came_from_idx))
        came_from_idx = length(came_from[nextnode])
    else #if path exists prior, then we use the (nonempty) pointer
        pointer_idx = path_pointer[1]
        came_from_idx = pointer_idx #label now has index for came_from 
    end


    #correct gen....
    PL = new_label.pathlength
    gen_pointer = findall(x -> x == [new_label.gen_bool, new_label._hold_gen_track_prior], gen_track[new_label.pathlength])
    if isempty(gen_pointer)
        push!(gen_track[PL], (new_label.gen_bool, new_label._hold_gen_track_prior))
        gentrack_idx = length(gen_track[PL])
    else
        pointer_idx = gen_pointer[1]
        gentrack_idx = pointer_idx #label now has index for gen_track
    end

    label_updated = L(
        gcost=new_label.gcost,
        fcost=new_label.fcost,
        hcost=new_label.hcost,
        node_idx=new_label.node_idx,
        prior_node_idx=new_label.prior_node_idx,
        _hold_came_from_prior=new_label._hold_came_from_prior,
        came_from_idx=came_from_idx,
        pathlength=new_label.pathlength,
        _hold_gen_track_prior=new_label._hold_gen_track_prior,
        gentrack_idx=gentrack_idx,
        gen_bool=new_label.gen_bool,
        batt_state=new_label.batt_state,
        gen_state=new_label.gen_state,
    )
    return label_updated
end


function get_path(label::Vector{Int64}, came_from::Vector{Vector{Vector{Int64}}}, start::Int64)
    #get path from recursive data struct
    path = Int64[]
    here = label[4]
    cf_idx_here = label[6]
    
    here == start && return [start]
    push!(path, here)
    while here != start
        next = came_from[here][cf_idx_here][1]
        cf_idx_next = came_from[here][cf_idx_here][2]
        push!(path, next)

        here = copy(next)
        cf_idx_here = copy(cf_idx_next)
    end
    return reverse(path)
end
function get_gen(label::Vector{Int64}, gen_track::Vector{Vector{Vector{Int64}}}) 
    #get generator pattern from recursive data struct
    genOut = Bool[]
    PL = label[7] #path length of tracked label (in # of nodes) ... so if PL=1 then no edges, 
    gt_idx = label[8]
    while PL > 1
        gen_now = gen_track[PL][gt_idx][1]
        push!(genOut, gen_now)
        gt_idx = gen_track[PL][gt_idx][2]
        
        PL -= 1
    end
    return reverse(genOut)
end


function get_heur_astar(E, locs, Fvec::Vector{Float64})
    function heur_astar(i)
        if isnan(Fvec[i])
            f = norm(locs[i,:] - locs[E,:])
        else
            f = Fvec[i]
        end
        return f
    end
    return heur_astar 
end
function get_heur_astar(E, locs, Fvec::Vector{Int64})
    function heur_astar(i)
        if isnan(Fvec[i]) || Fvec[i] == 2^63 - 1
            f = Int(floor(norm(locs[i,:] - locs[E,:])))
        else
            f = Fvec[i]
        end
        return f
    end
    return heur_astar 
end


function get_heur_label(Fvec::Vector{Float64}, graph, C, E, heur_astar)
    
    function heur_label!(i)
        #if needed, run astar.  then, add all nodes in path to Fvec. return only f cost for this label, while altering all of Fvec
        if isnan(Fvec[i])
            astar_out = a_star(graph, i, E, C, heur_astar) 
            if isempty(astar_out)
                path_list = [i] 
                printstyled("No path to goal!!", color=:light_red)
                LB_vec = [Inf]
            else
                path_list, LB_vec, cost = astar_proc(astar_out, C)
            end
            add_to_Fvec!(LB_vec, path_list, Fvec)
        end
        return Fvec[i]
    end
    
    return heur_label!
end
function get_heur_label(Fvec::Vector{Int64}, graph, C::SparseMatrixCSC{Int64, Int64}, E, heur_astar) #def if discrete version
    function heur_label!(i)
        #if needed, run astar.  then, add all nodes in path to Fvec. return only f cost for this label, while altering all of Fvec
        if Fvec[i] == 	2^63 - 1
            astar_out = a_star(graph, i, E, C, heur_astar)
            if isempty(astar_out)
                path_list = [i] 
                printstyled("No path to goal!!", color=:light_red)
                LB_vec = [	2^63 - 1]
            else
                path_list, LB_vec, cost = astar_proc(astar_out, C)
            end
            LB_vec = Int.(LB_vec)
            add_to_Fvec!(LB_vec, path_list, Fvec)
        end
        return Fvec[i]
    end
    
    return heur_label!
end



function get_heur_label_euc(Fvec::Vector{Float64}, locs, E)
    function heur_label!(i)
        if isnan(Fvec[i])
            heur = norm(locs[i,:] - locs[E,:]) 
            Fvec[i] = heur
        end
        return Fvec[i]
    end
    return heur_label!
end
function get_heur_label_euc(Fvec::Vector{Int64}, locs, E) #def if we're using euc distance
    function heur_label!(i)
        if Fvec[i] == 2^63 - 1
            heur = Int(floor(norm(locs[i,:] - locs[E,:]))) 
            Fvec[i] = heur
        end
        return Fvec[i]
    end
    return heur_label!
end



function add_to_Fvec!(LB_vec, path, Fvec)
    for idx in 1:length(path)
        i = path[idx]
        Fvec[i] = LB_vec[idx]
    end
end



function merge_3D(X::Vector{Vector{Int64}}, Y::Vector{Vector{Int64}}) #O(N^2 / 2 + N / 2 + NlogN)  -> slightly better than O(N^2)
    Mtemp::Vector{Vector{Int64}} = vcat(X, Y)
    for i in 1:length(Mtemp)
        Mtemp[i][2:3] .*= -1
    end
    sort!(Mtemp) #sorts by 1st el then 2nd el then 3rd el....
    boolv = ones(Bool, length(Mtemp))
    #now go through each element 
    for i in eachindex(Mtemp)
        label =  Mtemp[i]
        for ii in 1:(i-1)
            boolv[ii] == 0 && continue #if already dommed don't need to compare
            if dom_min(Mtemp[ii], label)
                boolv[i] = 0
                break
            end
     end
    end
    M_out = Mtemp[boolv] 
    for i in 1:length(M_out)
        M_out[i][2:3] .*= -1
    end
    return M_out
end

function dom_min(X::Vector{Int64},Y::Vector{Int64})
    #returns true if X >> Y
    bool = false
    (X[1] <= Y[1] && X[2] <= Y[2]  && X[3] <= Y[3]) && (bool = true)  
    return bool 
end
function dom(X::Vector{Int64},Y::Vector{Int64})
    #returns true if X >> Y
    bool = false
    (X[1] <= Y[1] && X[2] >= Y[2]  && X[3] >= Y[3]) && (bool = true)
    return bool
end
function EFF_list(Γ::Vector{Vector{Int64}}, new::Vector{Int64})  #should not need this, but may be faster if we have less? like a merge sort?
    EFF_bool = true
    for i in 1:length(Γ)  #γ in Γ
        γ = Γ[i]
        if dom(γ, new)
            EFF_bool = false
            return EFF_bool
        end
    end
    return EFF_bool
end


function EFF_heap(Q::MutableBinaryMinHeap{L}, label_new::L) where {L<:Label}
    isempty(Q) && (return true)
    node_map_copy = Q.node_map
    for k in 1:length(node_map_copy)
        node_map_copy[k] == 0 && continue
        Q[k].node_idx != label_new.node_idx && continue #if they are different nodes, then skip...

        (Q[k].gcost <= label_new.gcost && Q[k].batt_state >= label_new.batt_state && Q[k].gen_state >= label_new.gen_state) && (return false)
    end
    return true #if all this passes, then return true (is efficient)
end

#P should be vector of vector of labels
function EFF_P(P::Vector{Vector{L}}, label_new::L) where {L<:Label}
    #loop through P_i and return 0 if dominated or
    i = label_new.node_idx
    isempty(P[i]) && (return true)
    for label_closed in P[i]
        (label_new.gcost >= label_closed.gcost && label_new.batt_state <= label_closed.batt_state && label_new.gen_state <= label_closed.gen_state) && (return false) #then return false
    end
    return true #otherwise, return true.... Can we use a dictionary???????
end


function EFF_heap(Q::MutableBinaryMinHeap, label_new::Vector{T}) where T<:Number 
    isempty(Q) && (return true)
    node_map_copy = Q.node_map
    for k in 1:length(node_map_copy)
        node_map_copy[k] == 0 && continue
        Q[k][2][4] != label_new[4] && continue #if they are different nodes, then skip...
        #if label_new is dominated... return 0 and break
        (Q[k][2][1] <= label_new[1] && Q[k][2][2] >= label_new[2] && Q[k][2][3] >= label_new[3]) && (return false)
        #if label_new is equivalent... return 0 and break
        (Q[k][2][1] == label_new[1] && Q[k][2][2] == label_new[2] && Q[k][2][3] == label_new[3]) && (return false)
    end
    return true
end
function EFF_P(P::Vector{Matrix{T}}, label_new::Vector{T}) where T<:Number
    #loop through P_i and return 0 if dominated or
    i = Int(label_new[4])
    isempty(P[i]) && (return true)
    for k in 1:size(P[i],1)
        (P[i][k,1] <= label_new[1] && P[i][k,2] >= label_new[2] && P[i][k,3] >= label_new[3]) && return false
        (P[i][k,1] == label_new[1] && P[i][k,2] == label_new[2] && P[i][k,3] == label_new[3]) && return false
    end
    return true
end







function load_full_def(Dim, k)
    xx, yy, zz = Dim[1], Dim[2], Dim[3]
    @load "Problems\\grid_probs3D\\$(xx)x$(yy)x$(zz)_$(k)" prob
    @load "Problems\\grid_probs3D\\GraphDef_$(xx)x$(yy)x$(zz)" graphdef

    full_def = FullDef3D(prob.S, prob.E, graphdef.Alist, graphdef.A, graphdef.F, graphdef.C, .!prob.GFlipped, graphdef.Z, prob.B0, prob.Q0 + 500, prob.Bmax, prob.StartCost, prob.anchor_list, prob.Dim)
    return full_def
end



function IP_to_vec(xsol, gsol, S, E)
    path = Vector{Int}(undef, 0)
    gen = Vector{Bool}(undef, 0) #hold gen info
    k = S
    prev = S
    push!(path, S)
    while true
        next = findall(x->x>.98,xsol[k,:])[1]
        genbool = gsol[prev, next]
        push!(path, next)
        push!(gen, genbool)
        k = next
        next == E && return path, gen
    end
end
# 3D Utils Functions
# xx = num rows
# yy = num cols
# zz = num slives (height)
function get_single_idx3D(idx, Dim)
    xx, yy, zz = Dim[1], Dim[2], Dim[3]

    i,j,k = idx[1], idx[2], idx[3]
    out = xx*yy*(k-1) + yy*(i-1) + j # CHECK I J K order....
    return Int(out)
end


function get_3D_idx(node, Dim)
    xx, yy, zz = Dim[1], Dim[2], Dim[3]
    k = Int(floor(node/(xx*yy+.0001)) + 1)
    i = Int( floor((node - (k-1)*xx*yy)/(yy+.0001))+1 ) #row
    j = node - (k-1)*xx*yy - (i-1)*yy  #col


    return [i,j,k]
end
function manhattan(j,E, Dim, anchor)
    j_ij = get_3D_idx(j, Dim)
    E_ij = get_3D_idx(E, Dim)
    dist = sum(abs.(E_ij - j_ij))
end
function manhattan(j,E, Dim)
    j_ij = get_3D_idx(j, Dim)
    E_ij = get_3D_idx(E, Dim)
    dist = sum(abs.(E_ij - j_ij))


end
function euc_dist_plain3D(j,E, Dim)
    j_ij = get_3D_idx(j, Dim)
    E_ij = get_3D_idx(E, Dim)
    dist = norm(E_ij - j_ij)
    # dist = sum(abs.(E_ij - j_ij))
    return dist
end

function euc_dist_eucgraph(j,E, locs)
    dist = norm(locs[j,:] - locs[E, :])
    return dist
end

function euc_dist3D(j, E, Dim, anchor)
    if E <= prod(Dim) #n*m
        return euc_dist_plain3D(j,E, Dim)
    elseif E > prod(Dim)
        return euc_dist_plain3D(j,anchor,Dim) + 1
    end
end

#2D Utility Functions
function path_cost(path, C)
    cost = 0
    for idx in 2:length(path)
        i,j = path[idx-1], path[idx]
        Cij = C[i,j]
        cost += Cij
    end
    return cost
end
function get_single_idx(idx, Dim)
    n = Dim[1] #number of rows
    m = Dim[2] #lengthof rows

    row = idx[1]
    col = idx[2]

    out = (row-1)*m
    out += col
    return Int(out)
end


function get_2d_idx(node, Dim)
    nn = Dim[1]
    row = Int(floor(node/(nn+.0001)) + 1)
    col = node - nn*(row - 1)
    return [row, col]
end

function euc_dist_plain(j,E, Dim)
    j_ij = get_2d_idx(j, Dim)
    E_ij = get_2d_idx(E, Dim)
    dist = norm(E_ij - j_ij)
end


function euc_dist(j, E, Dim, anchor)
    if E <= Dim[1]*Dim[2] #n*m
        return euc_dist_plain(j,E, Dim)
    elseif E > Dim[1]*Dim[2]
        return euc_dist_plain(j,anchor,Dim) + 1
    end
end

function make_graph(euc_inst::EucGraphInt)
    A = euc_inst.Alist
    N = length(A)
    g = SimpleWeightedDiGraph(N)
    weights = euc_inst.C
    for i in 1:N
        for j in A[i]
            add_edge!(g, i, j, weights[i,j])
        end
    end
    return g
end

function make_graph(euc_inst::EucGraphInt, weighted::Bool)
    A = euc_inst.Alist
    N = length(A)
    g = SimpleGraph(N)
    for i in 1:N
        for j in A[i]
            add_edge!(g, i, j)
        end
    end
    return g
end

function astar_proc(path, C)
    LB_vec = zeros(length(path) + 1)
    path_list = zeros(Int, length(path) + 1)

    for i in 1:length(path)
        e = path[i]
        j,k = e.src, e.dst
        path_list[i] = j
    end
    path_list[end] = path[end].dst #add dest to end for last node
    #now get total path cost
    cost = 0
    for e in path
        i,j = e.src, e.dst
        cost += C[i,j]
    end


    idx = length(path_list)
    LB_vec[end] = 0
    for e in reverse(path)
        i,j = e.dst, e.src
        LB_vec[idx-1] = LB_vec[idx] + C[i,j]
        idx -= 1 
    end
    return path_list, LB_vec, cost
end

function hybrid_with_LB_depth(def::EucGraph)
    S, E, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
    Bstart, Qstart = def.B0, def.Q0
    @time Cmin = minimum(nonzeros(C))
    # Bstart = mean(nonzeros(C)) * 4
    genpen = Cmin*0
    Bstart = 9999+mean(nonzeros(C)) + 5*maximum(C)
    Qstart = 9999
    Bmax = Bstart
    SC = def.StartCost
    SC = 0
    anchor = def.anchor_list  #if a grid problem, this is zeros(N)
    locs = def.locs
    # G = .!def.GFlipped
    Dim = def.locs           ####### pass third argument of heuristics as locs
   
    N = length(Alist)
    # println("Solving with $(N) Nodes || B0 = $(Bstart)")
    graph = make_graph(def)
    Fvec = fill(NaN, N)
    heur_astar = get_heur_astar(E, locs, Fvec)
    heur_astar(S)
    heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
    if heur_label!(S) == Inf     
        println("Infinite Start Heuristic...")
        return 0,0,0,0,0,0,0
    end
    
    Q = MutableBinaryMinHeap([   ((heur_label!(S),  ), [0.0, Bstart, Qstart, S, 1.0, heur_label!(S) + 0]) ] )
    P = [zeros(0,6) for k = 1:size(C,1)] #set of treated labels | new labels compare to this AND Q
    X = Vector{Int}[ [S] ]  #hold path info
    Y = Vector{Int}[  []   ] #hold gen info
    
    z=0
    look_stack = Int[]
    dist_stack = Float32[]
    while true #loop until get to end node, or Q is empty
        isempty(Q) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break)

        #pull minimum cost label....
        next = pop!(Q)
        label_treated = next[2]
        i = Int(label_treated[4])
        K = Int(label_treated[5])
        GenPrior = 1 #label_treated[7]
        append!(look_stack, i)
        append!(dist_stack, heur_label!(i))
        #now add it to P
        P[i] = vcat(P[i], reshape(label_treated,(1,6)))
        # println(i)
        #if we are at the end, then stop algo and return this label
        if i == E
            opt_cost =  label_treated[1]
            opt_path =  X[K]
            opt_gen =   Y[K]
            return opt_cost, opt_path, opt_gen, length(X), look_stack, dist_stack
        end

         for j in Alist[i]
            j==i && continue
            j∈X[K] && continue
            gc = label_treated[1] + C[i,j]
            h =  heur_label!(j) 
            Fbin = F[i,j] # if we can glide, then ALWAYS glide
            label = []
            # No Gen
            if label_treated[2] + Z[i,j] < Bmax && def.GFlipped[i,j] == 0
                        # Gen On   IF allowed...
                if def.GFlipped[i,j] == 0 && label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior) >= 0 && label_treated[3]-Z[i,j] >= 0 && label_treated[2]-C[i,j] + Z[i,j] <= Bmax
                    label =  [gc + genpen, label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior), label_treated[3]-Z[i,j],  j, Int(length(X)+1), gc+h]
                    # println("label, $(i)->$(j) GEN")

                    if EFF_heap(Q, label) && EFF_P(P, label)
                        push!(Q, (gc+h + genpen,label))

                        push!(X, [X[K]; j ])
                        push!(Y, [Y[K]; 1 ])
                    end
                end
            else
                if label_treated[2]-C[i,j]*(1-Fbin) >= 0
                    label =  [gc, label_treated[2]-C[i,j]*(1-Fbin), label_treated[3],  j, Int(length(X)+1), gc+h]
                    # println("label, $(i)->$(j)")
                    
                    if EFF_heap(Q, label) && EFF_P(P, label)
                        push!(Q, (gc+h,label))
                        push!(X, [X[K]; j ])
                        push!(Y, [Y[K]; 0 ])
                    end
                end
            end
        end
        z+=1
        z== 50000 && (printstyled("Z BREAK... \n", color=:light_red); break)
    end
    #let's just return the "best" path we get so far...
    next = pop!(Q)
    label_treated = next[2]
    i = Int(label_treated[4]); append!(look_stack, i)
    K = Int(label_treated[5])
    out_cost =  label_treated[1]
    out_path =  X[K]
    out_gen =   Y[K]
    return out_cost, out_path, out_gen, length(X), look_stack, dist_stack
    
end
