module TestPkg


import LinearAlgebra.norm
using DataStructures
using Graphs
using JLD2
using Compose
using GLM
using ProgressMeter
using DataFrames
using SparseArrays
using DataStructures
using JLD2
using Gadfly
using Cairo
using Compose
using GraphPlot
import Statistics.mean
using Colors 
using Interpolations
using JuMP
using CPLEX
using Ipopt 
using DifferentialEquations
using Interpolations        





include("struct_defs.jl")
include("label_label_sel.jl")
include("label_node_sel.jl")


export hybrid_label_selection
export hybrid_node_selection
export EucGraphInt

export get_path
export get_gen


# export get_heur_astar
# export manhattan
# export get_heur_label

# export euc_dist
# export euc_dist3D
# export euc_dist_eucgraph
# export get_heur_label
# export get_heur_label_euc
# export get_heur_label_euc_disc

export plot_euc_graph
export plot_euc_graph_labeled

export get_sol_vec
export get_pathlength_vec
export get_avg_layer
export get_OCV_func
export get_one_by_Vsp
export get_OCV_func
export get_one_by_VspGLM
export get_OCV_func_LiPo
export get_OCV_table
export get_one_by_V
export riemman
export rieman2
export V_model2
export riemman7
export simps
#add comment


# export  merge_3D
# export EFF_list
# export EFF_heap 
# export EFF_P

export OptControlProb
export OptControlSolution
export MILP_to_opt_ctrl

#should put these in separate files.... and maybe using's? 

# 2) Labeling Utils
# 3) Labeling funcs (from labeling.jl)
# 4) Plotting Functions
# 5) Battery Functions
# 6) MILP Functions
# 7) Optimal Control utils
# 7) Merge Algo and Utils
# 7) Merge testing....



###########################################
## 2: Label Utils


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
        gt_idx = gen_track[PL][gt_idx][2]
        push!(genOut, gen_now)
        
        PL -= 1
    end
    return genOut
end



function get_heur_astar(E, locs, Fvec)
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

function get_heur_label(Fvec, graph, C, E, heur_astar)
    
    function heur_label!(i)
        #if needed, run astar.  then, add all nodes in path to Fvec. return only f cost for this label, while altering all of Fvec
        if isnan(Fvec[i])
            astar_out = a_star(graph, i, E, C, heur_astar)
            if isempty(astar_out)
                path_list = [i] 
                printstyled("No path to goal!! \n", color=:light_red)
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

function get_heur_label_euc(Fvec, locs, E)
    function heur_label!(i)
        if isnan(Fvec[i])
            heur = norm(locs[i,:] - locs[E,:]) 
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
    #
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

function make_graph(euc_inst)
    A = euc_inst.Alist
    N = length(A)
    g = SimpleGraph(N)
    for i in 1:N
        for j in A[i]
            add_edge!(g, i, j)
            # add_edge!(g, j, i)
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


###########################################
## 3: Label functions from labeling.jl files
function get_heur_astar(E, locs, Fvec)
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

function get_heur_label(Fvec, graph, C::SparseMatrixCSC{Int64, Int64}, E, heur_astar) #def if discrete version
    function heur_label!(i)
        #if needed, run astar.  then, add all nodes in path to Fvec. return only f cost for this label, while altering all of Fvec
        if Fvec[i] == 	2^63 - 1
            astar_out = a_star(graph, i, E, C, heur_astar)
            if isempty(astar_out)
                path_list = [i] 
                printstyled("No path to goal!! \n", color=:light_red)
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

function get_heur_label_euc_disc(Fvec, locs, E) #def if we're using euc distance
    function heur_label!(i)
        if Fvec[i] == 2^63 - 1
            heur = Int(floor(norm(locs[i,:] - locs[E,:]))) 
            Fvec[i] = heur
        end
        return Fvec[i]
    end
    return heur_label!
end


function get_heur_label(Fvec, graph, C, E, heur_astar) #definition if we are getting a star value
    
    function heur_label!(i)
        #if needed, run astar.  then, add all nodes in path to Fvec. return only f cost for this label, while altering all of Fvec
        if isnan(Fvec[i])
            astar_out = a_star(graph, i, E, C, heur_astar)
            if isempty(astar_out)
                path_list = [i] 
                printstyled("No path to goal!! \n", color=:light_red)
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


function get_heur_label_euc(Fvec, locs, E) #def if we're using euc distance
    function heur_label!(i)
        if isnan(Fvec[i])
            heur = norm(locs[i,:] - locs[E,:]) 
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



###########################################
## 4: Plotting Functions 
function plot_euc_graph(euc_inst; path = [], gen = [], color_ends = true)
    #make Graph() then just graph plot
    set_default_plot_size(20cm, 20cm)
    g = make_graph(euc_inst)
    N = length(euc_inst.Alist)
    edge_index(e::Edge) = edgemap[e]
    nE = ne(g)
    edge_colors = [colorant"gray" for i in 1:nE]
    node_colors = [colorant"gray", colorant"cyan", colorant"magenta", colorant"indianred3", colorant"magenta"]
    node_labels = ones(Int, N)
    nodesize = ones(N)
    if !isempty(path)
        edgelist = collect(edges(g))
        edgemap = Dict{Edge, Int}()
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
        end
        Gsummed = sum(euc_inst.GFlipped, dims = 1)
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
            src, dst = e.src, e.dst
            if Gsummed[src] >= 1 || Gsummed[dst] >= 1
                edge_colors[i] = colorant"red"
                if Gsummed[src] >= 1
                    node_labels[src] = 4
                end
                if Gsummed[dst] >= 1
                    node_labels[dst] = 4
                end
            end

        end
        for k in 2:length(path)
            i = path[k-1]
            j = path[k]
            edge_idx = edge_index(Edge(i,j))
            nodesize[i] = 2
            nodesize[j] = 2
            if gen[k-1] == 0
                edge_colors[edge_idx] = node_colors[2]
                k == 2 && (node_labels[i] = 2)
                node_labels[j] = 2
            elseif gen[k-1] == 1
                edge_colors[edge_idx] = node_colors[3]
                k == 2 && (node_labels[i] = 3)
                node_labels[j] = 3
            end
        end
        
    else
        edgelist = collect(edges(g))
        edgemap = Dict{Edge, Int}()
        Gsummed = sum(euc_inst.GFlipped, dims = 1)
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
            src, dst = e.src, e.dst
            if Gsummed[src] >= 1 || Gsummed[dst] >= 1
                edge_colors[i] = colorant"red"
                if Gsummed[src] >= 1
                    node_labels[src] = 4
                end
                if Gsummed[dst] >= 1
                    node_labels[dst] = 4
                end
            end

        end
        #get vector of edges, label them 1 (grey) or 3 (red)
        
    end
    if color_ends
        node_labels[euc_inst.S] = 5
        node_labels[euc_inst.E] = 5
    end
    nodefillc2 = node_colors[node_labels]
    # path_plt = gplot(g, x_locs, y_locs, edgestrokec = edge_colors, nodefillc = nodefillc2)
    
    plt = gplot(g, euc_inst.locs[:,1], euc_inst.locs[:,2], edgestrokec = edge_colors, nodefillc = nodefillc2, nodesize=nodesize)
    return plt

    #also get a legend object.... to add to graph plot

end
function plot_euc_graph_labeled(euc_inst::EucGraph; label_strings::Vector{String}, label_units::Vector{String} = [""], label_vals::Matrix{Float64}, label_edge_idxs::Vector{Int}, path = [], gen = [], color_ends = true)
    #make Graph() then just graph plot
    set_default_plot_size(20cm, 20cm)
    g = make_graph(euc_inst)
    N = length(euc_inst.Alist)
    edge_index(e::Edge) = edgemap[e]
    nE = ne(g)
    edge_colors = [colorant"gray" for i in 1:nE]
    node_colors = [colorant"gray", colorant"cyan", colorant"magenta", colorant"indianred3", colorant"magenta"]
    node_labels = ones(Int, N)
    nodesize = ones(N)
    edge_labels = fill("", nE)
    if !isempty(path)
        edgelist = collect(edges(g))
        edgemap = Dict{Edge, Int}()
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
        end
        Gsummed = sum(euc_inst.GFlipped, dims = 1)
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
            src, dst = e.src, e.dst
            if Gsummed[src] >= 1 || Gsummed[dst] >= 1
                edge_colors[i] = colorant"red"
                if Gsummed[src] >= 1
                    node_labels[src] = 4
                end
                if Gsummed[dst] >= 1
                    node_labels[dst] = 4
                end
            end

        end
        for k in 2:length(path)
            i = path[k-1]
            j = path[k]
            edge_idx = edge_index(Edge(i,j))
            nodesize[i] = 2 
            nodesize[j] = 2
            if gen[k-1] == 0
                edge_colors[edge_idx] = node_colors[2]
                k == 2 && (node_labels[i] = 2)
                node_labels[j] = 2
            elseif gen[k-1] == 1
                edge_colors[edge_idx] = node_colors[3]
                k == 2 && (node_labels[i] = 3)
                node_labels[j] = 3
            end
        end
        
    else
        edgelist = collect(edges(g))
        edgemap = Dict{Edge, Int}()
        Gsummed = sum(euc_inst.GFlipped, dims = 1)
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
            src, dst = e.src, e.dst
            if Gsummed[src] >= 1 || Gsummed[dst] >= 1
                edge_colors[i] = colorant"red"
                if Gsummed[src] >= 1
                    node_labels[src] = 4
                end
                if Gsummed[dst] >= 1
                    node_labels[dst] = 4
                end
            end

        end
    
        #get vector of edges, label them 1 (grey) or 3 (red)
     
    end
    if color_ends
        node_labels[euc_inst.S] = 5
        node_labels[euc_inst.E] = 5
    end
    nodefillc2 = node_colors[node_labels]
    for j in 1:length(label_edge_idxs)
        e = label_edge_idxs[j]
        stringy = ""
        Nlines = length(label_strings)
        for el in 1:Nlines
            stringy = string(stringy, string(label_strings[el], ": ", label_vals[j,el], " (", label_units[el], ")"  ))
            el == Nlines && continue
            stringy = string(stringy, "\n" )
        end
        edge_labels[e] = stringy
    end
    # path_plt = gplot(g, x_locs, y_locs, edgestrokec = edge_colors, nodefillc = nodefillc2)
    
    plt = gplot(g, euc_inst.locs[:,1], euc_inst.locs[:,2], edgestrokec = edge_colors, nodefillc = nodefillc2, nodesize=nodesize, edgelabel = edge_labels)
    return plt


end



function get_boxplot_plt(Nvec, times; color = "blue", xlabel="", xmax = Nvec[end])
    K = size(times,2)
    instance = collect(1:K)

    grouping = []
    for n in Nvec
        append!(grouping, fill(n, K))
    end


    
    df = DataFrame(x = repeat(instance, outer = [length(Nvec)]),
    time = vec(times'),
    group = grouping)
    set_default_plot_size(15cm, 7.5cm)
    plt = plot(df, x=:group, y = "time", Geom.boxplot, Scale.y_log10,
                Guide.xlabel(xlabel),
                Guide.ylabel("Time To Solve (s)"),
                Scale.x_continuous(format= :plain, maxvalue=xmax),
                Theme(default_color = color, highlight_width = 0pt, middle_color = fcolor, middle_width = 0.5pt))
    return plt
end
function fcolor(color)
    return colorant"white"
end


# make function to plot 3D graph (handmade with Gadfly... graphplot does not support this?)
# using Plots, PyPlot
function plot_euc_3D(euc_inst; path =[], gen = [])
    set_default_plot_size(90cm, 90cm)
    g = make_graph(euc_inst)
    N = length(euc_inst.Alist)
    locs = euc_inst.locs

end

#make function to plot a graph....
function plot_hybrid_soln(mapdef, path, genY)
    A = map_def.A
    N = size(A,1)
    g = SimpleGraph(N)

    G = map_def.G
    y_locs = -map_def.GPS_locs[:,1]
    x_locs = map_def.GPS_locs[:,2]

    for r in 1:size(A,1)
        row = A[r,:]
        for i in 1:length(row)
            if row[i] == 1
                add_edge!(g,r,i)
                add_edge!(g,i,r)

            end
        end
    end
    edgelist = collect(edges(g))
    edgemap = Dict{Edge, Int}()
    for (i,e) in enumerate(edgelist)
        edgemap[e] = i
        edgemap[reverse(e)] = i
    end

    edge_index(e::Edge) = edgemap[e]
    nE = ne(g)
    edge_colors = [colorant"lightgrey" for i in 1:nE]
    node_colors = [colorant"lightgrey", colorant"black", colorant"magenta"]
    node_labels = ones(Int, N)
    for k in 2:length(path)
        i = path[k-1]
        j = path[k]
        edge_idx = edge_index(Edge(i,j))
        if genY[k-1] == 0
            edge_colors[edge_idx] = colorant"black"
            k == 2 && (node_labels[i] = 2)
            node_labels[j] = 2
        elseif genY[k-1] == 1
            edge_colors[edge_idx] = colorant"magenta"
            k == 2 && (node_labels[i] = 3)
            node_labels[j] = 3
        end
    end
    nodefillc2 = node_colors[node_labels]
    path_plt = gplot(g, x_locs, y_locs, edgestrokec = edge_colors, nodefillc = nodefillc2)
    return path_plt

end


function get_sol_vec(Nvec, prob_title; K = 10, conn = "_4conn", type = "euc", algo = "", prob = "DP", heur = "" )
    times = zeros(length(Nvec),K)
    avg_times = zeros(length(Nvec))
    nidx = 0
    for n in Nvec
        nidx += 1
        for k = 1:K
            if prob == "MILP"
                @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)" tMILP
                time_i = tMILP
            else
                @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)$(heur)" tdp
                time_i = tdp
            end
            times[nidx,k] = time_i
        end
        avg_times[nidx] = mean(times[nidx,:])
    end

    return times, avg_times
end


function get_pathlength_vec(Nvec, prob_title; K = 10, conn = "_4conn", type = "euc", algo = "", prob = "DP", heur = "" )
    lengths = zeros(length(Nvec),K)
    avg_lengths = zeros(length(Nvec))
    nidx = 0
    for n in Nvec
        nidx += 1
        for k = 1:K
            @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)$(heur)" pathL
            lengths[nidx,k] = length(pathL)
        end
        avg_lengths[nidx] = mean(lengths[nidx,:])
    end

    return lengths, avg_lengths
end
function get_avg_layer(Ncec, avg_times; color_in = "grey")
    layer_out = layer(x = Nvec, y = avg_times, Geom.point, Theme(default_color = color_in))
    return layer_out
end

############################################# 
## 5: Battery...
# raw data from Fitting the OCV-SOC relationship of a battery lithium-ion using genetic algorithm method

function get_OCV_func()
    SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    OCV=[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];

    OCV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return OCV_interp, 0, 0
end
function get_OCV_func_LiPo()
    SOC=100 .*  [1.0,
    0.9668412021343039,
    0.934531323117733,
    0.9045785889184145,
    0.8724856321908501,
    0.8410823292377329,
    0.7019553954577109,
    0.5677038583756852,
    0.4366584105242912,
    0.30752419166017014,
    0.2511283417275686,
    0.19575747375201058,
    0.14133157690308176]
    
    OCV= [ 16.788396797180177,
    16.62169780731201,
    16.48680353164673,
    16.392482948303222,
    16.33413486480713,
    16.282863426208497,
    16.039941120147706,
    15.684522533416748,
    15.49577431678772,
    15.29106364250183,
    15.126658964157105,
    14.949568367004394,
    14.800385189056396]

    OCV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return OCV_interp, minimum(SOC), maximum(SOC)
end
function get_OCV_func(SOC, OCV)
    # SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    # OCV=[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];

    OCV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return OCV_interp
end
function get_OCV_table()
    SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    OCV=[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];


    for i in 1:length(SOC)
        println("$(round(SOC[i],digits=4)) & $(round(OCV[i],digits=4))  \\\\ \\hline"   )
    end
end
function get_one_by_V()
    SOC=100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    OCV=1 ./[4.1617,4.0913,4.0749,4.0606,4.0153,3.9592,3.9164,3.8687,3.8163,3.7735,3.7317,3.6892,3.6396,3.5677,3.5208,3.4712,3.3860,3.2880,3.2037,3.0747];
    onebyV_interp = LinearInterpolation(reverse(SOC), reverse(OCV))
    return onebyV_interp
end



function get_one_by_Vsp()
    OCV, min, max = get_OCV_func()
    V(S,P) = (OCV(S) + sqrt( OCV(S)^2 - 4*P*70/1000 )  )/2

    S = 100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563]
    P = range(0, stop=20, length = length(S))
    S = reverse(S)
    Vdata = zeros(length(S), length(P))
    for i = 1:size(Vdata,1), j = 1:size(Vdata,2)
        Vdata[i,j] = 1/V(S[i], P[j])
    end
    
    itp = LinearInterpolation((S,P),Vdata)
    return itp

end
function get_one_by_Vsp(OCV, Svec, R, Pmax)
    V(S,P) = (OCV(S) + sqrt( OCV(S)^2 - 4*P*R )  )/2

    # S = 100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    P = range(0, stop=Pmax, length = length(Svec))
    Svec = reverse(Svec)
    Vdata = zeros(length(Svec), length(P))
    for i = 1:size(Vdata,1), j = 1:size(Vdata,2)
        Vdata[i,j] = 1/V(Svec[i], P[j])
    end
    
    itp = LinearInterpolation((Svec,P),Vdata)
    return itp

end

function get_one_by_VspGLM(OCV, Svec, R, Pmax)
    V(S,P) = (OCV(S) + sqrt( OCV(S)^2 - 4*P*R )  )/2

    # S = 100 .*[1,0.9503,0.9007,0.8510,0.8013,0.7517,0.7020,0.6524,0.6027,0.5530,0.5034,0.4537,0.4040,0.3543,0.3046,0.2550,0.2053,0.1556,0.1059,0.0563];
    P = range(0, stop=Pmax, length = length(Svec))
    Svec = reverse(Svec)
    Vdata = zeros(length(Svec), length(P))
    for i = 1:size(Vdata,1), j = 1:size(Vdata,2)
        Vdata[i,j] = 1/V(Svec[i], P[j])
    end
    data = DataFrame(X=repeat(Svec, outer = length(P)), Y=repeat(P, outer=length(Svec)), Z =vec(Vdata))
    gm1 = lm(@formula(Z ~ X + Y), data)

    return gm1
end

function model(t0, tf; SOC0 = 90)
    tv = [t0]
    SOCv = [SOC0]
    SOC = SOC0
    Vv = [V0]
    V = V0
    Δ = (tf - t0)/500
    t = t0
    while t < tf

    end
end
function SOC_dyn(t, SOC)
    SOC += Δ*It/Q
    x2  += Δ*(I - 35*Dp*x2/Rp^2)
    return SOC, x2
end
# Make funciton for Eulers Method propogating
function riemman(V, tvec, S0, P, Cmax; N = 50)
    S = S0
    t0, tf = tvec[1], tvec[2]
    Δ =  (tf - t0)/N
    t = t0
    S_vec = Float64[]
    t_vec = Float64[]
    append!(S_vec,S)
    append!(t_vec, t0)
    while true
        S -= P/(Cmax*V(S))*Δ*100 #times 100 for "percent"
        t += Δ

        append!(t_vec, t)
        append!(S_vec, S)
        t >= tf && return S_vec, t_vec
    end
    return 0
end
#function for ohmic drop...
function riemman2(Vf,OCV, tvec, S0, P, Cmax, Req; N = 50)
    S = S0
    t0, tf = tvec[1], tvec[2]
    Δ =  (tf - t0)/N
    t = t0
    S_vec = Float64[]
    t_vec = Float64[]
    append!(S_vec,S)
    append!(t_vec, t0)
    while true
        V = Vf(OCV, S, Req, P)
        S -= P/(Cmax*V)*Δ*100 #times 100 for "percent"
        t += Δ
        append!(t_vec, t)
        append!(S_vec, S)
        t >= tf && return S_vec, t_vec
    end
    return 0
end

function V_model2(OCV, SOC, Req, P)
    V = (OCV(SOC) + sqrt(OCV(SOC)^2 - 4*Req*P))/2
    return V
end

function riemman7(OCV, tvec, S0, P, Cmax, Req, C; N = 50)
    S = S0
    t0, tf = tvec[1], tvec[2]
    Δ =  (tf - t0)/N
    t = t0
    τ = Req*C
    S_vec = Float64[]
    t_vec = Float64[]
    V_vec = Float64[]
    U_vec = Float64[]
    Vk = OCV(S)
    Uk = exp(-Δ/τ)*0 +Req*(1-exp(-Δ/τ))*P/Vk
    Vk1 = 0
    Uk1 = 0
    append!(S_vec,S)
    append!(t_vec, t0)
    append!(V_vec, Vk)
    append!(U_vec, Uk)

    while true
        Uk1 = exp(-Δ/τ)*Uk + Req*(1 - exp(-Δ/τ))*P/Vk
        Vk1 = OCV(S) - Uk1
        S -= P/(Cmax*Vk1)*Δ*100 #times 100 for "percent"
        t += Δ

        append!(t_vec, t)
        append!(S_vec, S)
        append!(V_vec, Vk1)
        append!(U_vec, Uk1)
        Uk = Uk1
        Vk = Vk1
        t >= tf && return S_vec, t_vec
    end
    return 0

end

function simps(xin::Vector, yin::Vector) #simpsons rule to get total charge used given current and time samples
    x,y = xin, yin
    n = length(y)-1

    n % 2 == 0 || (push!(x, x[end]  + (x[end] - x[end-1])/100); push!(y, 0); n += 1)
    length(x)-1 == n || error("`x` and `y` length must be equal")
    h = (x[end]-x[1])/n
    s = sum(y[1:2:n] + 4y[2:2:n] + y[3:2:n+1])
    return h/3 * s
end


###########################################
## 6: MILP Functions 
#trying edge notation... not possible to add fuel constraints....

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
    SC = 0*def.StartCost
    G = .! GFlipped

    Bstart = mean(nonzeros(def.C)) + maximum(def.C)
    Qstart = 9999
    Bmax = Bstart
    N = length(Alist)
    bigM = 99
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 1, "CPXPARAM_Emphasis_Memory" => 1, "CPXPARAM_WorkMem" => 12000, "CPX_PARAM_TILIM" => tlim, "CPX_PARAM_WORKDIR" => "C:\\Users\\Student\\CPLEX_WDIR"))
    @variable(m, 1 >= x[i=1:N,j=1:N] >= 0) #  , Bin)
    @variable(m, 1>= g[i=1:N,j=1:N] >= 0) # , Bin)
 
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
            if j ∉ Alist[i] == 0
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

###########################################
##  7: optimal control utils

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

function MILP_to_opt_ctrl(N, k; prob = "euc_probs2D", algo = "_label", conn = "", heur = "")
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


###########################################
##  8: Merge Function and Utils...
function array_to_list(X::Matrix{Float64})
    Xout = Vector{Float64}[]
    for i in 1:size(X,1)
        push!(Xout, X[i,:])
    end
    return Xout
end
function array_to_list(X::Matrix{Int64})
    Xout = Vector{Int64}[]
    for i in 1:size(X,1)
        push!(Xout, X[i,:])
    end
    return Xout
end
function dom2(X,Y)
    #returns true if arg1 >> arg2
    bool = false
    (X[1] <= Y[1] && X[2] <= Y[2]  ) && (bool = true) #&& X[3] <= Y[3]
    return bool 
end

function dom(X::Vector,Y::Vector)
    #returns true if X >> Y
    bool = false
    (X[1] <= Y[1] && X[2] >= Y[2]  && X[3] >= Y[3]) && (bool = true)
    return bool
end

function dom_min(X::Vector,Y::Vector)
    #returns true if X >> Y
    bool = false
    (X[1] <= Y[1] && X[2] <= Y[2]  && X[3] <= Y[3]) && (bool = true)  
    return bool 
end

function make_eff_list(N; maxes = [10,10,10])
    dim = length(maxes)
    M = zeros(N, )
end

function mergeXY(X::Vector{Vector{Float64}}, Y::Vector{Vector{Float64}})
    M = Vector{Float64}[]
    i,j = 1,1
    p,q = length(X), length(Y)
    while true
        if i>p
            M = vcat(M,  Y[j:end])
            return M
        end
        if j > q
            M = vcat(M, X[i:end])
            return M
        end
        Xi, Yj = X[i], Y[j]
        #3a
        if Xi[1] < Yj[1]
            while Xi[2] ≤ Yj[2] && Xi[1] < Yj[1]
                j += 1
                j > q && break 
                Yj = Y[j]
            end
            push!(M, Xi)
            i += 1
        #3b
        elseif Yj[1] < Xi[1]
            while Yj[2] ≤ Xi[2] && Yj[1] < Xi[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        #3c
        elseif Xi[2] < Yj[2]
            while Yj[2] ≥ Xi[2] && Xi[1] == Yj[1]
                j += 1
                j > q && break
                Yj = Y[j]
            end
            push!(M, Xi)
            i +=1   
        #3d
        else
            while Xi[2] ≥ Yj[2] && Xi[1] == Yj[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        end
    end
    return 0
end
function merge_sort3D(Temp, Perm)
    X, Y = Temp, Perm
    M = []
    i,j = 1,1
    p,q = length(X), length(Y)
    while true
        #2a
        if i>p
            M = vcat(M,  Y[j:end])
            return M
        #2b
        elseif j > q
            M = vcat(M, X[i:end])
            return M
        end
        Xi, Yj = X[i], Y[j]
        #3a
        if Xi[1] < Yj[1]
            while Xi[2] ≤ Yj[2] && Xi[1] < Yj[1]
                j += 1
                j > q && break 
                Yj = Y[j]
            end
            push!(M, Xi)
            i += 1
        #3b
        elseif Yj[1] < Xi[1]
            while Yj[2] ≤ Xi[2] && Yj[1] < Xi[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        #3c
        elseif Xi[2] < Yj[2]
            while Yj[2] ≥ Xi[2] && Xi[1] == Yj[1]
                j += 1
                j > q && break
                Yj = Y[j]
            end
            push!(M, Xi)
            i +=1   
        #3d
        elseif Xi[2] > Yj[2]
            while Xi[2] ≥ Yj[2] && Xi[1] == Yj[1]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
            
        # extra for third dimension: if xi1 == yj1 && xi2 == yj2
        elseif Xi[3] < Yj[3]
            while Yj[3] ≥ Xi[3] && Xi[1] == Yj[1] && Xi[2] == Yj[2]
                j += 1
                j > q && break
                Yj = Y[j]
            end
            push!(M, Xi)
            i += 1     
        
        else 
            while Xi[3] ≥ Yj[3] && Xi[1] == Yj[1] && Xi[2] == Yj[2]
                i += 1
                i > p && break
                Xi = X[i]
            end
            push!(M, Yj)
            j += 1
        end
    end
    return 0
end

function brute_merge_sort(X, Y)  #need to adjust for equivalent labels....
    #for equiv labels: know that x == y cannot have duplicates in X or Y.. so just remove 1...
    M = []
    eq_keep = []
    for y in Y
        add_y = true
        for x in X
            if dom_min(x,y)
                add_y = false
                if x == y #if equal, may need to keep 1 
                    #add to y to eq_keep, add_y = true
                    if y ∉ eq_keep #if not in eq_keep, then we have not seen this label before so add it to eq_keep and set bool to true
                        push!(eq_keep, y)
                        add_y = true
                    else #else we have seen in beore, so set add_y to false 
                        add_y = false
                    end
                end
            end
        end
        add_y && push!(M, y) 
    end
    #we do not need to do equality corrections when looping through other "direction"
    #..because we only need to keep 1 of each equivalent set
    for x in X
        keep_x = true 
        for y in Y
            if dom_min(y,x)
                keep_x = false
            end
        end
        keep_x && push!(M, x)
    end
    return M
end


function merge_3D_drew(X, Y) #BROKEN LOGIC!!!!!!
    #need to eventually change this to X::heap and Y::heap
    #1) Sort [X;Y] by C
    M = vcat(X, Y)
    sort!(M, by = x -> x[1])
    println("X cat Y:")
    display(M)
    #2) remove dominated labels
    Bmin = Inf
    Qmin = Inf
    bool = ones(Bool, size(M,1))
    done = zeros(Bool, size(M,1)) #use this to skip in main loop after looping through equal C values
    for idx in 1:length(M)
        done[idx] == 1 && continue
        label = M[idx]
        Bi, Qi = label[2], label[3]
        idx2 = idx + 1
        if idx2 <= size(M,1) && M[idx2][1] == label[1]
            while true
                if idx2 == size(M,1) || M[idx2+1,1] != label[1]
                    break
                end
                idx2 += 1
            end
            println("here")
            temp = M[idx:idx2]
            booltemp = EFF_X(temp, dom_min)
            #set bool slice to bool2 vals
            bool[idx:idx2] = booltemp
            done[idx:idx2] .= 1

            #need to grab new Bmin Qmin...
        elseif Bi >= Bmin && Qi >= Qmin #if inneficient
            #then remove this
            bool[idx] = 0
        else #else is dominated
            Bmin = minimum([Bi, Bmin])
            Qmin = minimum([Qi, Qmin]) 
        end

    end
    display(bool)
    return M[bool]
end


function merge_3D_old(X, Y) #BROKEN LOGIC!!!! FAILS IF EQUAL B AND Q BUT PRIOR C EXISTS
    Mtemp = vcat(X, Y)
    sort!(Mtemp) #sorts by 1st el then 2nd el then 3rd el....
    boolv = ones(Bool, length(Mtemp))
    #now go through each element 
    prior = (Inf, Inf, Inf)
    for i in eachindex(Mtemp)
        label =  Mtemp[i]
        if prior[1] == label[1] 
            if prior[2] == label[2] #we sorted, so if equal 1's and 2's then 3piror <= 3i
                boolv[i] = 0
            else #if 2's are not equal, then can just compare 
                if label[3] < prior[3] #if this label is efficient then update _prior_
                    #1'equal, 2's not.  We sorted, so if _label_[3] < prior[3] then we are eff
                    prior = copy(label)
                else #else if 3rds are equal or _label_3 is more, then set _boolv_[i] to 0 and keep prior as is
                    boolv[i] = 0
                end
            end
        else #if not equal C's, loop through every prior _i_ and check if dommed
            for ii in 1:(i-1)
                boolv[ii] == 0 && continue #if already dommed don't need to compare
                if dom_min(Mtemp[ii], label)
                    boolv[i] = 0
                    break
                end
            end
            if boolv[i] ==  1 #if this label is effieicnet across all others, then this is new _prior_
                prior = copy(label)
            end
        end
    end
    return Mtemp[boolv]

end

function brute_EFF_merge(X, Y)
    M = vcat(X, Y)
    sort!(M)
    boolM = EFF_X(M, dom_min)
    M = M[boolM]
    return M
end
function EFF_X(X, dommy)
    Xout = Vector{Number}[]
    eq_remove = Int[]
    eq_keep = Int[]
    unique!(X)
    bool_vec = ones(Bool, length(X))
    hard_dom = zeros(Bool, length(X))
    for i in 1:length(X)
        x = X[i]
        boolx = true
        for j in 1:length(X)
            xx = X[j] 
            i == j && continue
            if dommy(xx, x) && x != xx #if dominated AND not equal, set as hard dom
                boolx = false
                bool_vec[i] = false
                hard_dom[i] = 1
            end
            if x == xx #if they are equal label values, may need to turn boolx back on
                if i ∉ eq_keep && i ∉ eq_remove && hard_dom[i] == 0 #if we have never seen i before AND not hard dommed  then keep Xi
                    boolx = true; bool_vec[i] = true
                    println("first: ",i)
                    push!(eq_keep, i)
                    push!(eq_remove, j)
                elseif i ∈ eq_keep && hard_dom == 0
                    boolx = true; bool_vec[i] = true
                    push!(eq_remove, j)
                else 
                    boolx = false
                    bool_vec[i] = false
                end
            end
        end
    end
    return bool_vec
end
function EFF_list(Γ::BinaryMinHeap{Tuple{Float64, Vector{Float64}}}, new::Vector{Float64})  #should not need this, but may be faster if we have less? like a merge sort?
    EFF_bool = true
    map_Γ = Γ.node_map
    for i in 1:length(map_Γ)  #γ in Γ
        map_Γ[i] == 0 && continue
        γ = Γ[i][2]
        if dom(γ, new)
            EFF_bool = false
            return EFF_bool
        end
    end
    return EFF_bool
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




#-------END OF MODULE---------------
end
#-------END OF MODULE---------------
