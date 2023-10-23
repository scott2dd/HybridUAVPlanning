#we want to wrtie a labeling algo which gets SPP from current node to goal as lower bound....
# "expensive"  to compute, but will this reduce the search space enough to make it worth it?
#as we go along, can we re-use Astar computations as we evaluate new nodes?.

#Opt path from S -> E (SPP) also contains optimal path from every node in that path (i,j,k) -> e
#so longer SPP computatiosn (farther from goal) will give us more info, can reuse.
# using ProgressMeter


##############################################
"""
    input EucGraphInt instance
    heur - optional argument: "euc" or "astar" (manhattan not working right now)
"""
function hybrid_label_selection(def::EucGraphInt; heur::String = "astar") 
    S, E, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
    Bstart, Qstart = def.B0, def.Q0
    Cmin = minimum(nonzeros(C))
    # Bstart = mean(nonzeros(C)) * 4
    genpen = Cmin *0
    Bstart = Int(floor(5*mean(nonzeros(C))))
    Qstart = 9999
    Bmax = Bstart
    SC = def.StartCost
    SC = 0
    anchor = def.anchor_list  #if a grid problem, this is zeros(N)
    locs = def.locs
    Dim = def.locs 
    
    N = length(Alist)
    graph = make_graph(def)
    Fvec = fill(2^63 - 1, N)
    heur_astar = get_heur_astar(E, locs, Fvec)
    heur_astar(S)
    if heur == "astar"
        heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
    elseif heur == "euc"
        heur_label! = get_heur_label_euc(Fvec, locs, E)
    else
        println("invalid heuristic... breaking")
        return return 0,[0], Bool.([0])
    end
    if heur_label!(S) == 2^63 - 1     
        return -1,[0], Bool.([0])
    end
    
    Q = MutableBinaryMinHeap([   (heur_label!(S) + 0, [0, Bstart, Qstart, S,    S,1,    1,1,   heur_label!(S) + 0]) ] )

    P = [zeros(Int, 0,9) for k = 1:size(C,1)] #set of treated labels | new labels compare to this AND Q

    came_from = [Vector{Int64}[] for _ in 1:N]
    push!(came_from[S], [0, 9999]) #S is start... so we just add a dummy entry to [S]
    gen_track = [Vector{Int64}[] for _ in 1:N]  #data struct to track genertor patterns 
    z=0  #will not have a path larger than N
    # prog = ProgressUnknown("Working hard...", spinner=true)
    while true #loop until get to end node, or Q is empty
        isempty(Q) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break)

        #pull minimum cost label....
        next = pop!(Q)
        label_treated = next[2]
        i = label_treated[4]
        # K = label_treated[5]
        GenPrior = 1 #label_treated[7]
        P[i] = vcat(P[i], reshape(label_treated,(1,9)))

        if i == E
            opt_cost =  label_treated[1]
            opt_path = get_path(label_treated, came_from, def.S)
            opt_gen = get_gen(label_treated, gen_track)
            return opt_cost, opt_path, opt_gen
        end
        pathi = get_path(label_treated, came_from, def.S)
        for j in Alist[i]
            j==i && continue
            j∈pathi && continue
            hj =  heur_label!(j) 
            Fbin = F[i,j] # if we can glide, then ALWAYS glide
            
            #GEN ON
            if def.GFlipped[i,j] == 0 && label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior) >= 0 && label_treated[3]-Z[i,j] >= 0 && label_treated[2]-C[i,j] + Z[i,j] <= Bmax
                label = Vector{Int64}(undef, 9)  
                label[1] = label_treated[1] + C[i,j]
                label[2] = label_treated[2] - C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior)
                label[3] = label_treated[3] -Z[i,j]
                label[4] = j
                label[5] = i
                label[6] = label_treated[6] #correct later
                label[7] = label_treated[7] + 1 #path length
                label[8] = label_treated[8] #correct later
                label[9] = 1 #gen bool, correct later
                
                
                if EFF_heap(Q, label) && EFF_P(P, label)
                    #correct path...
                    path_pointer = findall(x-> x == [label[5], label[6]], came_from[j])
                    if isempty(path_pointer) #if no other label has used this same path...
                        push!(came_from[j], [label[5], label[6]])   #[prior is i, in idx _k_ in came_from[i]]
                        label[6] = length(came_from[j])
                    else #if another label already used this same path, then we use that same pointer
                        pointer_idx = path_pointer[1]
                        label[6] = pointer_idx #label now has index for path in came_from 
                        priori = label[5]
                    end
                    #correct gen....
                    PL = label[7]
                    gen_pointer = findall(x -> x == [label[9], label[8]], gen_track[PL])  #label[9] is holding gen on/off
                    if isempty(gen_pointer)
                        push!(gen_track[PL], [label[9], label[8]])
                        label[8] = length(gen_track[PL])
                    else
                        pointer_idx = gen_pointer[1]
                        label[8] = pointer_idx #label now has index for generator pattern in gen_track
                    end

                    #correct hj
                    label[9] = label[1] + hj
                    #add to Q..
                    push!(Q, (label[1]+hj,label))

                end
            end
            if label_treated[2]-C[i,j]*(1-Fbin) >= 0 #GEN OFF 
                # label =  [gc, label_treated[2]-C[i,j]*(1-Fbin), label_treated[3],  j, Int(length(X)+1), gc+h]
                label = Vector{Int64}(undef, 9)  
                label[1] = label_treated[1] + C[i,j]
                label[2] = label_treated[2] - C[i,j]*(1-Fbin)
                label[3] = label_treated[3] + 0
                label[4] = j
                label[5] = i
                label[6] = label_treated[6] #correct later
                label[7] = label_treated[7] + 1 
                label[8] = label_treated[8] #correct later
                label[9] = 0 #gen bool, correct later


                if EFF_heap(Q, label) && EFF_P(P, label)
                    #correct path...
                    path_pointer = findall(x-> x == [label[5], label[6]], came_from[j])
                    if isempty(path_pointer) #if no other label has used this same path...
                        push!(came_from[j], [label[5], label[6]])   #[prior is i, in idx _k_ in came_from[i]]
                        label[6] = length(came_from[j])
                    else #if another label already used this same path, then we use that same pointer
                        pointer_idx = path_pointer[1]
                        label[6] = pointer_idx #label now has index for path in came_from 
                        priori = label[5]
                    end
                    #correct gen....
                    PL = label[7]
                    gen_pointer = findall(x -> x == [label[9], label[8]], gen_track[PL])  #label[9] is holding gen on/off
                    if isempty(gen_pointer)
                        push!(gen_track[PL], [label[9], label[8]])
                        label[8] = length(gen_track[PL])
                    else
                        pointer_idx = gen_pointer[1]
                        label[8] = pointer_idx #label now has index for generator pattern in gen_track
                    end

                    #correct hj
                    label[9] = label[1] + hj
                    push!(Q, (label[9],label))
                end
            end
        end
        z+=1
        z == 200_000 && (printstyled("ZBREAK@$(z)", color=:light_red); break)
        # z%500 == 0 && ProgressMeter.next!(prog)
    end
    return 0,[0], Bool.([0])
end


##
# @load "Problems\\euc_probs_disc\\550_4conn_9" euc_inst
# @load "Problems\\lattice_probs_2D\\50_1" lattice_inst
# @time hybrid_label_selection(euc_inst)

##
# @load "Problems\\lattic/e_probs\\20_1" lattice_inst
# @load "Problems\\euc_probs\\50_4conn_1" euc_inst
# @time costL, pathL, genL, lX, look, dist = hybrid_with_LB(lattice_inst, heur = "astar")
# set_default_plot_size(15cm, 10cm)
# plot(x = 1:length(dist), y = dist, Geom.line, Guide.xlabel("Iterations"), Guide.ylabel("Distance to Goal"))
# Alist = euc_inst.Alist
# N = length(Alist)
# # println("Solving with $(N) Nodes || B0 = $(Bstart)")
# graph = make_graph(def)
# Fvec = fill(NaN, N)
# heur_astar = get_heur_astar(E, locs, Fvec)
# heur_astar(euc_inst.S)
# heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
# heur_label!(euc_inst.S)
# ##
# def = euc_inst
# S, E, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
# Bstart, Qstart = def.B0, def.Q0
# Cmin = minimum(nonzeros(C))
# # Bstart = mean(nonzeros(C)) * 4
# genpen = Cmin*.0
# Bstart = 9999+mean(nonzeros(C)) + 5*maximum(C)
# Qstart = 9999
# Bmax = Bstart
# SC = def.StartCost
# SC = 0
# anchor = def.anchor_list  #if a grid problem, this is zeros(N)
# locs = def.locs
# # G = .!def.GFlipped
# Dim = def.locs           ####### pass third argument of heuristics as locs

# N = length(Alist)
# # println("Solving with $(N) Nodes || B0 = $(Bstart)")
# graph = make_graph(def)
# Fvec = fill(NaN, N)
# heur_astar = get_heur_astar(E, locs, Fvec)
# heur_astar(S)
# heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
# if heur_label!(S) == Inf     
#     println("Infinite Start Heuristic...")
#     return 0,0,0,0,0,0,0
# end
##now test labeling algo with astar LB.....

# @load "Problems\\euc_probs\\1990_8conn" euc_inst
# @time out = hybrid_with_LB(euc_inst)
# println("Labels made: $(out[end])")



# @load "Problems\\euc_probs2D\\50_4conn_1" euc_inst
# cost, path, gen, lX, look = hybrid_with_LB(euc_inst)
# @load "Solutions\\euc_probs2D\\50_4conn_1" tdp cost path gen lX
# pathA1 = a_star(graph, euc_inst.S, euc_inst.E, euc_inst.C)
# pathA = [pathA1[i].src for i in 1:length(pathA1)]
# append!(pathA, pathA1[end].dst)
# plot_euc_graph(euc_inst, path = path, gen = gen)
# plot_euc_graph(euc_inst, path = pathA, gen = gen)




######################################################################################


## first do an Astar test with Graphs.jl package
# @load "Problems\\euc_probs\\30_4conn" euc_inst

# Alist = euc_inst.Alist
# S, E = euc_inst.S, euc_inst.E
# C= euc_inst.C
# locs = euc_inst.locs
# Fvec = fill(NaN, length(Alist))
# plt = plot_euc_graph(euc_inst)
# heur_astar = get_heur_astar(E, locs, Fvec)
# graph = make_graph(euc_inst)
# plot_euc_graph(euc_inst)

# @time path_astar = a_star(graph, S, E, C, heur_astar)

# path_list, LB_vec, cost = astar_proc(path_astar, C)
# println("from $(S) to $(E)")
# println("cost: $(cost)")
# display(path_list)

# ##
# Fvec = fill(NaN, length(Alist))
# heur_astar = get_heur_astar(E, locs, Fvec)
# heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
# heur_label!(euc_inst.S)


# ## time getting SPP for every node...
# function get_total_heur_time(euc_inst)
#     Fvec = fill(NaN, length(euc_inst.Alist))
#     heur_astar = get_heur_astar(euc_inst.E, euc_inst.locs, Fvec); heur_astar(euc_inst.S)
#     graph = make_graph(euc_inst)
#     heur_label! = get_heur_label(Fvec, graph, euc_inst.C, euc_inst.E, heur_astar)
#     heur_label!(1)
#     total_time = 0
#     for i in 1:length(euc_inst.Alist)
#         tdp = @elapsed h = heur_label!(i)
#         total_time += tdp
        
#     end
#     return total_time
# end


# Nvec = 30:40:2000
# tvec4 = zeros(length(Nvec))
# tvec8 = zeros(length(Nvec))

# cnt = 0
# for n in Nvec
#     println(n)
#     cnt += 1
    
#     conn = 4
#     @load "Problems\\euc_probs\\$(n)_$(conn)conn" euc_inst
#     t = get_total_heur_time(euc_inst)
#     tvec4[cnt] = t

#     conn = 8
#     @load "Problems\\euc_probs\\$(n)_$(conn)conn" euc_inst
#     t = get_total_heur_time(euc_inst)
#     tvec8[cnt] = t
    
# end
# plt = Gadfly.plot(
#     layer(x = Nvec, y = tvec4, Geom.line, Geom.point, Theme(default_color = "grey")),
#     layer(x = Nvec, y = tvec8, Geom.line, Geom.point, Theme(default_color = "red")),
#     Guide.manual_color_key("", ["branch = 4", "branch = 8"], ["grey", "red"]),
#     Guide.title("Time to get A-star LB for every Node in Graph"),
#     Guide.xlabel("Numer of Nodes"),
#     Guide.ylabel("Time (s)")
#            )

# ##  PLot the graph....
# graph = make_graph(euc_inst)
# plot_euc_graph(euc_inst)

