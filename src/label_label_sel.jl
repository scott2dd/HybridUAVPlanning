function hybrid_label_selection(def::EucGraphInt; heur::String="astar")
    start, goal, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
    Bstart, Qstart = def.B0, def.Q0
    Cmin = minimum(nonzeros(C))
    genpen = Cmin * 0
    Bstart = Int(floor(5 * mean(nonzeros(C))))
    Qstart = 9999
    Bmax = Bstart
    SC = def.StartCost
    SC = 0
    anchor = def.anchor_list  #if a grid problem, this is zeros(N)
    locs = def.locs
    Dim = def.locs

    N = length(Alist)
    graph = make_graph(def, false)
    Fvec = fill(2^63 - 1, N)
    heur_astar = get_heur_astar(goal, locs, Fvec)
    heur_astar(start)
    if heur == "astar"
        heur_label! = get_heur_label(Fvec, graph, C, goal, heur_astar)
    elseif heur == "euc"
        heur_label! = get_heur_label_euc(Fvec, locs, goal)
    else
        println("invalid heuristic... breaking")
        return return 0, [0], Bool.([0])
    end
    if heur_label!(start) == 2^63 - 1
        return -1, [0], Bool.([0])
    end

    label_init = MyLabel(
        batt_state=Bstart,
        gen_state=Qstart,
        gcost=0,
        fcost=heur_label!(start),
        hcost=heur_label!(start),
        node_idx=start,
        prior_node_idx=0,
        _hold_came_from_prior=0,
        came_from_idx=1,
        pathlength=1,
        _hold_gen_track_prior=0,
        gentrack_idx=1,
        gen_bool=0,
    )
    Q = MutableBinaryMinHeap{MyLabel}()
    push!(Q, label_init)

    P = [MyLabel[] for k = 1:N] #tuple of vectors of labels.... may change to hashmap???????? or set.... 

    came_from = [Tuple{Int64, Int64}[] for _ in 1:N]
    push!(came_from[start], (0, 9999)) #S is start... so we just add a dummy entry to [S]
    gen_track = [Tuple{Int64,Int64}[] for _ in 1:N]  #data struct to track genertor patterns 
    z = 0  #will not have a path larger than N
    while true 
        isempty(Q) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break)

        #pull minimum cost label....
        labelN = pop!(Q)
        nodei = labelN.node_idx
        GenPrior = 1 #label_treated[7]
        # P[nodei] = vcat(P[nodei], reshape(label_treated, (1, 9)))
        push!(P[nodei], labelN)

        if nodei == goal
            opt_cost = labelN.gcost
            opt_path = get_path(labelN, came_from, start)
            opt_gen = get_gen(labelN, gen_track)
            return opt_cost, opt_path, opt_gen
        end
        pathi = get_path(labelN, came_from, start)
        for nodej in Alist[nodei]
            nodej == nodei && continue
            nodej ∈ pathi && continue
            
            hj = heur_label!(nodej)
            Fbin = F[nodei, nodej] # if we can glide, then ALWAYS glide

            #GEN ON
            new_batt_state = labelN.batt_state - C[nodei, nodej] * (1 - Fbin) + Z[nodei, nodej] - SC * (1 - GenPrior)
            new_gen_state = labelN.gen_state - Z[nodei, nodej]
            Gij = def.GFlipped[nodei, nodej]
            if Gij == 0 && new_batt_state >= 0 && new_gen_state >= 0 && new_batt_state <= Bmax
                temp_new_label = MyLabel(
                    gcost = labelN.gcost + C[nodei, nodej],
                    fcost = labelN.gcost + C[nodei, nodej] + hj,
                    hcost = hj,
                    node_idx = nodej,
                    prior_node_idx = labelN.node_idx,
                    _hold_came_from_prior = labelN.came_from_idx,
                    came_from_idx = -1, #fill in later! (in update func)
                    pathlength = labelN.pathlength + 1,
                    _hold_gen_track_prior = labelN.gentrack_idx,
                    gentrack_idx = -1,  #fill in later!
                    gen_bool = 1,
                    batt_state = new_batt_state,
                    gen_state =  new_gen_state
                )
                if EFF_heap(Q, temp_new_label) && EFF_P(P, temp_new_label)
                    new_label = update_path_and_gen!(temp_new_label, came_from, gen_track)
                    push!(Q, new_label)
                end
            end

            #GEN OFF
            new_batt_state = labelN.batt_state - C[nodei, nodej] * (1 - Fbin)
            if new_batt_state >= 0 
                temp_new_label = MyLabel(
                    gcost = labelN.gcost + C[nodei, nodej],
                    fcost = labelN.gcost + C[nodei, nodej] + hj,
                    hcost = hj,
                    node_idx = nodej,
                    prior_node_idx = labelN.node_idx,
                    _hold_came_from_prior = labelN.came_from_idx,
                    came_from_idx = -1, #fill in later! (in update func)
                    pathlength = labelN.pathlength + 1,
                    _hold_gen_track_prior = labelN.gentrack_idx,
                    gentrack_idx = -1,  #fill in later!
                    gen_bool = 0,
                    batt_state = new_batt_state,
                    gen_state = labelN.gen_state
                )

                if EFF_heap(Q, temp_new_label) && EFF_P(P, temp_new_label)
                    new_label = update_path_and_gen!(temp_new_label, came_from, gen_track)
                    push!(Q, new_label)
                end
            end
        end
        z += 1
        z == 200_000 && (printstyled("ZBREAK@$(z)", color=:light_red); break)
    end
    return 0, [0], Bool.([0])
end



"""
    input EucGraphInt instance
    heur - optional argument: "euc" or "astar" (manhattan not working right now)
"""
function hybrid_label_selection_dumb(def::EucGraphInt; heur::String = "astar") 
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
    graph = make_graph(def, false)
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
