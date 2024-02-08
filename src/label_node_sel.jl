function hybrid_node_selection(def::EucGraphInt; heur::String = "astar")
    start, goal, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
    Bstart, Qstart = def.B0, def.Q0
    SC = def.StartCost
    Cmin = minimum(nonzeros(C))
    Qstart = 9999
    Bstart = Int(floor(5*mean(nonzeros(C))))
    Bmax = Bstart
    SC = 0
    locs = def.locs
    
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
        printstyled("invalid heuristic... breaking\n", color=:light_red)
        return 0,[0], Bool.([0])
    end
    if heur_label!(start) == 2^63 - 1     
        return -1,[0], Bool.([0])
    end
    node_queue = fill(Inf, N) # keep track of minimum candidate label
    node_queue[start] = 0+heur_label!(start)  #init node_queue with cost for S
    
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
        gen_bool=0
    )


    Q = [MyLabel[] for _ in 1:N]
    P = [Set{AbbreviatedLabel}()  for _ in 1:N]
    push!(Q[start], label_init)


    came_from = [Tuple{Int64, Int64}[] for _ in 1:N]
    push!(came_from[start], (0, 1)) #S is start... so we just add a dummy entry to [S]
    gen_track = [Tuple{Int64, Int64}[] for _ in 1:N]  #data struct to track genertor patterns 
    
    z=0
    while true #loop until get to end node, or Q is empty
        all(isempty.(Q)) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break) 

        nodei = argmin(node_queue)
        # print("nodei = $nodei\n")

        node_queue[nodei] = Inf
        if nodei == goal
            sort!(Q[goal])
            label_out = Q[goal][1]
            opt_cost =  label_out.gcost
            opt_path = get_path(label_out, came_from, start)
            opt_gen = get_gen(label_out, gen_track)
            return opt_cost, opt_path, opt_gen
        end
        Qi = copy(Q[nodei])
        Q[nodei] = MyLabel[]
        
        for label in Qi #dump Qi into closed set
            push!(P[nodei], abbreviated_label(label))
        end
        for nodej in Alist[nodei] #for each adjacent node, take all Q[i] and expand to j
            nodej==nodei && continue
            hj = heur_label!(nodej)
            new_labels = MyLabel[]

            #loop through Q[nodei]...
            for labeli in Qi
                path_to_i = get_path(labeli, came_from, def.S)
                nodej ∈ path_to_i  && continue 
                
                GenPrior = 1
                Fbin = F[nodei,nodej]

                #GEN ON
                new_batt_state = labeli.batt_state - C[nodei, nodej] * (1 - Fbin) + Z[nodei, nodej] - SC * (1 - GenPrior)
                new_gen_state = labeli.gen_state - Z[nodei, nodej]
                Gij = def.GFlipped[nodei, nodej]
                if Gij == 0 && new_batt_state >= 0 && new_gen_state >= 0 && new_batt_state <= Bmax
                    temp_new_label = MyLabel(
                        gcost=labeli.gcost + C[nodei, nodej],
                        fcost=labeli.gcost + C[nodei, nodej] + hj,
                        hcost=hj,
                        node_idx=nodej,
                        prior_node_idx=labeli.node_idx,
                        _hold_came_from_prior=labeli.came_from_idx,
                        came_from_idx=-1, #fill in later! (in update func)
                        pathlength=labeli.pathlength + 1,
                        _hold_gen_track_prior=labeli.gentrack_idx,
                        gentrack_idx=-1,  #fill in later!
                        gen_bool=1,
                        batt_state=new_batt_state,
                        gen_state=new_gen_state
                    )

                    if EFF_list(new_labels, temp_new_label) && EFF_P(P[nodej], abbreviated_label(temp_new_label))  #TODO need to rewrite this function
                        push!(new_labels, temp_new_label)
                    end
                end

                #GEN OFF
                new_batt_state = labeli.batt_state - C[nodei, nodej] * (1 - Fbin)
                if new_batt_state >= 0
                    temp_new_label = MyLabel(
                        gcost=labeli.gcost + C[nodei, nodej],
                        fcost=labeli.gcost + C[nodei, nodej] + hj,
                        hcost=hj,
                        node_idx=nodej,
                        prior_node_idx=labeli.node_idx,
                        _hold_came_from_prior=labeli.came_from_idx,
                        came_from_idx=-1, #fill in later! (in update func)
                        pathlength=labeli.pathlength + 1,
                        _hold_gen_track_prior=labeli.gentrack_idx,
                        gentrack_idx=-1,  #fill in later!
                        gen_bool=0,
                        batt_state=new_batt_state,
                        gen_state=labeli.gen_state
                    )
                    if EFF_list(new_labels, temp_new_label)  && EFF_P(P[nodej], abbreviated_label(temp_new_label)) 
                        push!(new_labels, temp_new_label)
                    end
                end
            end 
            
            isempty(new_labels) && continue
            #now correct paths and gens
            for lidx in 1:length(new_labels)
                temp_new_label = new_labels[lidx]
                new_label = update_path_and_gen!(temp_new_label, came_from, gen_track)
                new_labels[lidx] = new_label
            end
            Qj = merge_3D(Q[nodej], new_labels) #merge new labels into Q[nodej] 
            Q[nodej] = Qj

            #now change node_queue[j] value
            if !isempty(Q[nodej]) #may be no new labels...
                new_min = Q[nodej][1].fcost   #f cost of minimum cost label
                node_queue[nodej] = new_min
            elseif isempty(Q[nodej])
                node_queue[nodej] = Inf
            end
        end
        z+=1
        z ==200_000 && (printstyled("ZBREAK@$(z)", color=:light_red); break)
    end
    return 0,[0], Bool.([0])
end


"""
Old function.  No immutable structs...
"""
function hybrid_node_selection_dumb(def::EucGraphInt; heur::String = "astar")
    S, E, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
    Bstart, Qstart = def.B0, def.Q0
    SC = def.StartCost
    Cmin = minimum(nonzeros(C))
    Qstart = 9999
    Bstart = Int(floor(5*mean(nonzeros(C))))
    Bmax = Bstart
    SC = 0
    locs = def.locs
    
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
        printstyled("invalid heuristic... breaking\n", color=:light_red)
        return 0,[0], Bool.([0])
    end
    if heur_label!(S) == 2^63 - 1     
        return -1,[0], Bool.([0])
    end
    node_queue = fill(Inf, N) # keep track of minimum candidate label
    node_queue[S] = 0+heur_label!(S)  #init node_queue with cost for S
    i = Int(S)
    Q , P = [Vector{Int64}[] for _ in 1:N], [Vector{Int64}[] for _ in 1:N] 
    push!(Q[i], [0, Bstart, Qstart, S,     S, 1,     1, 1,     heur_label!(S) + 0])    
    # label [gc, B, Q, j, {{prior_i, came_from[j]pos}},{{|X|, gen_track pos}} f]
    
    # X = Vector{Int64}[ [S] ]    #hold path info
    # Y = Vector{Int64}[  []   ]  #hold gen info

    came_from = [Vector{Int64}[] for _ in 1:N]
    push!(came_from[S], [0, 1]) #S is start... so we just add a dummy entry to [S]

    gen_track = [Vector{Int64}[] for _ in 1:N]  #data struct to track genertor patterns 
    
    z=0
    # prog = ProgressUnknown("Working hard... $(title)", spinner=true)
    while true #loop until get to end node, or Q is empty
        all(isempty.(Q)) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break) 

        i = argmin(node_queue)
        # print("nodei = $i\n")
        node_queue[i] = Inf
        if i == E 
            sort!(Q[E])
            label_out = Q[E][1]
            opt_cost =  label_out[1]
            opt_path = get_path(label_out, came_from, def.S)
            opt_gen = get_gen(label_out, gen_track)
            return opt_cost, opt_path, opt_gen
        end
        Qi = copy(Q[i])
        Q[i] = Vector{Int64}[]
        
        for lidx in 1:length(Qi)
            push!(P[i], Qi[lidx])
        end
        for j in Alist[i] #for each adjacent node, take all Q[i] and expand to j
            j==i && continue
            hj = heur_label!(j)
            new_labels = Vector{Int64}[]

            #loop through Q[i]...
            for ii in 1:length(Qi)
                #for each q in Q, exapnd to j

                labeli = Qi[ii]

                fj = labeli[1] + C[i,j] + heur_label!(j)
                came_from_i = labeli[6] #not correct, need to adjust K index after the new_labels list is finalized
                path_to_i = get_path(labeli, came_from, def.S)

                j ∈ path_to_i  && continue 
                
                GenPrior = 1
                Fbin = F[i,j]
                if labeli[2]-C[i,j]*(1-Fbin) >= 0 #gen OFF
                    qj = Vector{Int64}(undef, 9)  
                    qj[1] = labeli[1] + C[i,j]
                    qj[2] = labeli[2] - C[i,j]*(1-Fbin)
                    qj[3] = labeli[3] + 0
                    qj[4] = j
                    qj[5] = i #prior i
                    qj[6] = came_from_i #this is came_from_iiiii not for _qj_, change this later in loop
                    qj[7] = labeli[7] + 1  #path length (in nodes)
                    qj[8] = labeli[8] #use later, then change to new gen_track_idx
                    qj[9] = 0 #save gen bool, correct later to be f cost
                    if EFF_list(new_labels, qj) 
                        push!(new_labels, qj)
                    end
                end
                if def.GFlipped[i,j] == 0 && labeli[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior) >= 0 && labeli[3]-Z[i,j] >= 0 && labeli[2]-C[i,j] + Z[i,j] <= Bmax #gen ON
                    qj = Vector{Int64}(undef, 9)  
                    qj[1] = labeli[1] + C[i,j] 
                    qj[2] = labeli[2] - C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior)
                    qj[3] = labeli[3] - Z[i,j]
                    qj[4] = j
                    qj[5] = i #came from i
                    qj[6] = came_from_i #this is came_from_iiii not for _qj_ --> correct this later in the loop after labels are pruned from _new_labels_
                    qj[7] = labeli[7] + 1 #path length (in nodes)
                    qj[8] = labeli[8] #use later, then change to new gen_track_idx
                    qj[9] = 1 #save genbool, replace with fj later...
                    if EFF_list(new_labels, qj) 
                        push!(new_labels, qj)
                    end
                end
            end 
            hj = heur_label!(j)
            isempty(new_labels) && continue

            for idx in 1:length(new_labels)
                qj = new_labels[idx]
                
                #correct path pointers:
                path_pointer = findall(x-> x == [qj[5], qj[6]], came_from[j])  
                if isempty(path_pointer) #if no other label has used this same path...
                    push!(came_from[j], [qj[5], qj[6]])   #[prior is i, in idx _k_ in came_from[i]]
                    new_labels[idx][6] = length(came_from[j])
                    #new_labels[idx][6] remains as prior i, but we don't need it...
                else #if another label already used this same path, then we use that same pointer
                    pointer_idx = path_pointer[1]
                    new_labels[idx][6] = pointer_idx #label now has index for path in came_from 
                    priori = qj[5]
                end
                
                #correct gen pointers to #[qj[9] (on/off), qj[8] (pointer back to gen_track[PL-1][poiter])
                PL = qj[7]
                gen_pointer = findall(x -> x == [qj[9], qj[8]], gen_track[PL])
                if isempty(gen_pointer)
                    push!(gen_track[PL], [qj[9], qj[8]])
                    new_labels[idx][8] = length(gen_track[PL])
                else
                    pointer_idx = gen_pointer[1]
                    qj[8] = pointer_idx #label now has index for generator pattern in gen_track
                end

                #correct fj 
                new_labels[idx][9] = new_labels[idx][1] + hj
            end
            Qj = merge_3D(Q[j], new_labels) #merge new labels into Q[j] 
            Q[j] = Qj

            #now change node_queue[j] value
            if !isempty(Q[j]) #may be no new labels...
                new_min = Q[j][1][9] #f cost of minimum cost label
                node_queue[j] = new_min
            elseif isempty(Q[j])
                node_queue[j] = Inf
            end
        end
        z+=1
        z == 200_000 && (printstyled("ZBREAK@$(z)", color=:light_red); break)
    end
    return 0,[0], Bool.([0])
end


