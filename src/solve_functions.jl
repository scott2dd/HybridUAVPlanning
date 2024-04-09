"""
Function to solve euclidean instances, need to specify type.
Input:  algo::"node" or "label:
        dims:: "2D" or "3D"
        heur:: "astar" or "euc"
        
        output: none, saves results to files.

        Need to incorporate a time limit to terminate... if problems are taking too long we just end it early (don't go to end of Nvec)
        Then save last size we finished at....
        
        adding multithreading (solving individual problems in parallel.....)
"""
function solve_euc(;algo::String = "label", dims::String="2D", heur::String = "astar", tlim = 3600, Nstart = 50, conn = 0)
    Nvec = [50:500:2000; 2000:1000:20000]
    Nvecwhole = [50:500:2000; 2000:1000:20000]
    nidx = 0
    Nstart > 1900 && (nidx = findall(x->x==Nstart, Nvec)[1] - 1; Nvec = Nstart:1000:20000)

    if algo == "label"
        algof = hybrid_label_selection
        algo_tag = "_label"
    elseif algo == "node"
        algof = hybrid_node_selection
        algo_tag = "_node"
    else
        error("invalid algorithm selection")
    end

    if dims == "2D"
        if conn == 0 
            conn = ""; prob = "euc_probs2D"
        elseif conn != 0
            conn = "_$(conn)conn"
            prob = "connectivity_expr"
            Nvec = [50:1000:2000; 2000:2000:20000]
            Nvecwhole = [50:1000:2000; 2000:2000:20000]
        end
    elseif dims == "3D"
        conn = "_4conn"
        prob = "euc_probs_disc"
    else
        error("invalid problem Dimensions.  Use \"2D\" or \"3D\" ")
    end

    if heur == "astar"
        heur_tag = ""
    elseif heur == "euc"
        heur_tag = "_eucLB"
    elseif heur == "manhattan"
        error("manhattan distance not working right now -_-")
    else 
        error("invalud heuristic.  Use \"astar\" or \"euc\" or \"manhattan\"")
    end

    #run 1 to compile
    if prob == "connectivity_expr"
        @load "Problems/$(prob)/50_1$(conn)" euc_inst
    else
        @load "Problems/$(prob)/50$(conn)_1" euc_inst
    end
    tdp = @elapsed cost, path, gen = algof(euc_inst, heur = heur)
    
    nidx = 0
    printstyled("\n Solving Euclidean $(dims) Problems  || h(i): $(heur) || $(algo) \n ", color=:light_green)
    for n in Nvec
        nidx += 1
        printstyled("N = $(n) ||", color = :light_green)
        times_vec = Float64[]
        Zbreak_count = 0
        Threads.@threads for k = 1:10
            print(" $(k)")
            if prob == "connectivity_expr"
                @load "Problems/$(prob)/$(n)_$(k)$(conn)" euc_inst
            else
                @load "Problems/$(prob)/$(n)$(conn)_$(k)" euc_inst
            end
            tdp = @elapsed cost, pathL, gen = algof(euc_inst, heur = heur)
            if prob == "connectivity_expr"
                @save "Solutions/$(prob)/$(n)_$(k)$(conn)$(algo_tag)$(heur_tag)" tdp cost pathL gen
            else
                @save "Solutions/$(prob)/$(n)$(conn)_$(k)$(algo_tag)$(heur_tag)" tdp cost pathL gen
            end
            push!(times_vec, tdp)
            cost == 0 && (Zbreak_count += 1)
            Zbreak_count > 3 && (nidx -= 1;  break) 
        end
        if mean(times_vec) > tlim || Zbreak_count > 3
            if Zbreak_count > 3
                printstyled("STOPPED EARLY -  Z break \n", color=:light_red) 
                nwholeidx = findall(x->x==n, Nvecwhole)[1] - 1 #minus 1, these dont count
            elseif mean(times_vec) > tlim 
                printstyled("STOPPED EARLY -  Z break \n", color=:light_red) #stop here because we solved them all
                nwholeidx = findall(x->x==n, Nvecwhole)[1] 
            end
            n = Nvecwhole[nwholeidx]
            @save "Solutions/END_$(prob)$(algo_tag)$(heur_tag)" n #save where we ended early.....
            return 0
        end
        meantts = round(mean(times_vec), digits = 2)
        println(" ||  mTTS (s) : $(meantts) ")
    end
    printstyled("\n SOLVED --- Euclidean $(dims) Problems  || h(i): $(heur) || $(algo) --- \n", color=:light_green)
    n = Nvecwhole[end]
    @save "Solutions/END_$(prob)$(algo_tag)$(heur_tag)" n
end 

"""
Function to solve lattice instances, need to specify type.
Input:  algo::"node" or "label:
        dims:: "2D" or "3D"
        heur:: "astar" or "euc"
        
        output: none, saves results to files.

        Time limit to terminate... if problems are taking too long we just end it early (don't go to end of Nvec)
        Save last size we finished at....
        adding threading support....
"""
function solve_lattice(;algo::String = "label", dims::String="2D", heur::String = "astar", tlim = 3600, Nstart = 5)
    Nvec = 5:50
    Nvecwhole = 5:50
    nidx = 0
    Nstart > 5 && (nidx = findall(x->x==Nstart, Nvec)[1] - 1; Nvec = Nstart:50)

    if algo == "label"
        algof = hybrid_label_selection
        algo_tag = "_label"
    elseif algo == "node"
        algof = hybrid_node_selection
        algo_tag = "_node"
    else
        error("invalid algorithm selection: \"node\" or \"label\"")
    end

    if dims == "2D"
        prob = "lattice_probs_2D"
    elseif dims == "3D"
        prob = "lattice_probs_disc"
    else
        error("invalid problem Dimensions.  Use \"2D\" or \"3D\" ")
    end

    if heur == "astar"
        heur_tag = ""
    elseif heur == "euc"
        heur_tag = "_eucLB"
    elseif heur == "manhattan"
        error("manhattan distance not working right now")
    else 
        error("invalud heuristic.  Use \"astar\" or \"euc\" or \"manhattan\"")
    end

    #run 1 to compile
    @load "Problems/$(prob)/5_1" lattice_inst
    tdp = @elapsed cost, path, gen = algof(lattice_inst, heur = heur)
    
    printstyled("\n Solving Lattice $(dims) Problems  || h(i): $(heur) || $(algo) \n ", color=:light_green)
    for n in Nvec
        printstyled("N = $(n) ||", color = :light_green)
        times_vec = Float64[]
        Zbreak_count = 0
        Threads.@threads for k = 1:10
            print(" $(k) ")
            @load "Problems/$(prob)/$(n)_$(k)" lattice_inst
            tdp = @elapsed cost, pathL, gen = algof(lattice_inst, heur = heur)
            # println(" $(tdp) ")
            @save "Solutions/$(prob)/$(n)_$(k)$(algo_tag)$(heur_tag)" tdp cost pathL gen
            push!(times_vec, tdp)
            cost == 0 && (Zbreak_count += 1)
            Zbreak_count > 3 && (nidx -= 1; break) 
        end
        if mean(times_vec) > tlim || Zbreak_count > 3
            printstyled("\n STOPPED EARLY || Lattice $(dims) || h(i): $(heur) || $(algo) \n", color=:light_red) 
            if Zbreak_count > 3
                nwholeidx = findall(x->x==n, Nvecwhole)[1] - 1
            elseif mean(times_vec) > tlim 
                nwholeidx = findall(x->x==n, Nvecwhole)[1] 
            end
            n = Nvecwhole[nwholeidx]
            @save "Solutions/END_$(prob)$(algo_tag)$(heur_tag)" n #save where we ended early.....
            return 0
        end
        meantts = round(mean(times_vec), digits = 2)
        println(" ||  mTTS (s) : $(meantts) ")
    end
    printstyled("\n SOLVED --- Lattice $(dims) Problems  || h(i): $(heur) || $(algo) --- \n", color=:light_green)
    n = Nvecwhole[end]
    @save "Solutions/END_$(prob)$(algo_tag)$(heur_tag)" n #save where we ended early.....

end 




