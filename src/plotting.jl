## COLOR Guide

#magenta - generator on
#grey    - generator off
#lightcoral - noise restricted zone
# green     - fuel level 
# sienna    - battery level (SOC)

function gen_plot(path, gen, euc_inst::EucGraphInt; max_time = 0, legendbool = true, uavi = 0)
    G = [!euc_inst.GFlipped[path[idx-1],path[idx]] for idx in 2:length(path)] #1 less than time_vec... 
    
    time = Vector(0:length(gen))
    timerep = repeat(time, inner = 2)
    popfirst!(timerep)
    pop!(timerep)
    genrep = repeat(gen, inner = 2)
    
    #normalize time vecs

    if max_time == 0
        time = time./time[end]
        timerep = timerep./timerep[end]
    else #need to normalize wrt to maximum time (0, maxhere/maxmax)
        time = time./max_time
        timerep = timerep./max_time
    end
    #else, we need to set xlims
    ylabel = "gen∈{0,1}"
    uavi != 0 && (ylabel = "UAV$(uavi)")
    if legendbool
        set_default_plot_size(11cm, 6cm)
        plt_gen = plot(
                    layer(x= timerep, y = genrep, Geom.line, Theme(default_color = "black")), 
                    layer(xmin=time[1:end-1], xmax = time[2:end].+ 0.01, Geom.vband, color = G), #max gen as ribbons....
                    Scale.color_discrete_manual("white", "lightcoral"),
                    Guide.colorkey(title="", labels = [""]),
                    # Guide.manual_color_key("",["Noise Resricted   ", "Gen Throttle"],["lightcoral", "black"]),
                    Guide.xlabel("Time (normalized)"),
                    Guide.ylabel(ylabel),
                    Guide.yticks(ticks = [0,1]),
                    Coord.cartesian(xmin = 0, xmax = 1),
                    # Theme(key_position=:top)
        )
    elseif !legendbool #if not bottom, remove legend and xlabel
        set_default_plot_size(20cm, 5cm)
        plt_gen = plot(
                    layer(x= timerep, y = genrep, Geom.line, Theme(default_color = "black")), 
                    layer(xmin=time[1:end-1], xmax = time[2:end].+ 0, Geom.vband, color = G), #max gen as ribbons....
                    Scale.color_discrete_manual("white", "lightcoral"),
                    Guide.colorkey(title="", labels = [""]),
                    # Guide.manual_color_key("",["Noise Resricted   ", "Gen Throttle"],["lightcoral", "black"]),
                    Guide.xlabel(""),
                    Guide.xticks(ticks = []),
                    Guide.ylabel(ylabel),
                    Guide.yticks(ticks = [0,1]),
                    Coord.cartesian(xmin = 0, xmax = 1),
                    # Theme(key_position=:top)
        )
    end
    return plt_gen

end
   


function plot_euc_graph(euc_inst; path = [], gen = [], color_ends = true)
    #make Graph() then just graph plot
    set_default_plot_size(20cm, 20cm)
    g = make_graph(euc_inst, false)
    N = length(euc_inst.Alist)
    edge_index(e::Edge) = edgemap[e]
    nE = ne(g)
    edge_colors = [colorant"gray" for i in 1:nE]
    node_colors = [colorant"darkgray", colorant"cyan", colorant"magenta", colorant"lightcoral", colorant"magenta"]
    node_labels = ones(Int, N)
    nodesize = 0.0008*ones(N)
    edgelinewidth = 0.125*ones(nE)
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
            if euc_inst.GFlipped[src,dst] == true
                edge_colors[i] = node_colors[4]
            elseif euc_inst.GFlipped[src,dst] == false
                edge_colors[i] = node_colors[1]
            # elseif Gsummed[src] >= 1 || Gsummed[dst] >= 1
            #     println("here")
            #     edge_colors[i] = colorant"red"
            #     if Gsummed[src] >= 1
            #         node_labels[src] = 4
            #     end
            #     if Gsummed[dst] >= 1
            #         node_labels[dst] = 4
            #     end
            end

        end
        for k in 2:length(path)
            i = path[k-1]
            j = path[k]
            edge_idx = edgemap[Edge(i,j)]
            edgelinewidth[edge_idx] = 0.5
            # nodesize[i] = 0.011
            # nodesize[j] = 0.011
            if gen[k-1] == 0
                edge_colors[edge_idx] = node_colors[2]
                # k == 2 && (node_labels[i] = 2)
                # node_labels[j] = 2
            elseif gen[k-1] == 1
                edge_colors[edge_idx] = node_colors[3]
                k == 2 && (node_labels[i] = 3)
                # node_labels[j] = 3
            end
        end
        
    else
        edgelist = collect(edges(g))
        edgemap = Dict{Edge, Int}()
        Gsummed = sum(euc_inst.GFlipped, dims = 1)
        edgelinewidth = 0.25*ones(nE)
        for (i,e) in enumerate(edgelist)
            edgemap[e] = i
            edgemap[reverse(e)] = i
            src, dst = e.src, e.dst
            if euc_inst.GFlipped[src,dst] == true
                edge_colors[i] = node_colors[4]
            elseif euc_inst.GFlipped[src,dst] == false
                edge_colors[i] = node_colors[1]
            elseif Gsummed[src] >= 1 || Gsummed[dst] >= 1
                println("here")
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
    
    plt = gplot(g, euc_inst.locs[:,1], euc_inst.locs[:,2], edgestrokec = edge_colors, nodefillc = nodefillc2, nodesize=nodesize, edgelinewidth = edgelinewidth, NODESIZE = maximum(nodesize), EDGELINEWIDTH = maximum(edgelinewidth))
    return plt
end

function plot_MAPF(euc_inst, paths; color_ends = true, NR = false)
    #for each path/gen, plot the solution....
    #make Graph() then just graph plot
    set_default_plot_size(20cm, 20cm)
    g = make_graph(euc_inst, false)
    N = length(euc_inst.Alist)
    edge_index(e::Edge) = edgemap[e]
    nE = ne(g)
    edge_colors = [colorant"gray" for i in 1:nE]
    node_colors = [colorant"darkgray", colorant"cyan", colorant"magenta", colorant"red", colorant"magenta", colorant"black"]
    path_colors = [colorant"cyan", colorant"green", colorant"magenta", colorant"orange", colorant"chartreuse1", colorant"sienna"]
    node_labels = ones(Int, N)
    nodesize = 0.0008*ones(N)
    edgelinewidth = 0.125*ones(nE)
    agent_colors = fill(colorant"black", length(paths))
    
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
        if euc_inst.GFlipped[src,dst] == true && NR == true
            edge_colors[i] = node_colors[4]
        elseif euc_inst.GFlipped[src,dst] == false
            edge_colors[i] = node_colors[1]
        end

    end
    for agent_idx in 1:length(paths)
        path = paths[agent_idx]
        gen = zeros(length(path)-1)
        if agent_idx <= 6
            agent_color = path_colors[agent_idx]
        else 
            #pick random color
            agent_color = rand(path_colors)
        end
        for k in 2:length(path)
            i = path[k-1]
            j = path[k]
            edge_idx = edgemap[Edge(i,j)]
            edgelinewidth[edge_idx] = 0.5
            edge_colors[edge_idx] = agent_color
            if color_ends
                nodesize[path[1]] = 0.033
                nodesize[path[end]] = 0.033
                node_labels[path[1]] = 6
                node_labels[path[end]] = 4
            end            
        end
    end
    nodefillc2 = node_colors[node_labels]
    plt = gplot(g, euc_inst.locs[:,1], euc_inst.locs[:,2], edgestrokec = edge_colors, nodefillc = nodefillc2, nodesize=nodesize, edgelinewidth = edgelinewidth, NODESIZE = maximum(nodesize), EDGELINEWIDTH = maximum(edgelinewidth), )
    return plt
end


function plot_euc_graph_solution(euc_inst::EucGraph; label_strings::Vector{String}, label_units::Vector{String} = [""], label_vals::Matrix{Float64}, label_edge_idxs::Vector{Int}, path::Vector{Int64} = [], gen::Vector{Int64} = [], color_ends::Bool = true)
    #make Graph() then just graph plot
    set_default_plot_size(20cm, 20cm)
    g = make_graph(euc_inst, false)
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
        edgemap = Dict{Edge, Int}()g
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

function get_stacked_plots(Nvecs, times, avg_times; labels = ["NODE", "LABEL"], colors = ["red", "grey"], plot_sizes = [20cm, 30cm])
    xmax = maximum([Nvecs[1][end], Nvecs[2][end]])
    avg_plt = plot(
        layer(x = Nvecs[1], y = avg_times[1], Geom.line, Geom.point, Theme(default_color = colors[1])),
        layer(x = Nvecs[2], y = avg_times[2], Geom.line, Geom.point, Theme(default_color = colors[2])),
        Scale.y_log10,
        Scale.x_continuous(format= :plain, maxvalue = xmax),
        Guide.xlabel(""),
        Guide.ylabel(" Mean Time to Solve (s)"),
        Theme(key_position = :top),
        Guide.manual_color_key("", [labels[1], labels[2]], [colors[1], colors[2]]))

    bplt1 = get_boxplot_plt(Nvecs[1], times[1], color = colors[1], xlabel = "", xmax = xmax) 
    bplt2 = get_boxplot_plt(Nvecs[2], times[2], color = colors[2], xlabel = "Number of Nodes", xmax = xmax) 
    set_default_plot_size(plot_sizes[1], plot_sizes[2])
    stacked = vstack(avg_plt, bplt1, bplt2)
    return stacked
end

function get_boxplot_plt(Nvec::Vector{Int64}, times::Matrix{Float64}; color::String = "blue", xlabel::String="", xmax::Int64 = Nvec[end])
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
    g = make_graph(euc_inst, false)
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

function get_sol_vec(prob_type, prob_title; K = 10, conn = "_4conn", type = "euc", algo = "", prob = "DP", heur = "" )
    #changing this from orig so we pull the END problem size, don't have to worry about loading old solutions from prior runs that went farther accidently....
    nEND = Int64
    try #stinky hack for not having saved nENDS for runs to the complete end....
        # if prob_title == "euc_probs_2D"
        #     # println("Solutions\\END_euc_probs2D_$(algo)$(heur)")
        #     @load "Solutions\\END_euc_probs2D_$(algo)$(heur)" n
        #     nEND = n + 0    
        # else
        @load "Solutions\\END_$(prob_title)$(algo)$(heur)" n
        nEND = n + 0
        # end
    catch #stinky hack... catching this and just assuming we solved to the end....
        if prob_type == "euc"
            nEND = 20000
        elseif prob_type == "lattice"
            nEND = 50
        end
    end
    if prob_type == "euc"
        if nEND > 2000
            Nvec = [50:500:2000; 2000:1000:nEND]
            if prob_title == "connectivity_expr"
                Nvec = [50:1000:2000; 2000:2000:20000]
            end
        else
            Nvec = Vector(50:500:nEND)
            if prob_title == "connectivity_expr"
                Nvec = Vector(50:1000:nEND)
            end
        end
    elseif prob_type == "lattice"
        Nvec = Vector(5:nEND)
    end

    times = zeros(length(Nvec),K)
    avg_times = zeros(length(Nvec))
    
    nidx = 0
    for n in Nvec
        nidx += 1
        for k = 1:K
            try
                if prob == "MILP"
                    # println("Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)")
                    @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)" tMILP
                    time_i = tMILP
                elseif prob_title == "connectivity_expr"
                    @load "Solutions/$(prob_title)/$(n)_$(k)$(conn)$(algo)$(heur)" tdp
                    time_i = tdp
                else
                    @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)$(heur)" tdp
                    time_i = tdp
                end
                times[nidx,k] = time_i
            catch #if here, then we are at the end of saved prolems... return up to the prior Nidx....
                return times[1:nidx-1, :], avg_times[1:nidx-1], Nvec[1:nidx-1]
            end
        end
        avg_times[nidx] = mean(times[nidx,:])
    end

    return times, avg_times, Nvec
end

function get_sol_vec_conn_expr(conn::Int64, algo::String = "_label")
#write code to load solutions iterating Nvec.  If doesn't exist, assume we didnt solve
    Nvec = [50:1000:2000; 2000:2000:20000]
    times = zeros(length(Nvec), 10)
    avg_times = zeros(length(Nvec))

    nidx = 1
    for n in Nvec
        for k = 1:10
            try
                @load "Solutions/connectivity_expr/$(n)_$(k)_$(conn)conn$(algo)" tdp;
                time_i = tdp
                times[nidx,k] = time_i
            catch #if here, then we are at the end of saved prolems... return up to the prior Nidx....
                return times[1:nidx-1, :], avg_times[1:nidx-1], Nvec[1:nidx-1]
            end
        end
        avg_times[nidx] = mean(times[nidx,:])
        nidx += 1
    end


    return times, avg_times, Nvec
end

function get_sol_vec_old(prob_type, prob_title; K = 10, conn = "_4conn", type = "euc", algo = "", prob = "DP", heur = "" )
    if prob_type == "euc"
        Nvec = [50:500:2000; 2000:1000:20000]
    elseif prob_type == "lattice"
        Nvec = Vector(5:50)
    end
    times = zeros(length(Nvec),K)
    avg_times = zeros(length(Nvec))

    nidx = 0
    for n in Nvec
        nidx += 1
        for k = 1:K
            try
                if prob == "MILP"
                    @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)" tMILP
                    time_i = tMILP
                else
                    @load "Solutions/$(prob_title)/$(n)$(conn)_$(k)$(algo)$(heur)" tdp
                    time_i = tdp
                end
                times[nidx,k] = time_i
            catch #if here, then we are at the end of saved prolems... return up to the prior Nidx....
                return times[1:nidx-1, :], avg_times[1:nidx-1], Nvec[1:nidx-1]
            end
        end
        avg_times[nidx] = mean(times[nidx,:])
    end

    return times, avg_times, Nvec
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

