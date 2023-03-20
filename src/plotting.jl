## COLOR Guide

#magenta - generator on
#grey    - generator off
#lightcoral - noise restricted zone
# green     - fuel level 
# sienna    - battery level (SOC)

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
function plot_euc_graph_solution(euc_inst::EucGraph; label_strings::Vector{String}, label_units::Vector{String} = [""], label_vals::Matrix{Float64}, label_edge_idxs::Vector{Int}, path::Vector{Int64} = [], gen::Vector{Int64} = [], color_ends::Bool = true)
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
    avg_plt = plot(
        layer(x = Nvecs[1], y = avg_times[1], Geom.line, Geom.point, Theme(default_color = colors[1])),
        layer(x = Nvecs[2], y = avg_times[2], Geom.line, Geom.point, Theme(default_color = colors[2])),
        Scale.y_log10,
        Scale.x_continuous(format= :plain),
        Guide.xlabel(""),
        Guide.ylabel(" Mean Time to Solve (s)"),
        Theme(key_position = :right),
        Guide.manual_color_key("", [labels[1], labels[2]], [colors[1], colors[2]]))

    bplt1 = get_boxplot_plt(Nvecs[1], times[1], color = colors[1], xlabel = "") 
    bplt2 = get_boxplot_plt(Nvecs[2], times[2], color = colors[2], xlabel = "") 
    set_default_plot_size(plot_sizes[1], plot_sizes[2])
    stacked = vstack(avg_plt, bplt1, bplt2)
    return stacked
end

function get_boxplot_plt(Nvec::Vector{Int64}, times::Matrix{Float64}; color::String = "blue", xlabel::String="", xmax::Float64 = float(Nvec[end]))
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


function get_sol_vec(prob_type, prob_title; K = 10, conn = "_4conn", type = "euc", algo = "", prob = "DP", heur = "" )
    if prob_type == "euc"
        Nvec = [50:500:2000; 2000:1000:20000]
    elseif prob_type == "lattice"
        Nvec = 5:50
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
                return times[1:nidx-1], avg_times[1:nidx-1], Nvec[1:nidx-1]
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

