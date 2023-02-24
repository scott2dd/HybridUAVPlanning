
## OLD LABEL_SEL_INT FUNCTION!!!!! 
#this is a copy prior to changing to recursive data structs for path and generator
# Function for EucFloat type
# function hybrid_with_LB_OLD(def::EucGraph; heur = "astar")
#     S, E, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
#     Bstart, Qstart = def.B0, def.Q0
#     Cmin = minimum(nonzeros(C))
#     # Bstart = mean(nonzeros(C)) * 4
#     genpen = Cmin*.0
#     Bstart = 4*mean(nonzeros(C))
#     Qstart = 9999
#     Bmax = Bstart
#     SC = def.StartCost
#     SC = 0
#     anchor = def.anchor_list  #if a grid problem, this is zeros(N)
#     locs = def.locs
#     # G = .!def.GFlipped
#     Dim = def.locs           ####### pass third argument of heuristics as locs
   
#     N = length(Alist)
#     graph = make_graph(def)
#     Fvec = fill(NaN, N)
#     heur_astar = get_heur_astar(E, locs, Fvec)
#     heur_astar(S)
#     if heur == "astar"
#         heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
#     elseif heur == "euc"
#         heur_label! = get_heur_label_euc(Fvec, locs, E)
#     elseif heur == "manhattan"

#     else
#         println("invalid heuristic... breaking")
#         return 0,0,0,0,0,0,0
#     end
#     if heur_label!(S) == Inf     
#         println("Infinite Start Heuristic...")
#         return 0,0,0,0,0,0,0
#     end
    
#     Q = MutableBinaryMinHeap([   (heur_label!(S) + 0, [0.0, Bstart, Qstart, S, 1.0, heur_label!(S) + 0]) ] )
#     P = [zeros(0,6) for k = 1:size(C,1)] #set of treated labels | new labels compare to this AND Q
#     X = Vector{Int}[ [S] ]  #hold path info
#     Y = Vector{Int}[  []   ] #hold gen info
    
#     z=0
#     look_stack = Int[]
#     dist_stack = Float32[]
#     f_stack = Float32[]
#     prog = ProgressUnknown("Working hard...", spinner=true)
#     while true #loop until get to end node, or Q is empty
#         isempty(Q) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break)

#         #pull minimum cost label....
#         next = pop!(Q)
#         label_treated = next[2]
#         i = Int(label_treated[4])
#         K = Int(label_treated[5])
#         GenPrior = 1 #label_treated[7]
#         append!(look_stack, i)
#         append!(dist_stack, heur_label!(i))
#         append!(f_stack, next[1])
#         # println(f_stack[end])
#         #now add it to P
#         P[i] = vcat(P[i], reshape(label_treated,(1,6)))
#         # println(i)
#         #if we are at the end, then stop algo and return this label
#         if i == E
#             opt_cost =  label_treated[1]
#             opt_path =  X[K]
#             opt_gen =   Y[K]
#             return opt_cost, opt_path, opt_gen, length(X), look_stack, dist_stack
#         end
#         for j in Alist[i]
#             j==i && continue
#             j∈X[K] && continue
#             gc = label_treated[1] + C[i,j]
#             h =  heur_label!(j) 
#             Fbin = F[i,j] # if we can glide, then ALWAYS glide
#             label = []
#             # No Gen 
#             # if label_treated[2] + Z[i,j] < Bmax && def.GFlipped[i,j] == 0
#                         # Gen On   IF allowed...
#                 if def.GFlipped[i,j] == 0 && label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior) >= 0 && label_treated[3]-Z[i,j] >= 0 && label_treated[2]-C[i,j] + Z[i,j] <= Bmax
#                     label =  [gc + genpen, label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior), label_treated[3]-Z[i,j],  j, Int(length(X)+1), gc+h]
#                     # println("label, $(i)->$(j) GEN")

#                     if EFF_heap(Q, label) && EFF_P(P, label)
#                         push!(Q, (gc+h + genpen,label))

#                         push!(X, [X[K]; j ])
#                         push!(Y, [Y[K]; 1 ])
#                     end
#                 end
#             # else
#                 if label_treated[2]-C[i,j]*(1-Fbin) >= 0
#                     label =  [gc, label_treated[2]-C[i,j]*(1-Fbin), label_treated[3],  j, Int(length(X)+1), gc+h]
#                     # println("label, $(i)->$(j)")
                    
#                     if EFF_heap(Q, label) && EFF_P(P, label)
#                         push!(Q, (gc+h,label))
#                         push!(X, [X[K]; j ])
#                         push!(Y, [Y[K]; 0 ])
#                     end
#                 end
#             # end
#         end
#         z+=1
#         z == 200_000 && (printstyled("Z BREAK... @Z=$(z)\n", color=:light_red); break)
#         z%500 == 0 && ProgressMeter.next!(prog)

#     end
#     return 0,0,0,0,[0], [0]
# end
# function hybrid_with_LB_OLD(def::EucGraphInt; heur = "astar") 
#     S, E, Alist, F, C, Z = def.S, def.E, def.Alist, def.F, def.C, def.Z
#     Bstart, Qstart = def.B0, def.Q0
#     Cmin = minimum(nonzeros(C))
#     # Bstart = mean(nonzeros(C)) * 4
#     genpen = Cmin *0
#     Bstart = Int(floor(4*mean(nonzeros(C))))
#     Qstart = 9999
#     Bmax = Bstart
#     SC = def.StartCost
#     SC = 0
#     anchor = def.anchor_list  #if a grid problem, this is zeros(N)
#     locs = def.locs
#     Dim = def.locs 
   
#     N = length(Alist)
#     graph = make_graph(def)
#     Fvec = fill(2^63 - 1, N)
#     heur_astar = get_heur_astar(E, locs, Fvec)
#     heur_astar(S)
#     if heur == "astar"
#         heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
#     elseif heur == "euc"
#         heur_label! = get_heur_label_euc(Fvec, locs, E)
#     elseif heur == "manhattan"

#     else
#         println("invalid heuristic... breaking")
#         return 0,0,0,0,0,0,0
#     end
#     if heur_label!(S) == Inf     
#         println("Infinite Start Heuristic...")
#         return 0,0,0,0,0,0,0
#     end

#     Q = MutableBinaryMinHeap([   (heur_label!(S) + 0, [0, Bstart, Qstart, S, 1, heur_label!(S) + 0]) ] )
#     P = [zeros(Int, 0,6) for k = 1:size(C,1)] #set of treated labels | new labels compare to this AND Q
#     X = Vector{Int}[ [S] ]  #hold path info
#     Y = Vector{Int}[  []   ] #hold gen info
    
#     z=0
#     prog = ProgressUnknown("Working hard...", spinner=true)
#     while true #loop until get to end node, or Q is empty
#         isempty(Q) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break)

#         #pull minimum cost label....
#         next = pop!(Q)
#         label_treated = next[2]
#         i = label_treated[4]
#         K = label_treated[5]
#         GenPrior = 1 #label_treated[7]
#         P[i] = vcat(P[i], reshape(label_treated,(1,6)))

#         if i == E
#             opt_cost =  label_treated[1]
#             opt_path =  X[K]
#             opt_gen =   Y[K]
#             return opt_cost, opt_path, opt_gen
#         end
        
#         for j in Alist[i]
#             j==i && continue
#             j∈X[K] && continue
#             gc = label_treated[1] + C[i,j]
#             h =  heur_label!(j) 
#             Fbin = F[i,j] # if we can glide, then ALWAYS glide
#             label = []
#             # No Gen 
#             if def.GFlipped[i,j] == 0 && label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior) >= 0 && label_treated[3]-Z[i,j] >= 0 && label_treated[2]-C[i,j] + Z[i,j] <= Bmax
#                 label =  [gc + genpen, label_treated[2]-C[i,j]*(1-Fbin) + Z[i,j] - SC*(1-GenPrior), label_treated[3]-Z[i,j],  j, Int(length(X)+1), gc+h]
#                 if EFF_heap(Q, label) && EFF_P(P, label)
#                     push!(Q, (gc+h + genpen,label))

#                     push!(X, [X[K]; j ])
#                     push!(Y, [Y[K]; 1 ])
#                 end
#             end
#             if label_treated[2]-C[i,j]*(1-Fbin) >= 0
#                 label =  [gc, label_treated[2]-C[i,j]*(1-Fbin), label_treated[3],  j, Int(length(X)+1), gc+h]
                
#                 if EFF_heap(Q, label) && EFF_P(P, label)
#                     push!(Q, (gc+h,label))
#                     push!(X, [X[K]; j ])
#                     push!(Y, [Y[K]; 0 ])
#                 end
#             end
#         end
#         z+=1
#         z == 200_000 && (printstyled("Z BREAK... @Z=$(z)\n", color=:light_red); break)
#         z%500 == 0 && ProgressMeter.next!(prog)
#     end
#     return 0,0,0
# end


# statements_files=wef er er er erg erg er "/.vscode/trace_startup.jl",]
# execution_files = ["/.vscode/precompile_execution_file.jl",]


## Test for curve fitting a "square wave" for binary constraints on quiet zones
# function get_bin_data()
#     len = 10
#     x = zeros(len)
#     y = zeros(len)
#     switch01 = true
#     for i in 1:len
#         x[i] = i-1
#         if switch01
#             y[i] = 0
#             switch01= false
#         else
#             y[i] = 1
#             switch01 = true
#         end
#     end

#     #now fill in points between.....

#     xx = LinRange(x[1], x[end], 1000)
#     yy = zeros(length(xx))
#     R = x[1]
#     for i in 2:length(y)
#         yi = y[i-1]
#         L,R = R, x[i]
#         yy[L .< xx .< R] .= yi 
#     end
#     return xx,yy
# end
# x,y = get_bin_data()
# using Gadfly, CurveFit

# fit = curve_fit(LinearFit, x,y)
# @time (h,f,a,b) = fourierSeriesSampledReal(x, y);
# @time t,u = fourierSeriesSynthesisReal(f,a,b);
# plt = plot(layer(x = x, y = y, Geom.line, Theme(default_color="grey")),
#             layer(x = t, y = u, Geom.line), Theme(default_color="red"))


###################################################
## 7: Heap Testing...
# X = rand(Int.(1:10), 200,3)
# X = array_to_list(X)

# # if we use Mutable heap, can we remove el's?
# heap_mut = MutableBinaryMinHeap{}(X)
# new_list = array_to_list(rand(Int.(1:10), 100,3))

# @time eff_heap = merge_sort!(heap_mut, new_list)
# @time sorted_list = merge_3D(X, new_list)



#####################################################################
## 8: Junk Code from label_node_sel.jl

# function merge_brute(Qi::BinaryMinHeap{Tuple{Float64, Vector{Float64}}}, new_labels::BinaryMinHeap{Tuple{Float64, Vector{Float64}}})
#     isempty(Qi) && return new_labels
#     M = deepcopy(Qi)
#     map_labels, map_Q = new_labels.node_map, Qi.node_map
#     #loop through each element of Q and compare it to each element in new_labels
#     for iL in 1:length(map_labels)
#         map_labels[iL] == 0 && continue
#         new_labels_i = new_labels[iL]
#         label = new_labels_i[2]
#         keep = true
#         for iQ in 1:length(map_Q)
#             map_Q[iQ] == 0 && continue
#             q = Qi[iQ][2]
#             if dom(q, label)
#                 keep = false
#                 break
#             end
#         end
#         keep && push!(M, new_labels_i)
#     end
#     return M
# end

# function merge_sort( Qi::Vector{Vector{Int64}}, new_labels::Vector{Vector{Int64}} ) 
    
#     Mtemp = vcat(Qi, new_labels)
#     sort!(Mtemp)
#     boolv = ones(Bool, length(Mtemp))
#     #now go through each element
#     for i in eachindex(Mtemp)
#         label =  Mtemp[i]
#         for ii in 1:(i-1)
#             if dom_min(Mtemp[ii], label)
#                 boolv[i] = 0
#                 break
#             end
#         end
#     end
#     return Mtemp[boolv]
# end
# function merge_sort!(Qi::MutableBinaryMinHeap{Vector{Int64}}, new_labels::Vector{Vector{Int64}}) #can delete ineff elements.... but keep fast sorting with heap adds and pops

#     #add each new_label into Qi
#     [push!(Qi, i) for i in new_labels]

#     ineff_list = Int64[]

#     map_Qi = Qi.node_map
#     sorted = sortslices(hcat(map_Qi, 1:length(map_Qi)),dims=1)
#     for i in 1:length(map_Qi) 
#         idx = sorted[i,2]
#         map_Qi[idx] == 0 && continue
#         label = Qi[idx]
#         for ii in 1:(i-1)
#             idx_ii = sorted[ii,2]
#             map_Qi[idx_ii] == 0 && continue
#             if dom_min(Qi[idx_ii], label)
#                 push!(ineff_list, idx)
#                 break
#             end
#         end
#     end
#     for i in ineff_list 
#         delete!(Qi, i)
#     end
# end


# function EFF_X(X::Vector{Vector{Float64}})
#     Xout = Vector{Float64}[]
#     eq_remove = []
#     eq_keep = []
#     for i in 1:length(X)
#         x = X[i]
#         boolx = true
#         eq_before = false
#         for j in 1:length(X)
#             xx = X[j] 
#             i == j && continue
#             dom(xx, x) && (boolx = false)
#             if x == xx #if they are equal label values, may need to turn boolx back on
#                 if i ∉ eq_keep && i ∉ eq_remove #
#                     boolx = true
#                     push!(eq_keep, i)
#                     push!(eq_remove, j)
#                 elseif i ∈ eq_keep
#                     boolx = true
#                     push!(eq_remove, j)
#                 end
#                 # eq_before = true
#             end
#         end
#         boolx && push!(Xout, x)
#     end
#     return Xout
# end



## More Junk Code....    
## Plotting Testingusing GLMakie, SGtSNEpi, SNAPDatasets
# using GLMakie, SGtSNEpi, SNAPDatasets
# GLMakie.activate!()

# g = loadsnap(:as_caida)
# y = sgtsnepi(g);
# show_embedding(y;
#   A = adjacency_matrix(g),        # show edges on embedding
#   mrk_size = 1,                   # control node sizes
#   lwd_in = 0.01, lwd_out = 0.001, # control edge widths
#   edge_alpha = 0.03 )             # control edge transparency




# using Graphs
# using GraphPlot
# using Gadfly
# locs = rand(2,10)
# g, dists = euclidean_graph(locs)
# gplot(g, xlocs = locs[:,1], ylocs = locs[:,2])


# ## Print "animation" testing...
# z = 0
# count = 0
# strings = ["O o o","o O o","o o O"]
# println("\n\n\n")
# while true
#     z += 1
#     sleep(.05)
#     if z%10 == 0
#         count += 1
#         print("\r")
#         print(strings[count])
#         count == 3 && (count = 0)
#     end
#     z > 100 && (print("\rFinished"); break)
# end

# ## progress meter
# using ProgressMeter
# prog = ProgressUnknown("  Working hard:", spinner=true)

# while true
#     ProgressMeter.next!(prog)
#     rand(1:2*10^8) == 1 && break
#     sleep(1)
# end
# ProgressMeter.finish!(prog)

# Old Code...
# Eij = 1000
# SOCmax = 10000
# xy = [0.2 3.25; 0.5 3.3; 0.9 3.38]
# xyinv = [0.2 1/3.25; 0.5 1/3.3; 0.9 1/3.38]
# SOC = 25:0.5:90
# # fV(x) =0.0015116*x + 3.21686 #from regression
# fV(x) = 0.0018649*x + 3.21054
# onebyFV(x) = -0.00016961081*x + .31124 #from regression

# deltaSOCf(x) = Eij/(fV(x)*SOCmax)*100
# deltaSOClinearf(x) = Eij/(SOCmax)*onebyFV(x)*100
# Vhat = fV.(SOC)
# onybyVhat = onebyFV.(SOC)

# x = -100:.1:1000
# deltaSOCx = deltaSOCf.(x)
# deltaSOCx_linear = deltaSOClinearf.(x)
# deltaSOC  = deltaSOCf.(SOC)
# deltaSOC_linear = deltaSOClinearf.(SOC)


# plt_deltaSOC_full = plot(
# layer(x=x, y = deltaSOCx , Geom.line, Theme(default_color="red")),
# layer(x=x, y = deltaSOCx_linear, Geom.line, Theme(default_color="grey")),
# Guide.xlabel("SOC"), Guide.ylabel("ΔSOC for Fixed Flight Leg ") )





# #plots...
# plt_V = plot(x = SOC, y = Vhat , Geom.line, Scale.x_continuous(minvalue=0, maxvalue=100) ,Scale.y_continuous(minvalue=0, maxvalue=3.5), Guide.xlabel("SOC"), Guide.ylabel("Open Circuit Voltage") )
# plt_compare = plot(
# layer(x=SOC, y = deltaSOC , Geom.line, Theme(default_color="red")),
# layer(x=SOC, y = deltaSOC_linear, Geom.line, Theme(default_color="grey")),
# Scale.x_continuous(minvalue=0, maxvalue=100), Guide.xlabel("SOC"), Guide.ylabel("ΔSOC for Fixed Flight Leg "), Guide.manual_color_key("Legend", ["Actual", "Linear Fit"], ["red", "grey"]))


# set_default_plot_size(20cm, 15cm)
# horizontal = hstack(compose(context(0,0,1/2, 1.0), render(plt_V)),
#                    compose(context(0,0,1/2, 1.0), render(plt_compare)) )
# plt_stack = vstack(compose(context(0,0,1.0, 0.5), horizontal),
#                    compose(context(0,0,1.0, 0.5), render(plt_deltaSOC_full)) )
# display(plt_stack)
# --sysimage=C:\Users\Drew\.julia\config\sys_atom.so

#bonk