module TestPkg


using LinearAlgebra
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
# using CPLEX
using Ipopt 
using DifferentialEquations
using Interpolations        


include("struct_defs.jl")
include("label_label_sel.jl")
include("label_node_sel.jl")
include("label_utills.jl")
include("plotting.jl")
include("battery.jl")
include("opt_contrl.jl")


precompile(hybrid_label_selection)
precompile(hybrid_node_selection)
precompile(MILP_to_opt_ctrl)
precompile(plot_euc_graph)
precompile(plot_hybrid_soln)
# include("MILP_definition3D.jl")



export hybrid_label_selection
export hybrid_node_selection
export MILP_hybrid

export EucGraphInt

export get_path
export get_gen

export plot_euc_graph
export plot_euc_graph_solution

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


export OptControlProb
export OptControlSolution
export MILP_to_opt_ctrl



#-------END OF MODULE---------------
end
