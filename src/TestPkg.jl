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
# using CPLEX  -> update each system w/ CPLEX 22.11 binaries installed....
using Ipopt 
using DifferentialEquations
using Interpolations        

#=
workflow for this package:
1) test new functions and data types in main project folder
2) once tested, add here 
3) "up TestPkg" in main project folder

This is a julia package for "permanent" or well-tested and often-used functions in my project
Exports as-needed (many old functions and structs in here that I don't use anymore, may use in future then export when needed)


NOTE: CPLEX not installed in this package.  Dependencies would not work with old CPLEX.jl needed for IBMCPLEXv12.9
Need to update machine CPLEX binaries, update julia CPLEX_BINARY_PATH, then add new CPLEX.jl (unpinned) and fix code to use new CPLEX 
=#

include("struct_defs.jl")
include("label_utills.jl")
include("label_label_sel.jl")
include("label_node_sel.jl")
include("plotting.jl")
include("battery.jl")
include("opt_contrl.jl")
include("solve_functions.jl")
# include("MILP_definition3D.jl")  #later will update with correct CPLEX



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
export get_tag
export solve_gen_optimal_control


export solve_euc
export solve_lattice
#do precompiles later? This gives an error when doing auto-Julia-Pkg.jl precompilation.....
# precompile(hybrid_label_selection)
# precompile(hybrid_node_selection)
# precompile(MILP_to_opt_ctrl)
# precompile(plot_euc_graph)
# precompile(plot_hybrid_soln)
#-------END OF MODULE---------------
end
