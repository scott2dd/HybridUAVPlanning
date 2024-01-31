###########################################
## 1: Data Structs


abstract type Label end

@kwdef struct MyLabel <: Label
    gcost::Float64
    fcost::Float64
    hcost::Float64
    node_idx::Int64
    prior_node_idx::Int64
    _hold_came_from_prior::Int64
    came_from_idx::Int64
    pathlength::Int64
    _hold_gen_track_prior::Int64
    gentrack_idx::Int64
    gen_bool::Int64
    batt_state::Float64
    gen_state::Float64
end

Base.isless(a::Label, b::Label) = (a.fcost, a.gcost) < (b.fcost, a.gcost) #tie breaker is gcost.
#note, for lattice2D, gcost seems to be best tie breaker.
#TODO did not test for lattice 3D or eucs.

struct ProblemDef
    S::Int64
    E::Int64
    Alist::Vector{Vector{Any}}
    A::Array{Float64}
    C::Array{Float64}
    G::Array{Float64}
    Z::Array{Float64}
    B0::Float64
    Q0::Float64
    anchor_list::Array{Float64}
    Dim::Vector{Float64}
end

struct MapDef
    Alist::Vector{Vector{Any}}
    A::Array{Float64}
    C::Array{Float64}
    G::Array{Float64}
    Z::Array{Float64}
    anchor_list::Vector{Float64}
    Dim::Vector{Float64} #[n.m]
    nN::Int64 #number of houses (total N - n*m + nN)
    GPS_locs::Array{Float64}
    void_list::Vector{Float64}
end

struct MapProb
    S::Int64
    E::Int64
    B0::Int64
    Q0::Int64
    Bmax::Float64
end



# 3D Structs
#Define label:
#(C, B, G, i, k, F, GenOnPrior)
#In heap:
#(F, [C, B, G, i, k, F, GenOnPrior])

struct GridProb
    S::Int64
    E::Int64
    GFlipped::SparseMatrixCSC{Bool, Int64}
    B0::Float64
    Q0::Float64
    Bmax::Float64
    StartCost::Float64
    anchor_list::SparseVector{Float64,Int64}
    Dim::Vector{Float64}
end

struct GraphDef3D
    Alist::Vector{Vector{Any}}
    A::SparseMatrixCSC{Float64, Int64}
    F::SparseMatrixCSC{Float64, Int64}
    C::SparseMatrixCSC{Float64, Int64}
    Z::SparseMatrixCSC{Float64, Int64}
end
struct FullDef3D
    S::Int64
    E::Int64
    Alist::Vector{Vector{Any}}
    A::Array{Float64}
    F::Array{Float64}
    C::Array{Float64}
    G::Array{Bool}
    Z::Array{Float64}
    B0::Float64
    Q0::Float64
    Bmax::Float64
    StartCost::Float64
    anchor_list::Array{Float64}
    Dim::Vector{Float64}
end
struct MapDef3D
    Alist::Vector{Vector{Any}}
    A::Array{Float64}
    C::Array{Float64}
    F::Array{Float64}
    G::Array{Float64}
    Z::Array{Float64}
    anchor_list::Vector{Float64}
    Dim::Vector{Float64} #[n.m]
    nN::Int64 #number of "houses" (total N - n*m + nN)
    GPS_locs::Array{Float64}
    void_list::Vector{Float64}
end
struct EucGraph
    S::Int64
    E::Int64
    Alist::Vector{Vector{Any}}
    F::SparseMatrixCSC{Bool, Int64}
    C::SparseMatrixCSC{Float64, Int64}
    GFlipped::SparseMatrixCSC{Bool, Int64}
    Z::SparseMatrixCSC{Float64, Int64}
    B0::Float64
    Q0::Float64
    Bmax::Float64
    StartCost::Float64
    anchor_list::SparseVector{Int64}
    locs::Matrix{Float64}
end
struct EucGraphInt
    S::Int64
    E::Int64
    Alist::Vector{Vector{Any}}
    F::SparseMatrixCSC{Bool, Int64}
    C::SparseMatrixCSC{Int64}
    GFlipped::SparseMatrixCSC{Bool, Int64}
    Z::SparseMatrixCSC{Int64}
    B0::Int64
    Q0::Int64
    Bmax::Int64
    StartCost::Int64
    anchor_list::SparseVector{Int64}
    locs::Matrix{Float64} #leave this as float... need to deal with herustic cost.... can just round down?
end
struct Battery{N<:AbstractFloat, F<:Function}
    R::N
    Cmax::N
    C::N
    OCV::F
end
