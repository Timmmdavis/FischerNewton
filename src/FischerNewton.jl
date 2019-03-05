module FischerNewton
using IterativeSolvers
using MAT
using LinearAlgebra
using SparseArrays
using Test
using Profile
using InteractiveUtils

#New Guys 07/02/2019
include("fischer_newton.jl")
include("CompResiduals.jl")
include("WorkOnJ.jl")
include("WorkOnJ_FastBigMats.jl")

end
