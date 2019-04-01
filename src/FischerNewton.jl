module FischerNewton
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using DelimitedFiles
#using KrylovKit
#using Test
#using Krylov
#using MAT
#using Profile
#using InteractiveUtils



#New Guys 07/02/2019
include("fischer_newton.jl")
include("CompResiduals.jl")
include("WorkOnJ.jl")
include("WorkOnJ_FastBigMats.jl")

end
