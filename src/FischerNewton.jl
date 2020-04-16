module FischerNewton
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using DelimitedFiles
#using KrylovKit
#using Test
using Krylov
#using MAT
#using Profile
#using InteractiveUtils

struct Floats{T<:Float64}
    alpha::T
    beta::T
    gamma::T
    rho::T
    max_iter::T
    tol_rel::T
    tol_abs::T
    lambda::T
end

mutable struct Integers{T<:Int64}
    iter::T
    N::T
    flag::T
    useSparse::T
end

mutable struct Vectors{T<:Vector}
    y::T
    phi::T
    phi_k::T
    phi_l::T
    phiM::T
    dx::T
    absdx::T
    y_k::T
    xdxtau::T
    x_k::T
    x::T
end


mutable struct Arrays{T<:Array}
    err::T
    nabdx::T
    test::T
    grad_f::T
    f_k::T
    old_err::T
end

mutable struct Matricies{T<:Array}
    J::T
    Jsubs::T
    JJ::T
    ISml::T
    JSml::T
    VSml::T
    phiT::T
    phi_kT::T
    nabla_phi::T
end

struct Bools{T<:Array}
    I::T
end

struct IntArrays{T<:Array}
    II::T
end




#New Guys 07/02/2019
include("fischer_newton.jl")
include("CompResiduals.jl")
include("WorkOnJ.jl")
include("WorkOnJ_Sparse.jl")
include("InitArrays.jl")
include("ResetInitArrays.jl")


end
