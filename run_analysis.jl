using Pkg; Pkg.activate(".")
using MixedRepSinglets
using HiRepParsing
using LatticeUtils
using DelimitedFiles
using HDF5
using Plots
using LaTeXStrings
using Statistics
using LinearAlgebra

include("scripts/utils.jl")
include("scripts/write_hdf5.jl")
ispath(hdf5out) || mkpath(hdf5out)
write_correlator && main_write_correlator_matrices(NsmearFUN,NsmearAS,hdf5parse,hdf5out,ensemble)