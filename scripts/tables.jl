using Pkg; Pkg.activate(".")
using DelimitedFiles

parameters = readdlm("input/parameters_gevp.csv",';';skipstart=1)
parameters_fitting = readdlm("input/parameters_corrfitter.csv",';';skipstart=1)
corrfitter_results = readdlm("output/corrfitter_results.csv",';';skipstart=0)

#check that the number of datasets match
@assert first(size(parameters)) == first(size(parameters_fitting)) == first(size(corrfitter_results)) 

function _get_measurement_id(ensemble,channel,table)
    isensemble = contains.(parameters[:,1],ensemble)
    ischannel  = contains.(parameters[:,2],channel)
    @assert sum(isensemble .* ischannel) == 1
    id = findfirst(isequal(1),isensemble .* ischannel)
end


ensembles = unique(parameters[:,1])
channels = unique(parameters[:,2])

for ensemble in ensembles
    
    id = _get_measurement_id(ensemble,"g5_singlet",parameters)

    nops = parameters[id,6]
    β, T, L, mf, mas = corrfitter_results[id,7], corrfitter_results[id,3], corrfitter_results[id,4], corrfitter_results[id,5], corrfitter_results[id,6]

    ma, Δma = corrfitter_results[id,8],  corrfitter_results[id,9] 
    mη, Δmη = corrfitter_results[id,10], corrfitter_results[id,11]

    idπF = _get_measurement_id(ensemble,"g5_nonsinglet_FUN",parameters)
    idρF = _get_measurement_id(ensemble,"g1_nonsinglet_FUN",parameters)
    idπA = _get_measurement_id(ensemble,"g5_nonsinglet_AS",parameters)
    idρA = _get_measurement_id(ensemble,"g1_nonsinglet_AS",parameters)
    mπF, ΔπF = corrfitter_results[idπF,8],  corrfitter_results[idπF,9] 
    mπA, ΔπA = corrfitter_results[idπA,8],  corrfitter_results[idπA,9] 
    mρF, ΔρF = corrfitter_results[idρF,8],  corrfitter_results[idρF,9] 
    mρA, ΔρA = corrfitter_results[idρA,8],  corrfitter_results[idρA,9] 

end