function write_all_tables(Nsmear,paramterfile_gevp,paramterfile_corrfitter,corrfitterpath,tablepath)

    parameters = readdlm(paramterfile_gevp,';';skipstart=1)
    parameters_fitting = readdlm(paramterfile_corrfitter,';';skipstart=1)
    corrfitter_results = readdlm(joinpath(corrfitterpath,"corrfitter_results.csv"),';';skipstart=0)
    corrfitter_results_HR = readdlm(joinpath(corrfitterpath,"corrfitter_results_HR.csv"),';';skipstart=0)

    ispath(tablepath) || mkpath(tablepath)

    io_resultsMR = open(joinpath(tablepath,"table_results_MR.csv"),"w")
    io_results = open(joinpath(tablepath,"table_results.csv"),"w")
    io_fitting = open(joinpath(tablepath,"table_fitting.csv"),"w")
    io_gevp    = open(joinpath(tablepath,"table_gevp.csv"),"w")

    # a - eta' notation
    #write(io_resultsMR,L"Label;$\beta$;$N_t$;$N_l$;$am_0^{\rm f}$;$am_0^{\rm as}$;$am_a$;$\Delta am_a$;$am_{\eta'}$;$\Delta am_{\eta'};$am_{\pi}^{(f)}$;$\Delta am_{\pi}^{(f)};$am_{\pi}^{(as)}$;$\Delta am_{\pi}^{(as)}$;$am_{\rho}^{(f)}$;$\Delta am_{\rho}^{(f)};$am_{\rho}^{(as)}$;$\Delta am_{\rho}^{(as)}$","\n")
    #write(io_results,L"Label;$\beta$;$N_t$;$N_l$;$am_0^{\rm f}$;$am_0^{\rm as}$;$am_a$;$am_{\eta'}$;$am_{\rm PS}$;$am_{\rm ps}$;$am_{\rm V}$;$am_{\rm v}$","\n")
    #write(io_fitting,L"	Label;$I_{a}$;$I_{\eta'}$;$I_{\rm{PS}}$;$I_{\rm{ps}}$;$I_{\rm{V}}$;$I_{\rm{v}}$;$N_{\rm exp}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$\\&&&&&&&&$ a$&$\eta'$&${\rm PS}$&${\rm ps}$&$ {\rm V}$&$ {\rm v}$","\n")
    #write(io_gevp,L"Label;$t_0$;$N^{\rm smear}_{a,\eta'}$;$N^{\rm smear}_{\pi(f)}$;$N^{\rm smear}_{\pi(as)}$;$N^{\rm smear}_{\rho(f)}$;$N^{\rm smear}_{\rho(as)}$","\n")

    # eta'_l - eta'_h notation
    write(io_resultsMR,L"Label;$\beta$;$N_t$;$N_l$;$am_0^{\rm f}$;$am_0^{\rm as}$;$am_{\eta^{\prime}_l}$;$\Delta am_{\eta^{\prime}_l}$;$am_{\eta^{\prime}_h}$;$\Delta am_{\eta^{\prime}_h};$am_{\pi}^{(f)}$;$\Delta am_{\pi}^{(f)};$am_{\pi}^{(as)}$;$\Delta am_{\pi}^{(as)}$;$am_{\rho}^{(f)}$;$\Delta am_{\rho}^{(f)};$am_{\rho}^{(as)}$;$\Delta am_{\rho}^{(as)}$","\n")
    write(io_results,L"Label;$\beta$;$N_t$;$N_l$;$am_0^{\rm f}$;$am_0^{\rm as}$;$am_{\eta^{\prime}_l}$;$am_{\eta^{\prime}_h}$;$am_{\rm PS}$;$am_{\rm ps}$;$am_{\rm V}$;$am_{\rm v}$","\n")
    write(io_fitting,L"	Label;$I_{\eta^{\prime}_l}$;$I_{\eta^{\prime}_h}$;$I_{\rm{PS}}$;$I_{\rm{ps}}$;$I_{\rm{V}}$;$I_{\rm{v}}$;$N_{\rm exp}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$;$\chi^2 / N_{\rm d.o.f.}$\\&&&&&&&&$ \eta^{\prime}_l$&$\eta^{\prime}_h$&${\rm PS}$&${\rm ps}$&$ {\rm V}$&$ {\rm v}$","\n")
    write(io_gevp,L"Label;$t_0$;$N^{\rm smear}_{\eta^{\prime}_l,\eta^{\prime}_h}$;$N^{\rm smear}_{\pi(f)}$;$N^{\rm smear}_{\pi(as)}$;$N^{\rm smear}_{\rho(f)}$;$N^{\rm smear}_{\rho(as)}$","\n")

    #check that the number of datasets match
    @assert first(size(parameters)) == first(size(parameters_fitting)) == first(size(corrfitter_results)) 

    parse_smearing_indices(ops) = parse.(Int,split(replace(ops,r"[()]"=>""),','))
    function _get_measurement_id(ensemble,channel,table)
        isensemble = contains.(parameters[:,1],ensemble)
        ischannel  = contains.(parameters[:,2],channel)
        @assert sum(isensemble .* ischannel) == 1
        id = findfirst(isequal(1),isensemble .* ischannel)
        return id
    end

    ensembles = unique(parameters[:,1])
    channels = unique(parameters[:,2])

    for ensemble in ensembles
        
        id = _get_measurement_id(ensemble,"g5_singlet",parameters)

        nops = parameters[id,6]
        β, T, L, mf, mas = corrfitter_results[id,7], corrfitter_results[id,3], corrfitter_results[id,4], corrfitter_results[id,5], corrfitter_results[id,6]

        ma, Δma = corrfitter_results[id,8],  corrfitter_results[id,9] 
        mη, Δmη = corrfitter_results[id,10], corrfitter_results[id,11]
        maHR, mηHR = corrfitter_results_HR[id,8], corrfitter_results_HR[id,9]
        χ2dofa, χ2dofη = corrfitter_results[id,12], corrfitter_results[id,13]

        idπF = _get_measurement_id(ensemble,"g5_nonsinglet_FUN",parameters)
        idρF = _get_measurement_id(ensemble,"g1_nonsinglet_FUN",parameters)
        idπA = _get_measurement_id(ensemble,"g5_nonsinglet_AS",parameters)
        idρA = _get_measurement_id(ensemble,"g1_nonsinglet_AS",parameters)
        mπFHR, mπF, ΔmπF, χ2dofπF = corrfitter_results_HR[idπF,8], corrfitter_results[idπF,8],  corrfitter_results[idπF,9], corrfitter_results[idπF,12] 
        mπAHR, mπA, ΔmπA, χ2dofπA = corrfitter_results_HR[idπA,8], corrfitter_results[idπA,8],  corrfitter_results[idπA,9], corrfitter_results[idπA,12] 
        mρFHR, mρF, ΔmρF, χ2dofρF = corrfitter_results_HR[idρF,8], corrfitter_results[idρF,8],  corrfitter_results[idρF,9], corrfitter_results[idρF,12] 
        mρAHR, mρA, ΔmρA, χ2dofρA = corrfitter_results_HR[idρA,8], corrfitter_results[idρA,8],  corrfitter_results[idρA,9], corrfitter_results[idρA,12] 

        χ2dofπF = round(χ2dofπF,sigdigits=2)
        χ2dofπA = round(χ2dofπA,sigdigits=2)
        χ2dofρF = round(χ2dofρF,sigdigits=2)
        χ2dofρA = round(χ2dofρA,sigdigits=2)
        χ2dofa = round(χ2dofa,sigdigits=2)
        χ2dofη = round(χ2dofη,sigdigits=2)

        t0a, t0η, t1a, t1η, symmetry, Nexp = parameters_fitting[id,3:8]
        t0πF, t0πF, t1πF, t1πF, symmetryπF, NexpπF = parameters_fitting[idπF,3:8]
        t0πA, t0πA, t1πA, t1πA, symmetryπA, NexpπA = parameters_fitting[idπA,3:8]
        t0ρF, t0ρF, t1ρF, t1ρF, symmetryρF, NexpρF = parameters_fitting[idρF,3:8]
        t0ρA, t0ρA, t1ρA, t1ρA, symmetryρA, NexpρA = parameters_fitting[idρA,3:8]
        @assert NexpπF == NexpπA == NexpρF == NexpρA == Nexp
    
        t0_gevp_η,  binsize_η,  deriv_η,  ops_η = parameters[id,3:6]
        t0_gevp_πF, binsize_πF, deriv_πF, ops_πF = parameters[id,3:6]
        t0_gevp_ρF, binsize_ρF, deriv_ρF, ops_ρF = parameters[id,3:6]
        t0_gevp_πA, binsize_πA, deriv_πA, ops_πA = parameters[id,3:6]
        t0_gevp_ρA, binsize_ρA, deriv_ρA, ops_ρA = parameters[id,3:6]

        @assert t0_gevp_η == t0_gevp_πF == t0_gevp_ρF == t0_gevp_πA == t0_gevp_ρA
        @assert binsize_η == binsize_πF == binsize_ρF == binsize_πA == binsize_ρA
        t0_gevp = t0_gevp_ρA
        binsize = binsize_η

        ops_η  = Tuple(getindex(Nsmear,filter(x-> x <= length(Nsmear),parse_smearing_indices(ops_η))))
        ops_πF = Tuple(getindex(Nsmear,filter(x-> x <= length(Nsmear),parse_smearing_indices(ops_πF))))
        ops_ρF = Tuple(getindex(Nsmear,filter(x-> x <= length(Nsmear),parse_smearing_indices(ops_ρF))))
        ops_πA = Tuple(getindex(Nsmear,filter(x-> x <= length(Nsmear),parse_smearing_indices(ops_πA))))
        ops_ρA = Tuple(getindex(Nsmear,filter(x-> x <= length(Nsmear),parse_smearing_indices(ops_ρA))))

        write(io_resultsMR,"$ensemble;$β;$T;$L;$mf;$mas;$ma;$Δma;$mη;$Δmη;$mπF;$ΔmπF;$mπA;$ΔmπA;$mρF;$ΔmρF;$mρA;$ΔmρA\n")
        write(io_results,"$ensemble;$β;$T;$L;$mf;$mas;$maHR;$mηHR;$mπFHR;$mπAHR;$mρFHR;$mρAHR\n")
        write(io_fitting,"$ensemble;($t0a,$t1a);($t0η,$t1η);($t0πF,$t1πF);($t0πA,$t1πA);($t0ρF,$t1ρF);($t0ρA,$t1ρA);$Nexp;$χ2dofπF;$χ2dofπA;$χ2dofρF;$χ2dofρA;$χ2dofa;$χ2dofη\n")
        write(io_gevp,"$ensemble;$t0_gevp;$ops_η;$ops_πF;$ops_ρF;$ops_πA;$ops_ρA\n")
    end
    close(io_resultsMR)
    close(io_results)
    close(io_fitting)
    close(io_gevp)
end