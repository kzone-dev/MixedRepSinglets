function _fold_corr!(corr;sign=+1)
    N, T = size(corr)
    # Skip first entry since it does not to be averaged
    # Add 1 to all indices: julia is one-indexed
    for t in 1:div(T,2)
        t1 = t+1
        t2 = T-t+1
        tmp1 = (corr[:,t1] + sign*corr[:,t2])/2
        tmp2 = (sign*corr[:,t1] + corr[:,t2])/2
        corr[:,t1] = tmp1
        corr[:,t2] = tmp2
    end
end
function awi_corr(file)
    AP = h5read(file,"g0g5_g5")
    g5 = h5read(file,"g5")

    # symmetrise the correlator
    _fold_corr!(AP;sign=-1)
    _fold_corr!(g5)

    dAP = correlator_derivative(AP;t_dim=2)
    awi_corr = dAP ./ g5 ./ 2

    N, T = size(AP)

    AWI  = dropdims(mean(awi_corr,dims=1),dims=1)
    ΔAWI = dropdims(std(awi_corr,dims=1),dims=1)/sqrt(N)
    return AWI, ΔAWI
end
@. constfit(x,p) = p[1] + 0*x
function awi_fit(AWI,ΔAWI;tmin,tmax)
    T  = length(AWI)
    t0 = tmin:tmax
    # initial guess at T/2
    p0 = AWI[T÷2]
    fit  = curve_fit(constfit,t0,AWI[t0],ΔAWI[t0].^(-2),[p0])
    # extract mass and estimate error
    mAWI = fit.param[1]
    ΔmAWI = stderror(fit)[1]
    return mAWI, ΔmAWI
end
