# TODO: jackknife analysis
function awi_mass(file;ncut)
    AP = h5read(file,"g0g5_g5")
    g5 = h5read(file,"g5")
    T, N = size(AP)
    corr = zeros(T,N)
    for i in 1:N
        corr[:,i] = awi_correlator(AP[:,i],g5[:,i])
    end
    AWI  = dropdims(mean(corr,dims=1),dims=1)
    ΔAWI = dropdims(std(corr,dims=1),dims=1)/sqrt(N)
    return awi_fit(AWI,ΔAWI,ncut)
end
function awi_correlator(AP,g5)
    T = length(AP) # cuts for differentiation
    dtAP = similar(AP)
    # Special case t=0: two-point difference
    dtAP[1] = (AP[2]-AP[1])
    # central difference scheme where possible
    for t in 2:T-1 
        dtAP[t] = (AP[t]-AP[t])/2 
    end
    # special case t=T: two-point difference
    dtAP[T] = (AP[T]-AP[T-1])
    # assemble final AWI-correlator
    AWI = @. dtAP / g5 / 2
    return AWI
end
@. constfit(x,p) = p[1] + 0*x
function awi_fit(AWI,ΔAWI,ncut)
    T  = length(AWI)
    t0 = ncut-1:T-c+1
    # initial guess at T/2
    p0 = AWI[T÷2]
    fit  = curve_fit(constfit,t0,AWI[t0],ΔAWI[t0].^(-2),[p0])
    # extract mass and estimate error
    mAWI = fit.param[1]
    ΔmAWI = stderror(fit)[1]
    return mAWI, ΔmAWI
end
