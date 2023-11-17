model_1st(t,T,m,f;sign=+1,)  = (exp(-m*t)+sign*exp(-m*(T-t)))*f^2*m/2
model_2nd(t,T,m,f;sign=+1)  = (exp(-m*t)+sign*exp(-m*(T-t)))*f/2
fitmass(fit)=first(coef(fit))
function _aply_cut(ncut,T)
    ncut = isa(ncut,Int) ? (ncut,div(T,2)) : ncut
    @assert isa(ncut,Tuple)
    @assert length(ncut)==2
    return ncut[1]:ncut[2]
end
function fit_corr_1exp(c,w,ncut;sign=+1)
    T = length(c)
    t = collect(1:T)
    cut = _aply_cut(ncut,T)
    if ndims(w)==2
        weight = w[cut,cut]
        weight = (weight + weight')/2
    else
        weight = w[cut]
    end
    @. model(t,p) = model_1st(t-1,T,p[1],p[2];sign)
    fit = curve_fit(model,t[cut],c[cut],weight,ones(2))
    fit.param[2] = abs(fit.param[2])
    return fit, model
end
function fit_corr(c,cov,ncut;sign=+1)
    if ndims(cov)==2
        w  = inv(cov)
        w = (w + w')/2
    else
        w = inv.(cov)
    end
    return fit_corr_1exp(c,w,ncut;sign)
end
function fit_corr_bars(c,cov,ncut;sign=+1)
    if ndims(cov)==2
        Δc = sqrt.(diag(cov))
    else
        Δc = sqrt.(cov)
    end
    fit1, model = fit_corr(c,cov,ncut;sign)
    fit2, model = fit_corr(c+Δc,cov,ncut;sign)
    fit3, model = fit_corr(c-Δc,cov,ncut;sign)
    return fit1, fit2, fit3, model
end