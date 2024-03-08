using Plots
using DelimitedFiles
using LaTeXStrings
gr(legend=:topleft,frame=:box)

dat2, header2 = readdlm("./output/decayconstants.csv",',',header=true)
dat3, header3 = readdlm("./output/PCACmasses.csv",',',header=true)

mq  = dat3[:,6]
Δmq = dat3[:,7]

names = dat3[:,1]

GMOR  = dat2[:,11]
ΔGMOR = dat2[:,12]

as = collect(1:2:length(names))
f = collect(2:2:length(names))
plt_as  = plot(xlabel=L"m_q^{\rm PCAC}",ylabel=L"f_\pi^2 m_\pi^2") 
plt_fun = plot(xlabel=L"m_q^{\rm PCAC}",ylabel=L"f_\pi^2 m_\pi^2")
plt = plot(xlabel=L"m_q^{\rm PCAC}",ylabel=L"f_\pi^2 m_\pi^2")
scatter!(plt_as ,mq[as],GMOR[as],xerr=Δmq[as],yerr=ΔGMOR[as],label="anti-symmetric")
scatter!(plt_fun,mq[f], GMOR[f] ,xerr=Δmq[f] ,yerr=ΔGMOR[f],label="fundamental")
scatter!(plt,mq[as],GMOR[as],xerr=Δmq[as],yerr=ΔGMOR[as],label="anti-symmetric")
scatter!(plt,mq[f], GMOR[f] ,xerr=Δmq[f] ,yerr=ΔGMOR[f],label="fundamental")

plt_as
plt_fun