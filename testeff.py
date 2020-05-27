# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # DAMIC100 Efficiency notebook
# %% [markdown]
# This notebook sums up the study about the efficiency calculation. It includes the processing steps and the cuts made to compute the efficiency. 
# 
# %% [markdown]
# ## processing

# %%


# %% [markdown]
# ## analysis

# %%
# import part
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 

# %% [markdown]
# ### import of the data:
# %% [markdown]
# The data were processed in several runs, we import them one by one and combine them

# %%
dfrun1rec = pd.read_pickle('/Users/gaior/DAMIC/code/efficiency/20200515/run1/pkl/run1.pkl')
dfrun1sim = pd.read_pickle('/Users/gaior/DAMIC/code/efficiency/20200515/run1/pkl/run1_sim.pkl')

dfrun2rec = pd.read_pickle('/Users/gaior/DAMIC/code/efficiency/20200515/run2/pkl/run2.pkl')
dfrun2sim = pd.read_pickle('/Users/gaior/DAMIC/code/efficiency/20200515/run2/pkl/run2_sim.pkl')

dfrun3rec = pd.read_pickle('/Users/gaior/DAMIC/code/efficiency/20200515/run3/pkl/run3.pkl')
dfrun3sim = pd.read_pickle('/Users/gaior/DAMIC/code/efficiency/20200515/run3/pkl/run3_sim.pkl')


# %%
# merge the files:
dfall = pd.concat([dfrun1rec, dfrun2rec, dfrun3rec])
dfallsim = pd.concat([dfrun1sim, dfrun2sim, dfrun3sim])
#print (dfall.shape[0])


# %%
recpositioncut = 'centery < 42 & centery > 1 & (centerx < 8250 & centerx > 4400) & (simdistx == 0 | (simdistx > 5) & (simdisty > 2) )'
simpositioncut = 'simy < 42 & simy > 1 & (simx < 8250 & simx > 4400) & ((simdistx > 5) & (simdisty > 2) )'
maskcut = 'is_masked == 0 & touchmask == 0 & success ==1'
llcut = 'll_14 < 90'
qmaxcut = 'qmax/(ene1*1000./3.77) > 0.2'
basecuts = maskcut + ' & ' +  llcut+ ' & ' +  qmaxcut
radoncut = " (RUNID<2564 | RUNID>2566) &  (RUNID< 2902 | RUNID> 2903) & (RUNID<3267 | RUNID>3336) & (RUNID<3353 | RUNID>3419) & (RUNID<3654 | RUNID> 3657) & (RUNID<3764 | RUNID>3767) & (RUNID<3826 | RUNID>3853) & (RUNID<3868 | RUNID>3874) & (RUNID<3913 | RUNID>3921) & (RUNID<4003 | RUNID > 4007) & (RUNID<4207 | RUNID > 4212)"
badimage = 'RUNID!=2473 & RUNID!=2479 & RUNID!=2482 & RUNID!=2559 & RUNID!=2577 & RUNID!=2611 & RUNID!=2623  & RUNID!=2829 & RUNID!=2843 & RUNID!=2849 & RUNID!=2853 &  RUNID!=2902 & RUNID!=2927 & RUNID !=3003 & RUNID!=3011 & RUNID!=3018 & RUNID!=3020 & RUNID!=3059 & RUNID!=3112 & RUNID!=3203 & RUNID!=3250 & RUNID!=3332 & RUNID!=3345 & RUNID!=3417 & RUNID!=3453 & RUNID!=3473 & RUNID!=3483 & RUNID!=3536 & RUNID!=3545 & RUNID!=3584 & RUNID!=3634 & RUNID!=3636 & RUNID!=3637 & RUNID!=3638 & RUNID!=3639 & RUNID!=3654 & RUNID!=3655 & RUNID!=3656 & RUNID!=3685  & RUNID!=3698 & RUNID!=3707 & RUNID!=3751 & RUNID!=3807 & RUNID!=3852 & RUNID!=3905 & RUNID!=3931 & RUNID!=3958 & RUNID!=4011 & RUNID!=4044 & RUNID!=4074 & RUNID!=4083 & RUNID!=4127'


negativepixelimage = '(RUNID!=3024 | EXTID!=11) &(RUNID!=3029 | EXTID!=11) & (RUNID!=3029 | EXTID!=12)  & (RUNID!=3125 | EXTID!=11) & (RUNID!=3125 | EXTID!=12) &(RUNID!=3126 | EXTID!=11) & (RUNID!=3126 | EXTID!=12) & (RUNID!=3537 | EXTID!=2) & (RUNID!=3537 | EXTID!=12)  & (RUNID!=3538 | EXTID!=2) & (RUNID!=3538 | EXTID!=12)  & (RUNID!=3539 | EXTID!=2) & (RUNID!=3539 | EXTID!=12)  & (RUNID!=3540 | EXTID!=2) &(RUNID!=3540 | EXTID!=12) &(RUNID!=3541 | EXTID!=2) & (RUNID!=3541 | EXTID!=12) &  (RUNID!=3542 | EXTID!=2) & (RUNID!=3542 | EXTID!=12) & (RUNID!=3543 | EXTID!=2) & (RUNID!=3543 | EXTID!=12) & (RUNID!=3544 | EXTID!=2) & (RUNID!=3546 | EXTID!=12) & (RUNID!=3548 | EXTID!=12) & (RUNID!=3598 | EXTID!=12) &(RUNID!=3657 | EXTID!=4) & (RUNID!=3657| EXTID!=6) & (RUNID!=3657 | EXTID!=11) & (RUNID!=3657 | EXTID!=12)'


# %%
dfall.columns


# %%
dfall = dfall.query(negativepixelimage)
dfallsim = dfallsim.query(negativepixelimage)
dfall = dfall.query(radoncut)
dfallsim = dfallsim.query(radoncut)
dfall = dfall.query(badimage)
dfallsim = dfallsim.query(badimage)

dfall = dfall.query(basecuts)
dfall = dfall.query(recpositioncut + ' & multirows == 0')
dfallsim = dfallsim.query(simpositioncut)


# %%
# function definition of the efficiency calculation:
def geteff(dfrec, dfsim, bins, cutrec, cutsim):    
    dfsel = dfrec.query(cutrec)   
    dfsel_ri = dfsel.reset_index()
    dfsel_ri['ebinned'] = pd.cut(dfsel_ri['ene1'],bins)
    ebinned = dfsel_ri.groupby('ebinned') 
    dfselsim = dfsim.query(cutsim)
    dfselsim_ri = dfselsim.reset_index()

    dfselsim_ri['ebinned'] = pd.cut(dfselsim_ri['sime'],bins)
    ebinnedsim = dfselsim_ri.groupby('ebinned')

    recnr = ebinned.count().ene1.values
    simnr = ebinnedsim.count().sime.values
    eff = recnr/simnr
    err_eff = eff_errors(recnr, simnr, 'poisson')
    cbins = (bins[1:] + bins[:-1]) /2
    return [cbins,eff,err_eff]


def geteffsim(dfrec, dfsim, bins, cutrec, cutsim):    
    dfsel = dfrec.copy()
    dfselsim = dfsim.copy()
    dfsel = dfrec.query(cutrec + "& sime >0")
    dfsel = dfsel.reset_index()
    dfsel['ebinned'] = pd.cut(dfsel['sime'],bins)
    ebinned = dfsel.groupby('ebinned') 
    dfselsim = dfsim.query(cutsim)
    dfselsim = dfselsim.reset_index()

    dfselsim['ebinned'] = pd.cut(dfselsim['sime'],bins)
    ebinnedsim = dfselsim.groupby('ebinned')

    recnr = ebinned.count().ene1.values
    simnr = ebinnedsim.count().sime.values
    
    # eff = recnr/simnr
    # #print (eff)
    # err_eff = eff_errors(recnr, simnr, 'poisson')
    # cbins = (bins[1:] + bins[:-1]) /2
    # return [cbins,eff,err_eff]



def eff_errors(k, N, method):
    k = k.astype(float)
    N = N.astype(float)
    if method == 'binomial':
        return (1/N)*np.sqrt(k*(1-(k/N) ))
    if method == 'poisson':
        return np.sqrt( (k*(N+k))/N**3 )
    if method == 'bayes':
	    return np.sqrt( ((k+1)*(k+2)) / ((N+2)*(N+3)) - (k+1)**2/(N+2)**2 )


# %%
bins = np.linspace(0.03,1,50)
dfrec = dfall.copy()
dfsim = dfallsim.copy()
#[cbins, eff_p, err_eff_p] = geteff(dfrec,dfsim,bins,"RUNID!=0","RUNID!=0")
#[cbins, eff_r, err_eff_r] = geteffsim(dfrec,dfsim,bins,"RUNID!=0","RUNID!=0")
#print (dfrec.shape[0])


# %%
# plt.errorbar(cbins,eff_p,yerr=err_eff_p,label='rec')
# plt.errorbar(cbins,eff_r,yerr=err_eff_r,label='sim')
# plt.xlabel('energy ene1 [keV]')
# plt.ylabel('efficiency')
# plt.legend()
# plt.show()


# %%
def eff_function(x,a,b,c,d,e,f,g):
    eff = (1./(1+np.exp(-(x-a)*b))  - c ) *((d*x +e)+(f*np.exp(-g*x)))
    return eff
oldparam = [7.33799254e-02, 7.44618389e+01, 4.72489658e-02, -1.11708055e-03,9.28825553e-01,  4.59657114e-02,  2.15910338e+00]
print (oldparam)
#plt.plot(cbins,eff_function(cbins,*oldparam) )
#plt.show()


# %%
bins = np.linspace(0.03,4,150)
dfrec = dfall.copy()
dfsim = dfallsim.copy()

geteffsim(dfrec,dfsim,bins,"dll< -22","RUNID!=0")
#[cbins, eff_r, err_eff_r] = geteffsim(dfrec,dfsim,bins,"dll< -22","RUNID!=0")
# [cbins, eff_r2, err_eff_r2] = geteffsim(dfrec,dfsim,bins,"dll< -22  & npix1p6 >1","RUNID!=0")


# plt.errorbar(cbins,eff_r, yerr=err_eff_r,label='22')
# plt.errorbar(cbins,eff_r2, yerr=err_eff_r2,label='22')

# #plt.errorbar(cbins,eff_26,yerr=err_eff_26,label='dll < -22 npix1p6')
# plt.plot(cbins,eff_function(cbins,*oldparam) )
# plt.xlim(0,3)
# plt.ylim(0.8,0.95)


# %%



# %%


