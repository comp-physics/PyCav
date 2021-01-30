import numpy as np

def get_moments(sols):
    moms = []
    for i, sol in enumerate(sols):
        Nt = len(sol.times)
        mom = np.zeros([Nt,sol.state.Nmom])
        if sol.filter: 
            for j in range(Nt):
                jMin = max(j-sol.Nfilt,0)
                mom[j,:] = sol.state.get_quad(
                        vals=sol.save[jMin:j+1],
                        Nfilt=j+1-jMin)
        else:
            for j, vals in enumerate(sol.save):
                mom[j,:] = sol.state.get_quad(vals=vals)
        moms.append(mom)
        sol.moms = mom
        sols[i] = sol

    return sols
