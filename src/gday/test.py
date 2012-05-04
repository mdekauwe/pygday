import numpy as np
from mate import *
import sys



def calc_npp(M, control, params, state, fluxes, met_data):
    
    state.lai = (params.slainit * const.M2_AS_HA /
                            const.KG_AS_TONNES / params.cfracts *
                            state.shoot)
    
    # Specific LAI (m2 onesided/kg DW)
    state.sla = params.slainit
    
    
    year = str(control.startyear)
    month = str(control.startmonth)
    day = str(control.startday)
    datex = datetime.datetime.strptime((year + month + day), "%Y%m%d")
        
    npp = np.zeros(0)
    for project_day in xrange(365):
    
        state.shootnc = state.shootn / state.shoot
        state.ncontent = (state.shootnc * params.cfracts /
                                state.sla * const.KG_AS_G)
        daylen = day_length(datex, params.latitude)
        state.wtfac_root = 1.0
        #state.lai = laidata[project_day]
    
    
        if float_lt(state.lai, params.lai_cover):
            frac_gcover = state.lai / params.lai_cover
        else:
            frac_gcover = 1.0
    
        state.light_interception = ((1.0 - math.exp(-params.kext *
                                            state.lai / frac_gcover)) *
                                            frac_gcover)

        M.calculate_photosynthesis(project_day, daylen)
        
        npp = np.append(npp, fluxes.npp_gCm2)
        datex += datetime.timedelta(days=1)
    return npp.sum()
    
if __name__ == "__main__":


    import matplotlib.pyplot as plt
    from file_parser import initialise_model_data
    from utilities import float_lt, day_length
    import datetime
    
    jm = np.linspace(10, 200.0, 20)
    vc = jm/2.0
    
    npp = np.zeros(0)
    for i in xrange(len(jm)):
    
    
        fname = "/Users/mdekauwe/research/NCEAS_face/GDAY_duke_simulation/params/NCEAS_dk_youngforest.cfg"
        (control, params, state, files,
            fluxes, met_data,
                print_opts) = initialise_model_data(fname, DUMP=False)
        
        params.jmaxn = jm[i]
        params.vcmaxn = vc[i]
        print params.vcmaxn
        
        M = Mate(control, params, state, fluxes, met_data)
        
        sum_npp = calc_npp(M, control, params, state, fluxes, met_data)
        npp = np.append(npp, sum_npp)
    
        
        
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(jm, vc, npp, zdir='z', s=20, c='b')    
    ax.set_xlabel("Jmax")
    ax.set_ylabel("Vcmax")
    ax.set_zlabel("NPP")
    plt.show()
    print npp