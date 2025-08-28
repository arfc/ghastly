import openmc
import openmc.deplete
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import colormaps
from armi import configure
from armi.nucDirectory import nuclideBases



def plot_isos(title, ylabel, xlabel, fig_fname, dep_fname, mat_id, 
              nuc_list=['U235', 'U236', 
                        'Pu238', 'Pu239', 'Pu240', 'Pu241', 
                        'Cs134', 'Cs137'], 
              nuc_colors = ['turquoise', 'teal',
                            'plum', 'orchid', 'mediumorchid', 'purple',
                            'lightcoral', 'firebrick'], 
              t_units='d', m_units='g/cm3', dpi=600):
    '''
    plot isotope concentrations over time, using depletion results h5 file
    and materials xml file
    '''

    res = openmc.deplete.Results(dep_fname)
    t = res.get_times(time_units=t_units)

    comp = {}

    for i, n in enumerate(nuc_list):
        _, conc = res.get_mass(mat=mat_id, nuc=n, 
                                   mass_units=m_units, time_units=t_units)
        comp[n] = {}
        comp[n]['conc'] = conc
        comp[n]['color'] = nuc_colors[i]

    for k, v in list(comp.items()):
        plt.plot(t, v['conc'] , label = k, 
                 color = v['color'], ls='', marker='.')
    plt.legend()
    plt.title(title)
    plt.ylabel(ylabel + '['+m_units+']')
    plt.xlabel(xlabel + '['+t_units+']')
    plt.xlim(t[0], t[-1])
    plt.savefig(fname=fig_fname, dpi=dpi)
    plt.close()


def compare_isos(title, xlabel, fig_fname, 
                 dep_fname1, mat_id1, name1, dep_fname2, mat_id2, name2, 
                 nuc_list=['U235', 'U236', 
                           'Pu238', 'Pu239', 'Pu240', 'Pu241', 
                           'Cs134', 'Cs137'], 
                 nuc_colors = ['turquoise', 'teal',
                               'plum', 'orchid', 'mediumorchid', 'purple',
                               'lightcoral', 'firebrick'], 
                 t_units='d', m_units='g/cm3', dpi=600):
    '''
    compare two isotopic concentrations over time w/ 
    relative and absolute difference.  relative difference will
    divide by the first isotope in the input arguments.  assumes they
    have the same time steps
    '''

    res1 = openmc.deplete.Results(dep_fname1)
    res2 = openmc.deplete.Results(dep_fname2)
    t = res1.get_times(time_units=t_units)

    diff = {}

    for i, n in enumerate(nuc_list):
        _, conc1 = res1.get_mass(mat=mat_id1, nuc=n, 
                                   mass_units=m_units, time_units=t_units)
        _, conc2 = res2.get_mass(mat=mat_id2, nuc=n, 
                                   mass_units=m_units, time_units=t_units)
        diff[n] = {}
        diff[n]['abs'] = conc1-conc2
        diff[n]['rel'] = (conc1-conc2)/conc1
        diff[n]['color'] = nuc_colors[i]

    for k, v in list(diff.items()):
        plt.plot(t, v['abs'] , label = k, 
                 color = v['color'], ls='', marker='.')
    plt.legend()
    plt.suptitle(title)
    plt.title('Absolute Difference: ' + name1 + '-' + name2)
    plt.ylabel('Absolute Difference ' + '['+m_units+']')
    plt.xlabel(xlabel + '['+t_units+']')
    plt.xlim(t[0], t[-1])
    plt.savefig(fname=fig_fname+'-abs', dpi=dpi)
    plt.close()

    for k, v in list(diff.items()):
        plt.plot(t, 100*v['rel'] , label = k, 
                 color = v['color'], ls='', marker='.')
    plt.legend()
    plt.suptitle(title)
    plt.title('Relative Difference: (' +name1+ '-' +name2+')/('+name1+')')
    plt.ylabel('Relative Difference [%]')
    plt.xlabel(xlabel + '['+t_units+']')
    plt.xlim(t[0], t[-1])
    plt.savefig(fname=fig_fname+'-rel', dpi=dpi)
    plt.close()

    
    

def check_converge(title, ylabel, xlabel, fig_fname, dep_fnames, mat_ids,
                   step_colors,
              nuc_list=['U235', 'U236', 
                        'Pu238', 'Pu239', 'Pu240', 'Pu241', 
                        'Cs134', 'Cs137'],
              t_units='d', m_units='g/cm3', dpi=600):
    '''
    given a list of depletion results in order
    of ascending iteration and the nuclide to compare, plot
    nuc over time for each iteration, then provide information such as Cs ratio,
    MOL, EOL concentrations/relative differences compared to the most recent
    iteration

    assumes first iteration is 0, all t steps the same
    '''
    res = openmc.deplete.Results(dep_fnames[0])
    t = res.get_times(time_units = t_units)
    N = len(dep_fnames)

    isos = {}
    for i, n in enumerate(nuc_list):
        isos[n] = {}
        for j, dep_file in enumerate(dep_fnames):
            res = openmc.deplete.Results(dep_file)
            _, conc = res.get_mass(mat=mat_ids[j], nuc=n, 
                                   mass_units=m_units, time_units=t_units)
            isos[n][j] = {}
            isos[n][j]['conc'] = conc
        for l in range(N-1):
            isos[n][l]['uncert'] = ((isos[n][l]['conc']-isos[n][(N-1)]['conc'])
                                    /isos[n][(N-1)]['conc'])
    mid = int(len(isos['U235'][0]['conc'])/2)

    for m in range(N):
        print(str(m) + ":")
        print("step 0: " + 
              str(isos['Cs134'][m]['conc'][0]/isos['Cs137'][m]['conc'][0]))
        print("mid step: " + 
              str(isos['Cs134'][m]['conc'][mid]/isos['Cs137'][m]['conc'][mid]))
        print("end step: "+
              str(isos['Cs134'][m]['conc'][-1]/isos['Cs137'][m]['conc'][-1]))
        print()

    
    for nuc in nuc_list:
        for k in range(N):
            plt.plot(t, isos[nuc][k]['conc'], label = k, 
                     color = step_colors[k], ls='', marker='.')
        plt.legend()
        plt.suptitle(title)
        plt.title(nuc)
        plt.ylabel(ylabel + ' ['+m_units+']')
        plt.xlabel(xlabel + ' ['+t_units+']')
        plt.xlim(t[0], t[-1])
        plt.savefig(fname=fig_fname + nuc, dpi=dpi)
        plt.close()

    for nuc in nuc_list:
        for k in range(N-1):
            plt.plot(t, 100*isos[nuc][k]['uncert'], 
                     label = str(k)+' vs '+str(N-1), 
                     color = step_colors[k], ls='', marker='.')
        plt.legend()
        plt.suptitle(title)
        plt.title(nuc)
        plt.xlabel(xlabel + ' ['+t_units+']')
        plt.ylabel('Relative Difference [%]')
        plt.xlim(t[0], t[-1])
        plt.savefig(fname=fig_fname + nuc +'_rel', dpi=dpi)
        plt.close()


def flux_fission(sp_fname, shape, out_fname):
    '''
    plot flux, given an OpenMC statepoint file
    '''

    sp = openmc.StatePoint(sp_fname)
    tally = sp.get_tally(scores=['flux'])
    flux = tally.get_slice(scores=['flux'])
    flux.mean.shape = shape
    flux.std_dev.shape = shape
    fission = tally.get_slice(scores=['fission'])
    fission.mean.shape = shape
    fission.std_dev.shape = shape

    flux_zip = np.array(list([[list(k) for k in zip(*j)] 
                              for j in zip(*i)] 
                             for i in zip(flux.mean, flux.std_dev)))
    flux_zip.shape = (8, 2)
    fission_zip = np.array(list([[list(k) for k in zip(*j)] 
                                 for j in zip(*i)] 
                                for i in zip(fission.mean, fission.std_dev)))
    fission_zip.shape = (8, 2) 
    np.savetxt(out_fname+'_flux.csv', flux_zip)
    np.savetxt(out_fname+'_fission.csv', fission_zip)



def plot_on_nucchart(dep_name = 
                     'inf_lat_dep/passwise/iter4/depletion_results.h5',
                     out_name='_chart.png'):
    '''
    '''
    configure(permissive=True)
    res = openmc.deplete.Results(dep_name)
    time = res.get_times(time_units='d')

    comp = {}
    f = open('inf_lat_dep/passwise/nuc_list.txt', 'r')
    nuc_w_data = []
    for n in f.readlines():
        nuc_w_data.append(n.strip())
    nuc_in_res = list(res[0].index_nuc.keys())

    nuc_list = list(set(nuc_w_data) & set(nuc_in_res)) 
    for i, n in enumerate(nuc_list):
        _, conc = res.get_mass(mat='1', nuc=n, 
                                   mass_units='g', time_units='d')
        comp[n] = {}
        comp[n]['conc'] = conc
        zam = openmc.data.zam(n)
        comp[n]['N'] = zam[1] - zam[0]
        comp[n]['Z'] = zam[0]
    
    isotopics_by_day = {}
    for i, t in enumerate(time):
        isotopics_by_day[str(int(t))] = []
        for k, v in list(comp.items()):
            isotopics_by_day[str(int(t))].append((v['N'],v['Z'],v['conc'][i]))
    for day in list(isotopics_by_day.keys()):
        N_raw, Z_raw, C_raw = zip(*isotopics_by_day[day])
        mass = sum(C_raw)
        C_adj = [v if v > 1.00e-25 else 0.0 for v in C_raw/mass]
        N = []
        Z = []
        C = []
        for i, c in enumerate(C_adj):
            if c == 0.0:
                pass
            else:
                N.append(N_raw[i])
                Z.append(Z_raw[i])
                C.append(c*100)
        plt.figure(figsize=(12, 8))
        plt.scatter(N, Z, c=C, norm='log', marker="s", vmin=1e-23, vmax=1.0)
        plt.set_cmap('magma')
        plt.colorbar()
        plt.xlim(0, 170)
        plt.ylim(0, 110)
        plt.title("Passwise BCC, i4: Composition in wt% on day "+str(day))
        plt.xlabel("Number of neutrons (N)")
        plt.ylabel("Number of protons (Z)")
        plt.savefig(((4-len(day))*'0')+day+out_name, dpi=600)
        plt.close()

'''
# iteration 0 (other sims branch from here)
plot_isos("Depleting Pebble Isotopics vs Time: Initial BCC Model (All Fresh), i0", 'Concentration', 
          'Time', 'iteration-0', 
          'inf_lat_dep/iter0/depletion_results.h5', '1')

# distinct corners
#i1 through i4:
plot_isos("Depleting Pebble Isotopics: Passwise-BCC Model, i1", 
          'Concentration', 
          'Time', 'passwise-i1', 
          'inf_lat_dep/dep-center/passwise/iter1/depletion_results.h5',
          '13')

plot_isos("Depleting Pebble Isotopics: Passwise-BCC Model, i2", 
          'Concentration', 
          'Time', 'passwise-i2', 
          'inf_lat_dep/dep-center/passwise/iter2/depletion_results.h5',
          '13')

plot_isos("Depleting Pebble Isotopics: Passwise-BCC Model, i3", 
          'Concentration', 
          'Time', 'passwise-i3', 
          'inf_lat_dep/dep-center/passwise/iter3/depletion_results.h5',
          '13')

plot_isos("Depleting Pebble Isotopics: Passwise-BCC Model, i4", 
          'Concentration', 
          'Time', 'passwise-i4', 
          'inf_lat_dep/dep-center/passwise/iter4/depletion_results.h5',
          '13')

# avg corners
#i1 through i4:
plot_isos("Depleting Pebble Isotopics: Corewise-BCC Model, i1", 
          'Concentration', 
          'Time', 'corewise-i1', 
          'inf_lat_dep/dep-center/corewise/iter1/depletion_results.h5',
          '14')

plot_isos("Depleting Pebble Isotopics: Corewise-BCC Model, i2", 
          'Concentration', 
          'Time', 'corewise-i2', 
          'inf_lat_dep/dep-center/corewise/iter2/depletion_results.h5',
          '14')

plot_isos("Depleting Pebble Isotopics: Corewise-BCC Model, i3", 
          'Concentration', 
          'Time', 'corewise-i3', 
          'inf_lat_dep/dep-center/corewise/iter3/depletion_results.h5',
          '14')

plot_isos("Depleting Pebble Isotopics: Corewise-BCC Model, i4", 
          'Concentration', 
          'Time', 'corewise-i4', 
          'inf_lat_dep/dep-center/corewise/iter4/depletion_results.h5',
          '14')

#comparisons:

# passwise vs avg: i1
compare_isos("Passwise vs Corewise BCC Depleting Pebble Isotopics: i1",'Time', 
             'i1-passwise-vs-corewise', 
          'inf_lat_dep/dep-center/passwise/iter1/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/corewise/iter1/depletion_results.h5',
             '14', 'Corewise')

#passwise vs avg: i2
compare_isos("Passwise vs Corewise BCC Depleting Pebble Isotopics: i2",'Time', 
             'i2-passwise-vs-corewise', 
          'inf_lat_dep/dep-center/passwise/iter2/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/corewise/iter2/depletion_results.h5',
             '14', 'Corewise')
#i3
compare_isos("Passwise vs Corewise BCC Depleting Pebble Isotopics: i3",'Time', 
             'i3-passwise-vs-corewise', 
          'inf_lat_dep/dep-center/passwise/iter3/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/corewise/iter3/depletion_results.h5',
             '14', 'Corewise')
#i4

compare_isos("Passwise vs Corewise BCC Depleting Pebble Isotopics: i4",'Time', 
             'i4-passwise-vs-corewise', 
          'inf_lat_dep/dep-center/passwise/iter4/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/corewise/iter4/depletion_results.h5',
             '14', 'Corewise')

#check convergence

#passwise:

fnames1 = ['inf_lat_dep/iter0/depletion_results.h5', 
          'inf_lat_dep/dep-center/passwise/iter1/depletion_results.h5',
          'inf_lat_dep/dep-center/passwise/iter2/depletion_results.h5',
          'inf_lat_dep/dep-center/passwise/iter3/depletion_results.h5',
          'inf_lat_dep/dep-center/passwise/iter4/depletion_results.h5']
mat_ids1 = ['1', '13', '13', '13', '13']
stepcolors= ['turquoise', 'teal', 'orchid', 'purple', 'firebrick']

check_converge("Convergence Check using Depleting Pebble: Passwise BCC Model",
               'Concentration', 'Time', 'passwise-converge', fnames1, mat_ids1,
                   stepcolors)

fnames2 = ['inf_lat_dep/iter0/depletion_results.h5', 
          'inf_lat_dep/dep-center/corewise/iter1/depletion_results.h5',
          'inf_lat_dep/dep-center/corewise/iter2/depletion_results.h5',
          'inf_lat_dep/dep-center/corewise/iter3/depletion_results.h5',
          'inf_lat_dep/dep-center/corewise/iter4/depletion_results.h5']
mat_ids2 = ['1', '14', '14', '14', '14']

check_converge("Convergence Check Using Depleting Pebble: Corewise BCC Corners",
               'Concentration', 'Time', 'corewise-converge', fnames2, mat_ids2,
                   stepcolors)



corewise_sp_fnames = ['bcc_phys/corewise/0days/statepoint.h5', 
                      'bcc_phys/corewise/249days/statepoint.h5',
                      'bcc_phys/corewise/499days/statepoint.h5', 
                      'bcc_phys/corewise/749days/statepoint.h5',
                      'bcc_phys/corewise/999days/statepoint.h5', 
                      'bcc_phys/corewise/1249days/statepoint.h5',
                      'bcc_phys/corewise/1549days/statepoint.h5']
passwise_sp_fnames = ['bcc_phys/passwise/0days/statepoint.h5', 
                      'bcc_phys/passwise/249days/statepoint.h5',
                      'bcc_phys/passwise/499days/statepoint.h5', 
                      'bcc_phys/passwise/749days/statepoint.h5',
                      'bcc_phys/passwise/999days/statepoint.h5', 
                      'bcc_phys/passwise/1249days/statepoint.h5',
                      'bcc_phys/passwise/1549days/statepoint.h5']

out_fnames = ['0days', '249days', '499days', 
              '749days', '999days', '1249days', '1549days']

for i, fname in enumerate(corewise_sp_fnames):
    flux_fission(fname, (2, 2, 2), 'corewise-'+out_fnames[i])

for i, fname in enumerate(passwise_sp_fnames):
    flux_fission(fname, (2, 2, 2), 'passwise-'+out_fnames[i])
'''

plot_on_nucchart()
