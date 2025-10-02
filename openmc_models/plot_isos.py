import openmc
import openmc.deplete
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import colormaps



def plot_isos(title, ylabel, xlabel, fig_fname, dep_fname, mat_id, 
              nuc_list=['U235', 'U236', 
                        'Pu238', 'Pu239', 'Pu240', 'Pu241', 
                        'Cs134', 'Cs137'], 
              nuc_colors = ['turquoise', 'teal',
                            'plum', 'orchid', 'mediumorchid', 'purple',
                            'lightcoral', 'firebrick'], 
              t_units='d', m_units='g/cm3', dpi=600):
    '''
    Plot isotope concentrations over time, using depletion_results.h5 file
    and materials xml file

    Parameters:
    ----------
    title : str
        Title of figure.
    ylabel : str
        Matplotlib ylabel of figure.
    xlabel : str
        Matplotlib xlabel of figure.
    fig_fname : str
        Name generated figure is saved as.
    dep_fname : str
        Relative filepath to depletion_results.h5 file concentrations
        are pulled from.
    mat_id : str
        Openmc material_id of material to be plotted.
    nuc_list : list of str
        List of Openmc recognizable isotope codes to be plotted.
    nuc_colors : list of str
        List of built-in matplotlib colors, used for plotting.
    t_units : str
        Units code for time units used when reading results.h5.
    m_units : str
        Units code for mass units used when reading results.h5.
    dpi : int
        Desired dpi when saving figure.
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
    Compare two isotopic concentrations over time with
    relative and absolute difference.  Relative difference will
    divide by the average of the two concentrations.

    Parameters:
    ----------
    title : str
        Title of figure.
    ylabel : str
        Matplotlib ylabel of figure.
    xlabel : str
        Matplotlib xlabel of figure.
    fig_fname : str
        Name generated figure is saved as.
    dep_fname(i) : str
        Relative filepath to depletion_results.h5 file concentrations
        are pulled from.
    mat_id(i) : str
        Openmc material_id of material to be plotted.
    name(i) : str
        Label of material, for figure legend.
    nuc_list : list of str
        List of Openmc recognizable isotope codes to be plotted.
    nuc_colors : list of str
        List of built-in matplotlib colors, used for plotting.
    t_units : str
        Units code for time units used when reading results.h5.
    m_units : str
        Units code for mass units used when reading results.h5.
    dpi : int
        Desired dpi when saving figure.

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
        avg = (conc1+conc2)/2
        diff[n]['rel'] = (conc1-conc2)/avg
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
    plt.title('Relative Difference: (' +name1+ '-' +name2+')/(Avg)')
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

    Parameters:
    ----------
    title : str
        Title of figure.
    ylabel : str
        Matplotlib ylabel of figure.
    xlabel : str
        Matplotlib xlabel of figure.
    fig_fname : str
        Name generated figure is saved as.
    dep_fnames : list of str
        List of relative filepath to depletion_results.h5 files concentrations
        are pulled from.
    mat_ids : list of str
        List of Openmc material_ids of each material to be plotted.
    step_colors : list of str
        List of built-in color names to plot each iteration's concentrations
        with.
    nuc_list : list of str
        List of Openmc recognizable isotope codes to be plotted.
    t_units : str
        Units code for time units used when reading results.h5.
    m_units : str
        Units code for mass units used when reading results.h5.
    dpi : int
        Desired dpi when saving figure.

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

    Parameters:
    ----------
    sp_fname : str
        Relative filepath to OpenMC statepoint file.
    shape : tuple of int
        Dimensions of detector mesh.
    out_fname : str
        Name to save output files to.
    '''

    sp = openmc.StatePoint(sp_fname)
    spatial = sp.get_tally(name='spatial')
    flux = spatial.get_slice(scores=['flux'])
    flux.mean.shape = shape
    flux.std_dev.shape = shape
    fission = spatial.get_slice(scores=['fission'])
    fission.mean.shape = shape
    fission.std_dev.shape = shape

    heating = sp.get_tally(name='heating')
    H = heating.get_slice(scores=['heating'])
    H.mean.shape = (1)
    J = 1.602e-19
    P = 165e6 #watts
    bcc_l = 100/(2697**(1/3))
    V = (bcc_l**3)/9
    f = (P/(H.mean[0]*J))/V

    spectrum = sp.get_tally(name='spectrum')
    spectra = spectrum.get_slice(scores=['flux'])
    spectra.mean.shape = (500)
    spectra.std_dev.shape = (500)
    ene = np.logspace(np.log10(1e-5), np.log10(20.0e6), 501)

    flux_zip = np.array(list([[list(k) for k in zip(*j)] 
                              for j in zip(*i)] 
                             for i in zip(f*flux.mean, f*flux.std_dev)))
    flux_zip.shape = (27, 2)
    fission_zip = np.array(list([[list(k) for k in zip(*j)] 
                                 for j in zip(*i)] 
                                for i in zip(f*fission.mean, f*fission.std_dev)))
    fission_zip.shape = (27, 2)
    
    spectra_zip = np.array(list(zip(ene, f*spectra.mean, f*spectra.std_dev)))

    np.savetxt(out_fname+'_flux.csv', flux_zip)
    np.savetxt(out_fname+'_fission.csv', fission_zip)
    np.savetxt(out_fname+'_fluxspectrum.csv', spectra_zip)



def plot_on_nucchart(dep_name = 
                     'inf_lat_dep/corewise/iter4/depletion_results.h5',
                     out_name='_corewise_chart.png'):
    '''

    Parameters:
    ----------
    dep_name : str
        Relative filepath to depletion_results.h5 file.
    out_name : str
        Suffix used when naming auto-generated figures.
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
        plt.title("Corewise BCC, i4: Composition in wt% on day "+str(day))
        plt.xlabel("Number of neutrons (N)")
        plt.ylabel("Number of protons (Z)")
        plt.savefig(((4-len(day))*'0')+day+out_name, dpi=600)
        plt.close()

def compare_on_nucchart(core_dep = 
                        'inf_lat_dep/corewise/iter4/depletion_results.h5', 
                        pass_dep = 
                        'inf_lat_dep/passwise/iter4/depletion_results.h5',
                        out_name='_compare_chart.png'):
    '''
    Parameters:
    ----------
    core_dep : str
        Relative filepath to corewise depletion_results.h5 file.
    pass_dep : str
        Relative filepath to passwise depletion_results.h5 file.
    out_name : str
        Suffix used when naming auto-generated figures.

    '''
    configure(permissive=True)
    pass_res = openmc.deplete.Results(pass_dep)
    core_res = openmc.deplete.Results(core_dep)
    time = pass_res.get_times(time_units='d')

    pass_comp = {}
    core_comp = {}
    f = open('inf_lat_dep/passwise/nuc_list.txt', 'r')
    nuc_w_data = []
    for n in f.readlines():
        nuc_w_data.append(n.strip())
    nuc_in_pass = list(pass_res[0].index_nuc.keys())

    nuc_list = list(set(nuc_w_data) & set(nuc_in_pass))
    for i, n in enumerate(nuc_list):
        _, pass_conc = pass_res.get_mass(mat='1', nuc=n, 
                                         mass_units='g', time_units='d')
        _, core_conc = core_res.get_mass(mat='1', nuc=n,
                                         mass_units='g', time_units='d')
        zam = openmc.data.zam(n)
        pass_comp[n] = {}
        core_comp[n] = {}
        pass_comp[n]['conc'] = pass_conc
        core_comp[n]['conc'] = core_conc
        pass_comp[n]['N'] = zam[1] - zam[0]
        core_comp[n]['N'] = zam[1] - zam[0]
        pass_comp[n]['Z'] = zam[0]
        core_comp[n]['Z'] = zam[0]
    
    diff_by_day = {}
    for i, t in enumerate(time):
        diff_by_day[str(int(t))] = []
        for k, v in list(pass_comp.items()):
            avg_comp = (v['conc'][i]+core_comp[k]['conc'][i])/2
            if avg_comp == 0.0:
                rel_diff = 0.0
            else:
                rel_diff = 100*(v['conc'][i]-core_comp[k]['conc'][i])/avg_comp
            diff_by_day[str(int(t))].append((v['N'], v['Z'], rel_diff))
    
    for day in list(diff_by_day.keys()):
        N_raw, Z_raw, D_raw = zip(*diff_by_day[day])
        N = []
        Z = []
        D = []
        for i, d in enumerate(D_raw):
            if d==0.0:
                pass
            else:
                N.append(N_raw[i])
                Z.append(Z_raw[i])
                D.append(d)
        plt.figure(figsize=(12, 8))
        plt.scatter(N_raw, Z_raw, c=D_raw, norm='symlog', marker="s", vmin=-100, vmax=100)
        plt.set_cmap('magma')
        plt.colorbar()
        plt.xlim(0, 170)
        plt.ylim(0, 110)
        plt.title("i4 BCC Isotopics Relative Difference (Pass-Core)/Avg: day "+str(day))
        plt.xlabel("Number of neutrons (N)")
        plt.ylabel("Number of protons (Z)")
        plt.savefig(((4-len(day))*'0')+day+'_all_log'+out_name, dpi=600)
        plt.close()

        plt.figure(figsize=(12, 8))
        plt.scatter(N, Z, c=D, norm='linear', marker="s", vmin=-100, vmax=100)
        plt.set_cmap('magma')
        plt.colorbar()
        plt.xlim(0, 170)
        plt.ylim(0, 110)
        plt.title("i4 BCC Isotopics Relative Difference (Pass-Core)/Avg: day "+str(day))
        plt.xlabel("Number of neutrons (N)")
        plt.ylabel("Number of protons (Z)")
        plt.savefig(((4-len(day))*'0')+day+'_linear'+out_name, dpi=600)
        plt.close()
        
        half = int(len(N_raw)/2)
        plt.figure(figsize=(12, 8))
        plt.scatter(N_raw[:half], Z_raw[:half], c=D_raw[:half], 
                    norm='symlog', marker="s", vmin=-100, vmax=100)
        plt.set_cmap('magma')
        plt.colorbar()
        plt.xlim(0, 90)
        plt.ylim(0, 60)
        plt.title("i4 BCC Isotopics Relative Difference (Pass-Core)/Avg: day "+str(day))
        plt.xlabel("Number of neutrons (N)")
        plt.ylabel("Number of protons (Z)")
        plt.savefig(((4-len(day))*'0')+day+'_firsthalf_log'+out_name, dpi=600)
        plt.close()

        plt.figure(figsize=(12, 8))
        plt.scatter(N_raw[half:], Z_raw[half:], c=D_raw[half:], 
                    norm='symlog', marker="s", vmin=-100, vmax=100)
        plt.set_cmap('magma')
        plt.colorbar()
        plt.xlim(80, 170)
        plt.ylim(50, 110)
        plt.title("i4 BCC Isotopics Relative Difference (Pass-Core)/Avg: day "+str(day))
        plt.xlabel("Number of neutrons (N)")
        plt.ylabel("Number of protons (Z)")
        plt.savefig(((4-len(day))*'0')+day+'_secondhalf_log'+out_name, dpi=600)
        plt.close()


'''

#iteration 0 (other sims branch from here)
plot_isos("Depleting Pebble Isotopics vs Time: Initial BCC Model (All Fresh), i0", 
          'Concentration', 'Time', 'iteration-0', 
          'inf_lat_dep/iter0/depletion_results.h5', '1')

#  passwise
#i1 through i4:
passwise_files = ['inf_lat_dep/passwise/iter1/depletion_results.h5',
                  'inf_lat_dep/passwise/iter2/depletion_results.h5',
                  'inf_lat_dep/passwise/iter3/depletion_results.h5',
                  'inf_lat_dep/passwise/iter4/depletion_results.h5']
for i, file in enumerate(passwise_files):
    img_name = 'passwise-i'+str(i+1)
    title = 'Depleting Pebble Isotopics: Passwise-BCC Model, i'+str(i+1)
    plot_isos(title, 'Concentration', 'Time', img_name, file, '13')

# corewise
#i1 through i4:
corewise_files = ['inf_lat_dep/corewise/iter1/depletion_results.h5',
                  'inf_lat_dep/corewise/iter2/depletion_results.h5',
                  'inf_lat_dep/corewise/iter3/depletion_results.h5',
                  'inf_lat_dep/corewise/iter4/depletion_results.h5']
for i, file in enumerate(corewise_files):
    img_name = 'corewise-i'+str(i+1)
    title = 'Depleting Pebble Isotopics: Corewise-BCC Model, i'+str(i+1)
    plot_isos(title, 'Concentration', 'Time', img_name, file, '14')

#comparisons:

for i in range(4):
    title = 'Passwise vs Corewise BCC Depleting Pebble Isotopics: i'*str(i+1)
    img_name = 'i'+str(i+1)+'-passwise-vs-corewise'
    passfile = 'inf_lat_dep/passwise/iter'+str(i+1)+'/depletion_results.h5'
    corefile = 'inf_lat_dep/corewise/iter'+str(i+1)+'/depletion_results.h5'
    compare_isos(title, 'Time', img_name, 
                 passfile, '1', 'Passwise',
                 corefile, '1', 'Corewise')

#check convergence

#passwise:

fnames1 = ['inf_lat_dep/iter0/depletion_results.h5', 
          'inf_lat_dep/passwise/iter1/depletion_results.h5',
          'inf_lat_dep/passwise/iter2/depletion_results.h5',
          'inf_lat_dep/passwise/iter3/depletion_results.h5',
          'inf_lat_dep/passwise/iter4/depletion_results.h5']
mat_ids1 = ['1', '1', '1', '1', '1']
stepcolors= ['turquoise', 'teal', 'orchid', 'purple', 'firebrick']

check_converge("Convergence Check using Depleting Pebble: Passwise BCC Model",
               'Concentration', 'Time', 'passwise-converge', fnames1, mat_ids1,
                   stepcolors)

fnames2 = ['inf_lat_dep/iter0/depletion_results.h5', 
          'inf_lat_dep/corewise/iter1/depletion_results.h5',
          'inf_lat_dep/corewise/iter2/depletion_results.h5',
          'inf_lat_dep/corewise/iter3/depletion_results.h5',
          'inf_lat_dep/corewise/iter4/depletion_results.h5']
mat_ids2 = ['1', '1', '1', '1', '1']

check_converge("Convergence Check Using Depleting Pebble: Corewise BCC Corners",
               'Concentration', 'Time', 'corewise-converge', fnames2, mat_ids2,
                   stepcolors)


'''

corewise_sp_fnames = ['bcc_phys/corewise/0days/statepoint.80.h5', 
                      'bcc_phys/corewise/249days/statepoint.80.h5',
                      'bcc_phys/corewise/499days/statepoint.80.h5', 
                      'bcc_phys/corewise/749days/statepoint.80.h5',
                      'bcc_phys/corewise/999days/statepoint.80.h5', 
                      'bcc_phys/corewise/1249days/statepoint.80.h5',
                      'bcc_phys/corewise/1549days/statepoint.80.h5']
passwise_sp_fnames = ['bcc_phys/passwise/0days/statepoint.80.h5', 
                      'bcc_phys/passwise/249days/statepoint.80.h5',
                      'bcc_phys/passwise/499days/statepoint.80.h5', 
                      'bcc_phys/passwise/749days/statepoint.80.h5',
                      'bcc_phys/passwise/999days/statepoint.80.h5', 
                      'bcc_phys/passwise/1249days/statepoint.80.h5',
                      'bcc_phys/passwise/1549days/statepoint.80.h5']

out_fnames = ['0days', '249days', '499days', 
              '749days', '999days', '1249days', '1549days']

for i, fname in enumerate(corewise_sp_fnames):
    flux_fission(fname, (3, 3, 3), 'corewise-'+out_fnames[i])

for i, fname in enumerate(passwise_sp_fnames):
    flux_fission(fname, (3, 3, 3), 'passwise-'+out_fnames[i])


#plot_on_nucchart()

#compare_on_nucchart()
