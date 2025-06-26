import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt



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
        plt.plot(t, v['conc'] , label = k, color = v['color'])
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
        plt.plot(t, v['abs'] , label = k, color = v['color'])
    plt.legend()
    plt.suptitle(title)
    plt.title('Absolute Difference: ' + name1 + '-' + name2)
    plt.ylabel('Absolute Difference ' + '['+m_units+']')
    plt.xlabel(xlabel + '['+t_units+']')
    plt.xlim(t[0], t[-1])
    plt.savefig(fname=fig_fname+'-abs', dpi=dpi)
    plt.close()

    for k, v in list(diff.items()):
        plt.plot(t, 100*v['rel'] , label = k, color = v['color'])
    plt.legend()
    plt.suptitle(title)
    plt.title('Relative Difference: (' +name1+ '-' +name2+')/('+name1+')')
    plt.ylabel('Relative Difference ' + '[%]')
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

    isos = {}
    for i, n in enumerate(nuc_list):
        isos[n] = {}
        for j, dep_file in enumerate(dep_fnames):
            res = openmc.deplete.Results(dep_file)
            _, conc = res.get_mass(mat=mat_ids[j], nuc=n, 
                                   mass_units=m_units, time_units=t_units)
            isos[n][j] = conc


    for nuc in nuc_list:
        for k in range(len(dep_fnames)):
            plt.plot(t, isos[nuc][k], label = k, color = step_colors[k])
        plt.legend()
        plt.suptitle(title)
        plt.title(nuc)
        plt.ylabel(ylabel + ' ['+m_units+']')
        plt.xlabel(xlabel + ' ['+t_units+']')
        plt.xlim(t[0], t[-1])
        plt.savefig(fname=fig_fname + nuc, dpi=dpi)
        plt.close()





# iteration 0 (other sims branch from here)
plot_isos("Concentration vs Time: Iteration 0", 'Concentration [g/cc]', 
          'Time [days]', 'iter0-isos', 
          'inf_lat_dep/iter0/depletion_results.h5', '1')

# distinct corners
#i1 through i4:
plot_isos("Concentration vs Time: Passwise Corners, i1", 
          'Concentration', 
          'Time', 'iter1psw-isos', 
          'inf_lat_dep/dep-center/distinct_corners/iter1/depletion_results.h5',
          '13')

plot_isos("Concentration vs Time: Passwise Corners, i2", 
          'Concentration', 
          'Time', 'iter2psw-isos', 
          'inf_lat_dep/dep-center/distinct_corners/iter2/depletion_results.h5',
          '13')

plot_isos("Concentration vs Time: Passwise Corners, i3", 
          'Concentration', 
          'Time', 'iter3psw-isos', 
          'inf_lat_dep/dep-center/distinct_corners/iter3/depletion_results.h5',
          '13')

plot_isos("Concentration vs Time: Passwise Corners, i4", 
          'Concentration', 
          'Time', 'iter4psw-isos', 
          'inf_lat_dep/dep-center/distinct_corners/iter4/depletion_results.h5',
          '13')

# avg corners
#i1 through i4:
plot_isos("Concentration vs Time: Core-Average Corners, i1", 
          'Concentration', 
          'Time', 'iter1avg-isos', 
          'inf_lat_dep/dep-center/core_avg/iter1/depletion_results.h5',
          '14')

plot_isos("Concentration vs Time: Core-Average Corners, i2", 
          'Concentration', 
          'Time', 'iter2avg-isos', 
          'inf_lat_dep/dep-center/core_avg/iter2/depletion_results.h5',
          '14')

plot_isos("Concentration vs Time: Core-Average Corners, i3", 
          'Concentration', 
          'Time', 'iter3avg-isos', 
          'inf_lat_dep/dep-center/core_avg/iter3/depletion_results.h5',
          '14')

plot_isos("Concentration vs Time: Core-Average Corners, i4", 
          'Concentration', 
          'Time', 'iter4avg-isos', 
          'inf_lat_dep/dep-center/core_avg/iter4/depletion_results.h5',
          '14')

#comparisons:

# passwise vs avg: i1
compare_isos("Passwise vs Core-Averaged Corners: i1",'Time', 'i1-compare', 
          'inf_lat_dep/dep-center/distinct_corners/iter1/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/core_avg/iter1/depletion_results.h5',
             '14', 'Core_Averaged')

#passwise vs avg: i2
compare_isos("Passwise vs Core-Averaged Corners: i2",'Time', 'i2-compare', 
          'inf_lat_dep/dep-center/distinct_corners/iter2/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/core_avg/iter2/depletion_results.h5',
             '14', 'Core_Averaged')
#i3
compare_isos("Passwise vs Core-Averaged Corners: i3",'Time', 'i3-compare', 
          'inf_lat_dep/dep-center/distinct_corners/iter3/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/core_avg/iter3/depletion_results.h5',
             '14', 'Core_Averaged')
#i4

compare_isos("Passwise vs Core-Averaged Corners: i4",'Time', 'i4-compare', 
          'inf_lat_dep/dep-center/distinct_corners/iter4/depletion_results.h5',
          '13', 'Passwise',
             'inf_lat_dep/dep-center/core_avg/iter4/depletion_results.h5',
             '14', 'Core-_Averaged')

#check convergence

#passwise:

fnames1 = ['inf_lat_dep/iter0/depletion_results.h5', 
          'inf_lat_dep/dep-center/distinct_corners/iter1/depletion_results.h5',
          'inf_lat_dep/dep-center/distinct_corners/iter2/depletion_results.h5',
          'inf_lat_dep/dep-center/distinct_corners/iter3/depletion_results.h5',
          'inf_lat_dep/dep-center/distinct_corners/iter4/depletion_results.h5']
mat_ids1 = ['1', '13', '13', '13', '13']
stepcolors= ['turquoise', 'teal', 'orchid', 'purple', 'firebrick']

check_converge("Convergence Check: Passwise Corners",
               'Concentration', 'Time', 'passwise-converge', fnames1, mat_ids1,
                   stepcolors)

fnames2 = ['inf_lat_dep/iter0/depletion_results.h5', 
          'inf_lat_dep/dep-center/core_avg/iter1/depletion_results.h5',
          'inf_lat_dep/dep-center/core_avg/iter2/depletion_results.h5',
          'inf_lat_dep/dep-center/core_avg/iter3/depletion_results.h5',
          'inf_lat_dep/dep-center/core_avg/iter4/depletion_results.h5']
mat_ids2 = ['1', '14', '14', '14', '14']

check_converge("Convergence Check: Core-Averaged Corners",
               'Concentration', 'Time', 'avg-converge', fnames2, mat_ids2,
                   stepcolors)


