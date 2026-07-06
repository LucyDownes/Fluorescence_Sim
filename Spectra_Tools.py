import numpy
import matplotlib.pyplot as pyplot
import os
from arc import *
pyplot.rcParams.update(pyplot.rcParamsDefault) #reset matplotlib defaults as ARC overrides them!

def make_transrate_LUT(atom_type = 'Cs', n_max = 80, l_max = 5, temp = 350, save = True, printing = False):
    '''
    Creates a list of all atomic states within specified limits and evaluates the rate of transitions between them (in s^{-1})
    using the ARC `getTransitionRate` function. Since the method in ARC will error if the transition is not dipole allowed 
    (instead of just returning zero) the function first checks this is satisfied. 
    If the transition is not dipole allowed, the returned rate will be exactly zero.

    Args:
        atom_type (str `'Cs'` or `'Rb'`): The atomic species to consider. 
            Default is Cs (caesium), Rb is rubidium 87.
        n_max (float, optional): maximum value of the principal quantum number $n$ to be considered. 
            Default is 80.
        l_max (float, optional): maximum value of the angular momentum quantum number $l$ to be considered. 
            Default is 5, corresponding to H states.
        temp (float, optional): temperature in K. Affects the rates of transitions affected by blackbody radiation. 
            Default is 350 K. 
        save (bool, optional): whether the calculated array should be saved automatically. 
            Filename is "Transition_Rates_{atom_type}_nmax={n_max}_temp={temp}K.csv". Default is `True`.
        printing (bool, optional): Enables print statements to monitor progress. 
            Default is `False`.
        recalculate (bool, optional): Whether to overwrite a previously calculated LUT (if it exists).
            Default is `False`. 
    
    Returns:
        states (2darray of floats): array of all atomic states with n<=n_max and l<=l_max in format (n,l,j)
        rates (2darray of floats): array of rates (in s^{-1}) of each possible transition between states in `states` array. 
            Will be of size s x s for s states.'''

   
    atom_dict = {'Cs':Caesium(), 'Rb':Rubidium87()}
    atom = atom_dict[atom_type]
    ## Sets the ground state, in this case the state at which the Monte-Carle model stops
    ground_dict = {'Cs':([6,0,0.5]), 'Rb':([5,0,0.5])}
    ground_state = ground_dict[atom_type]

    ## States below the ground state spectroscopically but above energetically
    states = atom.extraLevels[:]
    ## Add the ground state
    states.append(ground_state)
    ## Now add all states up to n = n_max
    for n in range(atom.groundStateN,int(n_max+1)):
        for l in range(0,int(l_max)):
            if l==0:
                new_state = [n,l,0.5]
                if new_state not in states:
                    states.append(new_state)
            else:
                new_state = [n,l,l-0.5]
                if new_state not in states:
                    states.append(new_state)
                new_state = [n,l,l+0.5]
                if new_state not in states:
                    states.append(new_state)
    if printing:
        print('Max. n: {:.0f}'.format(n_max))
        print('Total states: {}'.format(int(len(states))))

    transitionRates = numpy.zeros((len(states),len(states)))
    for i in range(len(states)):
        if printing:
            print('Calculating transitions for state {:.0f} of {:.0f}'.format(i+1, len(states)), end = '\r')
        for j in range(len(states)):
            n1,l1,j1=states[i]
            n2,l2,j2=states[j]
            if abs(l2-l1)==1:
                if abs(j2-j1)<2-1e-8:
                    transitionRates[i,j] = atom.getTransitionRate(n1,l1,j1,n2,l2,j2,temperature=temp)
    if printing:
        print('\r')
    if save:
        table = numpy.hstack((states, transitionRates))
        numpy.savetxt('Transition_Rates_{}_nmax={}_temp={}K.csv'.format(atom_type, int(n_max), int(temp)), table, delimiter = ',')
    return states, transitionRates


def get_rates_from_LUT(path = None, atom_type = 'Cs', n_max = 80, l_max = 5, temp = 350, save = True, printing = False, recalculate_LUT = False):
    '''
    Check to see whether the required table of transition rates exists in the location specified by path. If not, calculate a new one for the given parameters.
    Args:
        path (str, optional): The path to the precalculated table. 
            Default is None, so only the current working directory is checked.
        atom_type (str `'Cs'` or `'Rb'`): The atomic species to consider. 
            Default is Cs (caesium), Rb is rubidium 87.
        n_max (float, optional): maximum value of the principal quantum number $n$ to be considered. 
            Default is 80.
        l_max (float, optional): maximum value of the angular momentum quantum number $l$ to be considered. 
            Default is 5, corresponding to H states.
        temp (float, optional): temperature in K. Affects the rates of transitions affected by blackbody radiation. 
            Default is 350 K. 
        save (bool, optional): whether the calculated array should be saved automatically. 
            Filename is "Transition_Rates_{atom_type}_nmax={n_max}_temp={temp}K.csv". Default is `True`.
        printing (bool, optional): Enables print statements to monitor progress. 
            Default is `False`.
        recalculate_LUT (bool, optional): Whether to overwrite a previously calculated LUT (if it exists).
            Default is `False`. 
    
    Returns:
        states (2darray of floats): array of all atomic states with n<=n_max and l<=l_max in format (n,l,j)
        rates (2darray of floats): array of rates (in s^{-1}) of each possible transition between states in `states` array. 
            Will be of size s x s for s states.'''
    
    if path != None:
        file = path+'\\Transition_Rates_{}_nmax={:.0f}_temp={:.0f}K'.format(atom_type, n_max, temp)
    else:
        file = 'Transition_Rates_{}_nmax={:.0f}_temp={:.0f}K'.format(atom_type, n_max, temp)
    filename = file+'.csv'
    if os.path.exists(filename):
        if not recalculate_LUT:
            if printing:
                print('Using existing LUT of transition rates. If you need to overwrite the existing file, set `recalculate_LUT=True` and run again.')
            LUT = numpy.genfromtxt(filename, delimiter = ',')
            states = LUT[:,:3]
            rates = LUT[:,3:]
    else:
        if printing:
            print('Calculating new transition rate LUT. This may take some time.')
        states, rates = make_transrate_LUT(atom_type = atom_type, n_max = n_max, l_max = l_max, temp = temp, save = save, printing = printing)
    return states, rates
        

def make_transwvl_LUT(atom_type = 'Cs', n_max = 80, l_max = 5, save = True, printing = False):
    '''
    Creates a list of all atomic states within specified limits and evaluates the rate of transitions between them (in s^{-1})
    using the ARC `getTransitionRate` function. Since the method in ARC will error if the transition is not dipole allowed 
    (instead of just returning zero) the function first checks this is satisfied. 
    If the transition is not dipole allowed, the returned rate will be exactly zero.

    Args:
        atom_type (str `'Cs'` or `'Rb'`): The atomic species to consider. 
            Default is Cs (caesium), Rb is rubidium 87.
        n_max (float, optional): maximum value of the principal quantum number $n$ to be considered. 
            Default is 80.
        l_max (float, optional): maximum value of the angular momentum quantum number $l$ to be considered. 
            Default is 5, corresponding to H states.
        save (bool, optional): whether the calculated array should be saved automatically. 
            Filename is "Transition_Wavelengths_{}_nmax={n_max}.csv". 
            Default is `True`.
        printing (bool, optional): Enables print statements to monitor progress. 
            Default is `False`.
        recalculate (bool, optional): Whether to overwrite a previously calculated LUT (if it exists).
            Default is `False`. 
    
    Returns:
        states (2darray of floats): array of all atomic states with n<=n_max and l<=l_max in format (n,l,j)
        rates (2darray of floats): array of wavelengths (in nm) of each possible transition between states in `states` array. 
            Will be of size s x s for s states.'''

    atom_dict = {'Cs':Caesium(), 'Rb':Rubidium87()}
    atom = atom_dict[atom_type]
    ## Sets the ground  (lowest energy) state
    ground_dict = {'Cs':([6,0,0.5]), 'Rb':([5,0,0.5])}
    ground_state = ground_dict[atom_type]

    ## States below the ground state spectroscopically but above energetically
    states = atom.extraLevels[:]
    ## Add the ground state
    states.append(ground_state)
    ## Now add all states up to n = n_max
    for n in range(atom.groundStateN,int(n_max+1)):
        for l in range(0,int(l_max)):
            if l==0:
                new_state = [n,l,0.5]
                if new_state not in states:
                    states.append(new_state)
            else:
                new_state = [n,l,l-0.5]
                if new_state not in states:
                    states.append(new_state)
                new_state = [n,l,l+0.5]
                if new_state not in states:
                    states.append(new_state)
    if printing:
        print('Max. n: {:.0f}'.format(n_max))
        print('Total states: {}'.format(int(len(states))))

    transitionWvls = numpy.zeros((len(states),len(states)))
    for i in range(len(states)):
        if printing:
            print('Calculating wavelengths for state {:.0f} of {:.0f}'.format(i+1, len(states)), end = '\r')
        for j in range(len(states)):
            n1,l1,j1=states[i]
            n2,l2,j2=states[j]
            if abs(l2-l1)==1:
                if abs(j2-j1)<2-1e-8:
                    transitionWvls[i,j] = abs(atom.getTransitionWavelength(n1,l1,j1,n2,l2,j2))*1e9 #nm
    if printing:
        print('\r')
    if save:
        table = numpy.hstack((states, transitionWvls))
        numpy.savetxt('Transition_Wavelengths_{}_nmax={}.csv'.format(atom_type, int(n_max)), table, delimiter = ',')
    return states, transitionWvls

def get_wvls_from_LUT(path = None, atom_type = 'Cs', n_max = 80, l_max = 5, save = True, printing = False, recalculate_LUT = False):
    '''
    Check to see whether the required table of transition wavelengths exists in the location specified by path. If not, calculate a new one for the given parameters.
        Args:
        path (str, optional): The path to the precalculated table. 
            Default is None, so only the current working directory is checked.
        atom_type (str `'Cs'` or `'Rb'`): The atomic species to consider. 
            Default is Cs (caesium), Rb is rubidium 87.
        n_max (float, optional): maximum value of the principal quantum number $n$ to be considered. 
            Default is 80.
        l_max (float, optional): maximum value of the angular momentum quantum number $l$ to be considered. 
            Default is 5, corresponding to H states.
        save (bool, optional): whether the calculated array should be saved automatically. 
            Filename is "Transition_Wavelengths_{}_nmax={n_max}.csv". 
            Default is `True`.
        printing (bool, optional): Enables print statements to monitor progress. 
            Default is `False`.
        recalculate_LUT (bool, optional): Whether to overwrite a previously calculated LUT (if it exists).
            Default is `False`. 
    
    Returns:
        states (2darray of floats): array of all atomic states with n<=n_max and l<=l_max in format (n,l,j)
        rates (2darray of floats): array of wavelengths (in nm) of each possible transition between states in `states` array. 
            Will be of size s x s for s states.'''
    
    if path != None:
        filename = path+'\\Transition_Wavelengths_{}_nmax={:.0f}.csv'.format(atom_type, n_max)
    else:
        filename = 'Transition_Wavelengths_{}_nmax={:.0f}.csv'.format(atom_type, n_max)
    if os.path.exists(filename):
        if not recalculate_LUT:
            if printing:
                print('Using existing LUT of transition wavelengths. If you need to overwrite the existing file, set `recalculate_LUT=True` and run again.')
            LUT = numpy.genfromtxt(filename, delimiter = ',')
            states = LUT[:,:3]
            wvls = LUT[:,3:]
    else:
        if printing:
            print('Calculating new wavelength LUT. This may take some time.')
        states, wvls = make_transwvl_LUT(atom_type = atom_type, n_max = n_max, l_max = l_max, save = save, printing = printing)
    return states, wvls

def simulate_spectrum_mc(state, atom_type = 'Cs', temp = 350, n_max = 80, iters = 5000, 
spectrum_range = (400,750), spectrum_resolution = 0.5, give_paths = False, give_pops = False, path_to_LUT = None, recalculate_LUT = False, save_LUT = True, printing = False):
    '''
    Simulates the fluorescence from a specified atomic state via a Monte-Carlo approach
    Args:
        state (tuple): the state to evaluate the fluorescence spectrum of, 
            in the format (n,l,j).
        atom_type (string, optional): `Rb` or `Cs`. Defaults to Cs (caesium), 
            Rb is Rubidium 87.
        temp (float, optional): temperature of the atomic ensemble in Kelvin. 
            *Must* match the temperature used to create the look-up-table. Default is 350 K.
        n_max (float, optional): maximum value of n to be considered. 
            *Must* match the value of n_max used to create the look-up-table. Default is 80.
        iters (float, optional): number of iterations to run the simulation for. 
            Default is 5000 (5e3) which balances convergence with simulation time.
        spectrum_range (tuple, optional): the wavelength range (in nm) over which 
            the spectrum will be evaluated. Any occurrences of wavelengths outside 
            of this range will not be saved. Defaults to (400, 750), corresponding to visible light.
        spectrum_resolution (float, optional): bin width (in nm) used to create 
            the histogram. Smaller values will yield more bins (hence a higher 
            resolution spectrum) but data files will be larger. Default is 0.5.
        give_paths (bool, optional): whether to return information on decay pathways 
            alongside the spectrum. Default is False. 
        give_pops (bool, optional): whether to return information on the frequency 
            with which each state was visited (population passing through each state). 
            Default is False. 
        path_to_LUT (string, optional): file path of the "Transition_Rates" 
            look-up table. Will be used to load the states and rates arrays. 
            Default is `None` (assumes the file is in the current working directory)
        recalculate_LUT (bool, optional): whether to recalculate LUTs even if they already exist
            Default is False.
        save_LUT (bool, optional): whether to save any LUTs calculated to the current working directory.
            Default is True.
        printing (bool, optional): whether to print status messages about code progress.
            Default is False.
 
    Returns:
        bin_edges (1darray): the right-hand edges of the histogram bins. 
            Essentially the wavelengths of the simulated emissions (in nm).
        hist/iters (1darray): the occurrences of each wavelength in `bin_edges`, 
            divided by the total number of iterations performed. Note that the 
            sum of this array will likely be greater than the number of iterations, 
            as more than one photon is likely emitted per decay.
 
    If give_pops, also returns:
        all_states (2darray): array of every state visited at least once during 
            simulation which resulted in emission of a photon within the specified range. 
            States will be in the form (n,l,j). Note that the states are 
            returned in no particular order.
        states_count (1darray): number of times each of the states in `all_states` 
            was visited.
     
    If give_paths, also returns:
        all_trans (2darray): array of every transition that resulted in an emission 
            of a photon in the specified range. In the format (n_1, l_1, j_1, n_2, l_2, j_2) 
            for transition between states 1 and 2. Note that the transitions are 
            returned in no particular order.
        trans_count (1darray): number of times each transition in `all_trans` 
            occurred in the simulation.'''
    
    (n1, l1, j1) = state
    states, transitionRates = get_rates_from_LUT(path_to_LUT, n_max = n_max, atom_type = atom_type, temp = temp, recalculate_LUT = recalculate_LUT, save=save_LUT, printing = printing)
    spdf = ['S', 'P', 'D', 'F', 'G', 'H']
    atom_dict = {'Cs':Caesium(), 'Rb87':Rubidium87()}
    atom = atom_dict[atom_type]
    ground_dict = {'Cs':([6,0,0.5]), 'Rb87':([5,0,0.5])}
    ground_state = ground_dict[atom_type]
    starting_index = numpy.where(numpy.isclose(states,[n1,l1,j1]).all(axis=1))[0][0]

    emitted_wavelengths = []
    pathways = []
    trans = []
    popn = []
    for i in range(int(iters)):
        state_index = starting_index
        while True:
            n1,l1,j1 = states[state_index]
            rates = transitionRates[state_index,:] # extract transition rates of all transitions out of current state
            probs = rates/sum(rates) # normalise
            probs_cum = numpy.cumsum(probs) # create cumulative probability array
            new_state_index = numpy.where(probs_cum > numpy.random.random())[0][0] # choose state to move to based on RNG
            n2,l2,j2 = states[new_state_index] # add new state 
            wvl = -atom.getTransitionWavelength(int(n1),int(l1),j1,int(n2),int(l2),j2) # calculate wavelength of photon emitted in transition
            if spectrum_range[0]<=wvl*1e9<=spectrum_range[1]:
                emitted_wavelengths.append(wvl) # if photon is in specified range, store
                pathway = '{}{}{}/2 -> {}{}{}/2'.format(int(n1),spdf[int(l1)],int(j1*2),int(n2),spdf[int(l2)],int(j2*2))
                pathways.append(pathway) # add details of transition pathway to array
                trans.append([n1,l1,j1,n2,l2,j2]) # add state pair 
                popn.append([int(n1),int(l1),j1]) # add state to track population
            state_index = new_state_index
            if (n2<=ground_state[0]) & (l2==ground_state[1]):
                popn.append([int(n2),int(l2),j2])
                break
                
    emitted_wavelengths = numpy.array(emitted_wavelengths)*1e9 # turn list of wavelengths into array and convert to nm
    hist,bin_edges = numpy.histogram(emitted_wavelengths,int((spectrum_range[1] - spectrum_range[0])/spectrum_resolution),spectrum_range)
    hist = numpy.array(hist,dtype=float) # photon 'counts' for each wavelength bin
    pops = numpy.asarray(popn)
    all_states, states_count = numpy.unique(popn, axis = 0, return_counts=True) # unique states visited, and number of times visited 
    all_trans, trans_counts = numpy.unique(trans, axis = 0, return_counts=True) # unique transitions that resulted in photon emission, and number of times the transition occurred
    
    if give_pops and give_paths:
        return bin_edges[:-1], hist/iters, all_states, states_count, all_trans, trans_counts
    elif give_pops:
        return bin_edges[:-1], hist/iters, all_states, states_count
    elif give_paths:
        return bin_edges[:-1], hist/iters, all_trans, trans_counts
    else:
        return bin_edges[:-1], hist/iters


def steady_state_popn(excitation_states, beam_intensities, atom_type = 'Cs', temp = 350, n_max = 80, path_to_LUT=None, recalculate_LUT = False, printing = False):
    '''
    Calculates the steady-state population fraction in each state for a specified excitation scheme and laser intensities
    Args:
        excitation_states (list of tuples): the states involved in the driving scheme in the format (n,l,j)
        beam_intensities (list of floats): intensities of the fields coupling the states in the driving scheme. 
            *Must* be in the same order as the specified excitation states, e.g. beam_intensities[i] 
            is the intensity of the field coupling excitation_states[i] to excitation_states[i+1].
        atom_type (string, optional): `Rb` or `Cs`. Defaults to Cs (caesium), 
            Rb is Rubidium 87.
        temp (float, optional): temperature of the atomic ensemble in Kelvin. 
            *Must* match the temperature used to create the look-up-table. Default is 350 K.
        n_max (float, optional): maximum value of n to be considered. 
            *Must* match the value of n_max used to create the look-up-table. Default is 80.
        path_to_LUT (string, optional): file path of the "Transition_Rates" 
            look-up table. Will be used to load the states and rates arrays. 
            Default is `None` (assumes the file is in the current working directory)
        recalculate_LUT (bool, optional): whether to recalculate LUTs even if they already exist
            Default is False.
        save_LUT (bool, optional): whether to save any LUTs calculated to the current working directory.
            Default is True.
        printing (bool, optional): whether to print status messages about code progress.
            Default is False.
         
    Returns:
        states (2darray): all of the states considered in the rate equation model in the format (n,l,j). 
            Will be set by the LUT for the transition rates.
        norm_Ns (1darray): the steady-state population fraction in each of the states returned by states.'''
    
    states, A_coeffs = get_rates_from_LUT(path_to_LUT, n_max = n_max, atom_type = atom_type, temp = temp, recalculate_LUT = recalculate_LUT, printing = printing)
    n_levels = len(states)

    M_spont = numpy.zeros((n_levels, n_levels))
    for i in range(n_levels):
        for j in range(n_levels):
            if i==j:
                M_spont[i,j] = -numpy.sum(A_coeffs[i,:])
            else:
                M_spont[i,j] = A_coeffs[j,i]
    
    indices = numpy.zeros(len(excitation_states))
    for i in range(len(excitation_states)):
        state = excitation_states[i]
        indices[i] = numpy.where(numpy.isclose(states, state).all(axis = 1))[0][0]
    M_stim = numpy.zeros(M_spont.shape)
    for i in range(len(beam_intensities)):
        i1, i2 = int(indices[i]), int(indices[i+1])
        W = beam_intensities[i]*A_coeffs[i2, i1]
        M_stim[i1, i1] += -W
        M_stim[i2, i2] += -W
        M_stim[i1, i2] += W
        M_stim[i2, i1] += W
    M_tot = M_spont + M_stim

    u, sig, v = numpy.linalg.svd(M_tot)
    abs_sig = numpy.abs(sig)
    minval, maxval = numpy.min(abs_sig), numpy.max(abs_sig)
    if minval>maxval*1e-15:
        # if there is no zero singular value (within floating point tolerance)
        # return an empty array
        print('ERROR - Matrix is non-singular')
        
    index = abs_sig.argmin() # find the position of the zero solution
    Ns = numpy.conjugate(v[index,:]) # extract the solution, which is the column of V at the position of the zero singular value
    norm_Ns = Ns/numpy.sum(Ns)
    return states, norm_Ns

def spectrum_from_pops(states, rates, trans_wvls, pop_frac, spectrum_range = (400,750), spectrum_resolution = 0.5):
    '''
    Creates simulated emission spectrum for a given list of atomic states and the steady-state population fraction in each.
    Returns the emitted photons per second per atom in the wavelength range specified by `spectrum_range` and `spectrum_resolution`.

    Args:
        states (2darray): n by 3 array containing all atomic states to consider, in the format n,l,j.
        rates (2darray): array of transition rates (in s^{-1}) of each possible transition between states in `states` array.
            *Must* be of shape n x n for n states.
        trans_wvls (2darray): array of transition wavelengths (in nm) of each possible transition between states in `states` array.
            *Must* be the same shape as `rates`.
        pop_frac (1darray): array containing the fractional steady-state population in each of the states in `states`.
            *Must* be the same length as `states`.
        spectrum_range (tuple, optional): the wavelength range (in nm) over which 
            the spectrum will be evaluated. Any occurrences of wavelengths outside 
            of this range will not be saved. Defaults to (400, 750), corresponding to visible light.
        spectrum_resolution (float, optional): step size (in nm) used to create 
            the spectrum. Smaller values will yield more bins (hence a higher 
            resolution spectrum) but the calculation will take longer. Default is 0.5.
    
    Returns:
        spectrum_wvls (1darray of floats): the regularly spaced array of wavelengths over which the spectrum has been calculated.
        spectrum (1darray of floats): the emission rate in photons per second per atom for each of the wavelength intervals in `spectrum_wvls`.'''
    
    n_levels = len(states)
    weights = numpy.zeros(rates.shape)
    for i in range(n_levels):
        # Calculate the `weight` of each emission line as the initial state population multiplied by the transition rate to the final state
        weights[:,i] = rates[:,i]*pop_frac
    
    spectrum_wvls = numpy.arange(spectrum_range[0], spectrum_range[1], spectrum_resolution)
    spectrum = numpy.zeros(spectrum_wvls.shape)
    for i in range(len(spectrum_wvls)):
        # Find all transitions that result in photons in the given wavelength range and store their locations
        all_trans = numpy.where(numpy.logical_and(trans_wvls>spectrum_wvls[i], trans_wvls<=spectrum_wvls[i]+spectrum_resolution))
        xs, ys = all_trans[0], all_trans[1]
        for j in range(len(xs)):
            # Add the weights of these transitions to the spectrum array
            spectrum[i] += weights[xs[j], ys[j]]
    return spectrum_wvls, spectrum

def simulate_spectrum_re(excitation_states, beam_intensities, spectrum_range = (400,750), spectrum_resolution = 0.5, atom_type = 'Cs', temp = 350, n_max = 80, 
                         path_to_LUTs = None, recalculate_LUTs = False, save_LUTs = True, printing = False):
    '''
    Simulates the emission spectrum from an atom using a rate equation approach. Requires the states used in the driving scheme and the beam powers in units of I_sat.

    Args:
        excitation_states (2darray): n by 3 array containing all atomic states involved in the driving scheme, in the format n,l,j.
        beam_intensities (1darray): the intensities of the driving fields in units of the saturation intensity. 
            *Must* be of length len(excitation_states) - 1
        spectrum_range (tuple, optional): the wavelength range (in nm) over which 
            the spectrum will be evaluated. Any occurrences of wavelengths outside 
            of this range will not be saved. Defaults to (400, 750), corresponding to visible light.
        spectrum_resolution (float, optional): step size (in nm) used to create 
            the spectrum. Smaller values will yield more bins (hence a higher 
            resolution spectrum) but the calculation will take longer. Default is 0.5.
        atom_type (string, optional): `Rb` or `Cs`. Defaults to Cs (caesium), 
            Rb is Rubidium 87.
        temp (float, optional): temperature of the atomic ensemble in Kelvin. 
            *Must* match the temperature used to create the look-up-tables. Default is 350 K.
        n_max (float, optional): maximum value of n to be considered. 
            *Must* match the value of n_max used to create the look-up-tables. Default is 80.
        path_to_LUTs (string, optional): file path of the "Transition_Rates" and "Transition Wavelengths" 
            look-up tables. Will be used to load the states and rates arrays. 
            Default is `None` (assumes the file is in the current working directory)
        recalculate_LUTs (bool, optional): whether to recalculate LUTs even if they already exist
            Default is False.
        save_LUTs (bool, optional): whether to save any LUTs calculated to the current working directory.
            Default is True.
        printing (bool, optional): whether to print status messages about code progress.
            Default is False.
    
    Returns:
        spectrum_wvls (1darray of floats): the regularly spaced array of wavelengths over which the spectrum has been calculated.
        spectrum (1darray of floats): the emission rate in photons per second per atom for each of the wavelength intervals in `spectrum_wvls`.'''
    
    states2, transwvls = get_wvls_from_LUT(path=path_to_LUTs, n_max = n_max, atom_type = atom_type, recalculate_LUT = recalculate_LUTs, save = save_LUTs, printing = printing)
    states3, transrates = get_rates_from_LUT(path=path_to_LUTs, n_max = n_max, atom_type = atom_type, temp = temp, recalculate_LUT = recalculate_LUTs, save = save_LUTs, printing = printing)
    states, norm_Ns = steady_state_popn(excitation_states, beam_intensities, path_to_LUT=path_to_LUTs, n_max = n_max, atom_type = atom_type, temp = temp, recalculate_LUT = False, printing=False)
    spectrum_wvls, spectrum = spectrum_from_pops(states, transrates, transwvls, norm_Ns, spectrum_range=spectrum_range, spectrum_resolution=spectrum_resolution)
    return spectrum_wvls, spectrum

    
