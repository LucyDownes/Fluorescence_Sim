import numpy
import matplotlib.pyplot as pyplot
from arc import *
pyplot.rcParams.update(pyplot.rcParamsDefault) #reset matplotlib defaults as ARC overrides them!

def make_transrate_LUT(atom_type, n_max = 80, l_max = 5, temp = 350, save = True, printing = False):
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
            Filename is "Transition_Rates_nmax={n_max}\_temp={temp}K_{atom_type}.csv". Default is `True`.
        printing (bool, optional): Enables print statements to monitor progress. Default is `False`.
    
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
    ## Now add all states up to n = nmax
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

    transitionRates = np.zeros((len(states),len(states)))
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
        numpy.savetxt('Transition_Rates_nmax={}_temp={}K_{}.csv'.format(int(n_max), int(temp), atom_type), table, delimiter = ',')
    return states, transitionRates

def get_rates_from_LUT(path=None, nmax = 100, atom_type = 'Cs', temp = 350):
    '''Loads the pre-calculated look-up table of transition rates from file for 
        use in the spectrum simulation'''
    if path != None:
        file = path+'\\Transition_Rates_nmax={:.0f}_temp={:.0f}K_'.format(nmax, temp)
    else:
        file = 'Transition_Rates_nmax={:.0f}_temp={:.0f}K_'.format(nmax, temp)
    filename = file+atom_type+'.csv'
    LUT = numpy.genfromtxt(filename, delimiter = ',')
    states = LUT[:,:3]
    rates = LUT[:,3:]
    return states, rates

def calculate_spectrum(state, path_to_LUT = None, atom_type = 'Cs', temp = 350, nmax = 80, iters = 50000, 
spectrum_range = (400,750), spectrum_resolution = 0.5, give_paths = False, give_pops = False):
    '''
    Simulates the fluorescence from a specified atomic state via a Monte-Carlo approach
    Args:
        state (tuple): the state to evaluate the fluorescence spectrum of, 
            in the format (n,l,j).
        path_to_LUT (string, optional): file path of the "Transition_Rates" 
            look-up table. Will be used to load the states and rates arrays. 
            Default is `None` (assumes the file is in the current working directory)
        atom_type (string, optional): `Rb` or `Cs`. Defaults to Cs (caesium), 
            Rb is Rubidium 87.
        temp (float, optional): temperature of the atomic ensemble in Kelvin. 
            *Must* match the temperature used to create the look-up-table. Default is 350 K.
        nmax (float, optional): maximum value of n to be considered. 
            *Must* match the value of n_max used to create the look-up-table. Default is 80.
        iters (float, optional): number of iterations to run the simulation for. 
            Default is 50,000 (5e4).
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
    states, transitionRates = get_rates_from_LUT(path_to_LUT, nmax = nmax, atom_type = atom_type, temp = temp)
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

