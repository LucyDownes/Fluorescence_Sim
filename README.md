# Rydberg Fluorescence Simulation

Many Rydberg states emit fluorescence in the visible, making them ideal for visualising external fields for example THz fields [[1,2,3]](#References). In order to inform the choice of atom or atomic transition it is helpful to be able to predict the spectrum of fluorescence emitted by a specific Rydberg state. The functions in the `Spectra_Tools` module enable a Monte-Carlo simulation of population randomly decaying out of the specified state to be performed. At every step any emitted photon that is within the specified range of interest is recorded and a spectrum returned. This can be quite time consuming and requires some set-up beforehand.

All of these codes are built on ARC [(Alkali Rydberg Calculator)](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/), and I cannot be held responsible for any errors therein. The functions are described below, examples of their use can be found in the Jupyter Notebook.

## Creating Look-Up Tables (LUTs)

To decide where the decays end up based on random number generation, we need to know the probabilities/rates of each possible decay pathway. To this end we create a look-up table that contains the rate of decay from every state to every other state (within limits of maximum $n$ and $l$ values). This is done via the `make_transrate_LUT` function. This file *must* exist in order to run any Monte-Carlo fluorescence models. Note that it can take some time to run (several minutes) depending on parameters.

The function creates a list of all atomic states within the specified limits and evaluates the rate of transitions between them (in $s^{-1}$) using the ARC `getTransitionRate` function. Since the method in ARC will error if the transition is not dipole allowed (instead of just returning zero) the function first checks this is satisfied. If the transition is not dipole allowed, the returned rate will be exactly zero.  

Args:
 - `atom_type` (str `'Cs'` or `'Rb'`): The atomic species to consider. Default is Cs (caesium), Rb is rubidium 87.
 - `n_max` (float, optional): maximum value of the principal quantum number $n$ to be considered. Default is 80.
 - `l_max` (float, optional): maximum value of the angular momentum quantum number $l$ to be considered. Default is 5, corresponding to $\rm{H}$ states.
 - `temp` (float, optional): temperature in K. Affects the rates of transitions affected by blackbody radiation. Default is 350 K. 
 - `save` (bool, optional): whether the calculated array should be saved automatically. Filename is "Transition_Rates_nmax={n_max}\_temp={temp}K_{atom_type}.csv". Default is `True`.
 - `printing` (bool, optional): Enables print statements to monitor progress. Default is `False`.
 
Returns:
 - `states` (2darray of floats): array of all atomic states with $n\leq n_{\rm{max}}$ and $l\leq l_{\rm{max}}$ in format ($n,l,j$)
 - `rates` (2darray of floats): array of rates (in $s^{-1}$) of each possible transition between states in `states` array. Will be of size $s\times s$ for $s$ states. 

The order in which the rates are given corresponds to the order that the states are in. For example for 3 states $1,2,3$, the rate array will be
$$\begin{pmatrix}
\Gamma_{11} & \Gamma_{12} & \Gamma_{13} \\
\Gamma_{21} & \Gamma_{22} & \Gamma_{23} \\
\Gamma_{31} & \Gamma_{32} & \Gamma_{33}
\end{pmatrix}$$
where $\Gamma_{ij}$ is the rate of the transition $i\rightarrow j$. The array should always be antidiagonal, since $\Gamma_{ii} = 0$. The file will be saved in the format 
$$\begin{pmatrix}
n_1 l_1 j_1 & \Gamma_{11} & \Gamma_{12} & \Gamma_{13} \\
n_2 l_2 j_2 & \Gamma_{21} & \Gamma_{22} & \Gamma_{23} \\
n_3 l_3 j_3 & \Gamma_{31} & \Gamma_{32} & \Gamma_{33}
\end{pmatrix}$$

The higher the principal quantum number that is considered the more accurate the simulation will be, however this comes at an increase in time to create the table. The bare minimum needed would be to set `n_max` to be a little higher than the principal quantum number of the state you are evaluating the spectrum for. So if you were wanting to know the spectrum of fluorescence emitted by the $14\rm{D}_{3/2}$ state, setting `n_max = 16` would give you a reasonable estimate. However, since this table only needs to be evaluated once, it is worth going to higher $n$ (80 or more) to save time overall. 

## Simulating Spectra

Once the look-up-table of transition rates has been created, the spectrum can be simulated using the Monte-Carlo approach. This is done via the `calculate_spectrum` function.

Args:
 - `state` (tuple): the state to evaluate the fluorescence spectrum of, in the format $(n,l,j)$.
 - `path_to_LUT` (string, optional): file path of the "Transition_Rates" look-up table. Will be used to load the states and rates arrays. Default is `None` (assumes the file is in the current working directory)
 - `atom_type` (string, optional): `Rb` or `Cs`. Defaults to Cs (caesium), Rb is Rubidium 87.
 - `temp` (float, optional): temperature of the atomic ensemble in Kelvin. *Must* match the temperature used to create the look-up-table. Defaults to 350 K.
 - `nmax` (float, optional): maximum value of n to be considered. *Must* match the value of n_max used to create the look-up-table. Default is 80.
 - `iters` (float, optional): number of iterations to run the simulation for. Default is 50,000 (5e4).
 - `spectrum_range` (tuple, optional): the wavelength range (in nm) over which the spectrum will be evaluated. Any occurrences of wavelengths outside of this range will not be saved. Defaults to (400, 750), corresponding to visible light.
 - `spectrum_resolution` (float, optional): bin width (in nm) used to create the histogram. Smaller values will yield more bins (hence a higher resolution spectrum) but data files will be larger. Default is 0.5.
 - `give_paths` (bool, optional): whether to return information on decay pathways alongside the spectrum. Default is False. Note that the pathways are in no particular order.
 - `give_pops` (bool, optional): whether to return information on the frequency with which each state was visited (population passing through each state). Default is False. 
 
Returns:
 - `bin_edges` (1darray): the right-hand edges of the histogram bins. Essentially the wavelengths of the simulated emissions.
 - `hist/iters` (1darray): the occurrences of each wavelength in `bin_edges`, divided by the total number of iterations performed. Hence if a photon was emitted on each iteration, the value would be 1. Note that the sum of this array will likely be greater than the number of iterations, as more than one photon is emitted per decay.
 
If `give_pops = True`, also returns:
 - `all_states` (2darray): array of every state visited at least once during simulation which resulted in emission of a photon within the specified range. States will be in the form $(n,l,j)$.
 - `states_count` (1darray): number of times each of the states in `all_states` was visited.
     
If `give_paths = True`, also returns:
 - `all_trans` (2darray): array of every transition that resulted in an emission of a photon in the specified range. In the format $(n_1, l_1, j_1, n_2, l_2, j_2)$ for transition between states 1 and 2.
 - `trans_count` (1darray): number of times each transition in `all_trans` occurred in the simulation.

### Notes

This method does not include any effects that the applied driving fields may have on the fluorescence, for example optical pumping. It also does not take into account any collisional effects that occur within the vapour which we do see in the experiment, for example we see evidence of population transfer between fine-structure states. The simulations can be made to more closely match the experimental results by 'manually' adding in decay from other fine-structure states, as illustrated in Fig. 4 of [[2]](#References). For example if we were wanting to simulate the spectrum for the $13\rm{D}_{5/2}$ state there would be some collisional transfer to the $13\rm{D}_{3/2}$ state, so we would need to simulate the spectrum from both of these states and then sum them in an appropriate ratio (usually around 80% target state, 20% other fine-structure state, but this varies a lot!).

## References

[1] Lucy A. Downes et al., *Full-Field Terahertz Imaging at Kilohertz Frame Rates Using Atomic Vapor* Phys. Rev. X **10** 011027 (2020) [https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011027](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011027)

[2] Lucy A. Downes, Lara Torralbo-Campo and Kevin J. Weatherill, *A practical guide to terahertz imaging using thermal atomic vapour* New J. Phys. **25** 035002 (2023) [https://iopscience.iop.org/article/10.1088/1367-2630/acb80c/meta](https://iopscience.iop.org/article/10.1088/1367-2630/acb80c/meta)

[3] Lucy A. Downes, *A High-speed THz Imaging System based on THz-to-optical Conversion in Atomic Vapour* PhD Thesis, Durham University (2020) [https://etheses.dur.ac.uk/13797/](https://etheses.dur.ac.uk/13797/)
