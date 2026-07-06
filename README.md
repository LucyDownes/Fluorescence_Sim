# Rydberg Fluorescence Simulation

Many Rydberg states emit fluorescence in the visible, making them ideal for visualising external fields for example THz fields [[1,2,3]](#References). In order to inform the choice of atom or atomic transition it is helpful to be able to predict the spectrum of fluorescence emitted by a specific Rydberg state. The functions in the `Spectra_Tools` module enable this to be done in two different ways: either through a Monte-Carlo simulation of decay from a target Rydberg state, or by solving the rate equations to find the steady-state population and using this to build the emission spectrum.

All of these codes are built on ARC [(Alkali Rydberg Calculator)](https://arc-alkali-rydberg-calculator.readthedocs.io/en/latest/), and I cannot be held responsible for any errors therein. The functions are described briefly below, more detailed explanations, benchmarking and examples of their use can be found in the accompanying Jupyter Notebooks.

## Creating Look-Up Tables (LUTs)

Both methods rely on some pre-calculated look-up tables (LUTs) of transition rates and transition wavelengths between all possible pairs of atomic states in a given range (maximum $n$ and $l$ values).
If the files are not found the code will generate them automatically using the default values ($n_{\rm{max}}=80$, $l_{\rm{max}} = 5$, temperature = $350\,\rm{K}$). Note that it can take some time to run (several minutes) depending on parameters.
The higher the principal quantum number that is considered the more accurate the simulation will be, however this comes at an increase in time to create the LUTs. The bare minimum needed would be to set `n_max` to be a little higher than the principal quantum number of the state you are evaluating the spectrum for. So if you are driving to the $14\rm{D}_{3/2}$ state, setting `n_max = 16` would give you a reasonable estimate. However, since this table only needs to be evaluated once, it is worth going to higher $n$ (80 or more) to save time overall. 

### Transition Rates
The transition rates look-up table is a list of all atomic states within the specified limits and the rate of transitions between them (in $s^{-1}$) found using the ARC `getTransitionRate` function. Since the method in ARC will error if the transition is not dipole allowed (instead of just returning zero) the function first checks this is satisfied. If the transition is not dipole allowed, the returned rate will be exactly zero.  

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

### Transition Wavelengths
The transition wavelengths look-up table is a list of all atomic states within the specified limits and the wavelength of transitions between them (in nm) using the ARC `getTransitionWavelength` function. Since the method in ARC will error if the transition is not dipole allowed (instead of just returning zero) the function first checks this is satisfied. If the transition is not dipole allowed, the returned wavelength will be exactly zero. Since ARC returns negative wavelength values for transitions from a higher energy state to a lower energy one, this function returns the absolute value.

Again, the order in which the wavelengths are given corresponds to the order that the states are in. For example for 3 states $1,2,3$, the array will be

$$\begin{pmatrix}
\lambda_{11} & \lambda_{12} & \lambda_{13} \\
\lambda_{21} & \lambda_{22} & \lambda_{23} \\
\lambda_{31} & \lambda_{32} & \lambda_{33}
\end{pmatrix}$$

where $\lambda_{ij}$ is the wavelength of the transition $i\rightarrow j$. The array should always be antidiagonal, since $\lambda_{ii} = 0$. The file will be saved in the format 

$$\begin{pmatrix}
n_1 l_1 j_1 & \lambda_{11} & \lambda_{12} & \lambda_{13} \\
n_2 l_2 j_2 & \lambda_{21} & \lambda_{22} & \lambda_{23} \\
n_3 l_3 j_3 & \lambda_{31} & \lambda_{32} & \lambda_{33}
\end{pmatrix}$$


## Simulating Spectra

Once the look-up-tables have been created, the spectrum can be simulated using either the Monte-Carlo (MC) or rate equation (RE) approach. This is done via the `simulate_spectrum_mc` or `simulate_spectrum_re` functions respectively. The different methods require slightly different input parameters: the Monte-Carlo method requires a single input state in the form $(n,l,j)$ which is the target state in the excitation. The rate-equation method requires all states involved in the excitation scheme (assuming a ladder scheme) and the corresponding intensities of the driving beams (as a fraction of the saturation intensity). Both functions return an array of evenly spaced wavelengths and the strength of the expected emission for each of the wavelength steps.


### Notes

These methods do not include any collisional effects that occur within the vapour which we do see in experiments, for example we see evidence of population transfer between fine-structure states. The simulations can be made to more closely match the experimental results by 'manually' adding in decay from other fine-structure states, as illustrated in Fig. 4 of [[2]](#References). For example if we were wanting to simulate the spectrum for the $13\mathrm{D_{5/2}}$ state there would be some collisional transfer to the $13\rm{D_{3/2}}$ state, so we would need to simulate the spectrum from both of these states and then sum them in an appropriate ratio (usually around 80% target state, 20% other fine-structure state, but this varies a lot!).

## References

[1] Lucy A. Downes et al., *Full-Field Terahertz Imaging at Kilohertz Frame Rates Using Atomic Vapor* Phys. Rev. X **10** 011027 (2020) [https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011027](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.011027)

[2] Lucy A. Downes, Lara Torralbo-Campo and Kevin J. Weatherill, *A practical guide to terahertz imaging using thermal atomic vapour* New J. Phys. **25** 035002 (2023) [https://iopscience.iop.org/article/10.1088/1367-2630/acb80c/meta](https://iopscience.iop.org/article/10.1088/1367-2630/acb80c/meta)

[3] Lucy A. Downes, *A High-speed THz Imaging System based on THz-to-optical Conversion in Atomic Vapour* PhD Thesis, Durham University (2020) [https://etheses.dur.ac.uk/13797/](https://etheses.dur.ac.uk/13797/)
