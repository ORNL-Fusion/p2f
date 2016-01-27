# p2f
particle to distribution function

The two example cases have "real" Maxwellian distrubtions, i.e., more particles at the thermal velocity rather than a uniform velocity space grid and variable weight.

These sample particle lists (the "f.nc" files) are created using an [IDL script](https://github.com/dlg0/sMC/blob/smc_cpp/idl/create_test_particle_f.pro) available in the sMC repository using the following command 

```
IDL> create_test_particle_f, /standard_maxwellian_3d, energy_keV=0.5, n_particles=10000, eqdsk='eqdsk', density_m3=1e13, Z=1.0, amu=1.0
```

This will create a uniform (flat) density / uniform temperature (Maxwellian) particle list with vTh=sqrt(3kT/m)

There are 3 smoothing options that p2f provides ...


* No smoothing at all, i.e., just a histogram.
```
DistributeAlongOrbit = .false.
gParticle = .false.
```

```
lap101336:Cmod_case dg6$ mpirun -n 4 ~/code/p2f/xp2f.lap101336
...
 Time taken:   0.243772984
TookMaxStepsBeforeBounce:       0.00%
   *** this means you need a larger MaxSteps
Wall:       0.00%
Bad:        0.00%
off_vGrid:   0.000%   *** only applicable for gParticle = .false.
badWeight:  0.00%
badEnergy:  0.00%
Suggested eNorm:  46.2 keV
max ( density ): 0.20E+14
```
We can then run the IDL post processing script to look at the result
```
IDL>plot_p2f
Total number of particles :   9.05868e+12
```
![No smoothing f0](https://github.com/dlg0/p2f/blob/master/example/Cmod_case/p2f_f0-0.png)
<img src="https://github.com/dlg0/p2f/blob/master/example/Cmod_case/p2f_profiles-0.png" width="600px">

* Distribute a particle along its guiding center orbit according to how long it spends in each velocity space bin
```
DistributeAlongOrbit = .true.
gParticle = .false.
```
```
lap101336:Cmod_case dg6$ mpirun -n 4 ~/code/p2f/xp2f.lap101336
...
 Time taken:    23.9451637
TookMaxStepsBeforeBounce:       0.00%
   *** this means you need a larger MaxSteps
Wall:       4.88%
Bad:        0.00%
off_vGrid:   0.000%   *** only applicable for gParticle = .false.
badWeight:  0.00%
badEnergy:  0.00%
Suggested eNorm:  46.2 keV
 Getting CPU time
max ( density ): 0.12E+14
```
```
IDL>plot_p2f
Total number of particles :   8.61599e+12
```
The number of total particles is LOWER here by the 4.69% that were lost due to their orbits going outside the LCFS, i.e., they should never have really been in the particle list in the first place.
![No smoothing](https://github.com/dlg0/p2f/blob/master/example/Cmod_case/p2f_f0-1.png)
<img src="https://github.com/dlg0/p2f/blob/master/example/Cmod_case/p2f_profiles-1.png" width="600px">

* Distribute and use a gaussian particle shape in velocity space
```
DistributeAlongOrbit = .true.
gParticle = .true.
particleSize = 4.0e-04
```
```
lap101336:Cmod_case dg6$ mpirun -n 4 ~/code/p2f/xp2f.lap101336
...
 Time taken:    221.696213
TookMaxStepsBeforeBounce:       0.00%
   *** this means you need a larger MaxSteps
Wall:       4.88%
Bad:        0.00%
off_vGrid:   0.000%   *** only applicable for gParticle = .false.
badWeight:  0.00%
badEnergy:  0.00%
Suggested eNorm:  46.2 keV
 Getting CPU time
max ( density ): 0.12E+14
```
```
IDL>plot_p2f
Total number of particles :   8.57338e+12
```
![No smoothing](https://github.com/dlg0/p2f/blob/master/example/Cmod_case/p2f_f0-2.png)
<img src="https://github.com/dlg0/p2f/blob/master/example/Cmod_case/p2f_profiles-2.png" width="600px">
