Modelling of the Main Sequence (Thorne-Zytkow Group)
====================================================

This code solve the physical properties for stars based on initial values (_T~c~_). Each solved star is stored as a _pickle_ file is the "_stars_" directory (a new directory must be made before solving for a set of stars). The main sequence/main sequences can then be plotted by running the _main_sequence()_ function.

Furthermore each star can be plotted for their individual properties: density, temperature, mass, etc... This can only be done after the star is solved.

Setup
-----

The range of initial value can be alternated with the _rng_T_ variable in "main.py". The parameters for the range are (starting T, ending T, number of stars), respectively. Lastly, ensure the target directory (stars) exists.

For the Thorne-Zytkow simulations, set up new variable: _M_bh_ and add it at the end of the call to the "_bisection_" method in _solveStar_. All other steps can be run the same way as with the Main Sequence. 

Building the Main Sequence
--------------------------

### With Spyder or other IDE

1. Set up initial conditions
2. Run "main.py" script
3. Call _solveAllStars_ with "_rng_T_" as parameter
4. Once stars are calculated change directory name to "MS_stars" for main sequence or "stars(1e-\_)" for varying mass conditions (Thorne-Zytkow simulations)
5. Main sequences can be individually plotted with "_main_sequence()_" or overplotted with "_main_sequence()_" (First parameter set to True to plot normal MS as well)
6. To plot multiple MS's change the array: _seqs_ to the corresponding exponents

##### Example

```sh
After running script...
>>> solveAllStars(rng_t)
>>> main_sequence()
>>> main_sequence(True)
```

### With python3 only

1. Set up initial conditions
2. Add desired function alls at the bottom of "main.py" script

Plotting Individual Stars
-------------------------

A star's property can be plotted after it has been solve. This uses the _plotAll()_ function. The star can be loaded using the _load_star()_ (takes the directory name and filename). For example,

```sh
>>> s = load_star("star(1e2)","30000000.0.p")
>>> plotAll(s)
```
