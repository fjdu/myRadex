**What's new**

- 2024-11-18: `Cython`, instead of `f2py`, is used to generate the python wrapper.  For this to work you will need to have a recent version of `Cython`, which is likely already installed if you use `anaconda`.  In case not please follow the [Installing Cython](https://cython.readthedocs.io/en/stable/src/quickstart/install.html).  Examples include a [jupyter notebook](2024-11-18-myradex-with-cython-example.ipynb) and two python scripts [example-1.py](example-1.py) and [example-2.py](example-2.py). The `f2py` version may still work, but sometimes it is hard to compile successfully for some unknown reasons.

- 2023-03-30: Now the `collisioPartnerCrit` parameter can be used to specify which collision partner (1-based) will have critical densities saved in the output.

- 2023-03-29: The critical densities (calculated in two manners) are included in the output file.  Note that different collisional partners have different sets of critical densities.  At present only the values for the first collisional partner are included in the output file.

- 2022-06-06: `myRadex` is now included in the Astrophysics Source Code Library as ["myRadex: Radex with a twist"](https://ascl.net/2205.011).

- 2021-11-25: Instead of being calculated from the energy levels, now the frequencies in the input file will be used by default.  Also the energy level numbers will not be subtracted by 1 by default.  For backward compatibility, two Boolean options are added: `recalculateFreqWithEupElow` and `iLevel_subtract_one` (both are `False` by default).
- 2020-04-05: A function for calculating critical densities of all the transitions of a molecule is included.

---

## Meaning of the output columns

- `iup`: index of upper level
- `ilow`: index of lower level
- `Eup`: energy of upper level
- `freq`: transition frequency
- `lam`: transition wavelength
- `Tex`: excitation temperature $$-\frac{E_{u} - E_{l}}{\ln\left((y_{u}g_{l})/(y_{l}g_{u})\right)}$$
- `tau`: optical depth $$\tau = \alpha \times L,$$ where $$\alpha = \frac{h\nu}{4\pi\delta\nu} n_\text{mol} (y_{l}B_{lu} - y_{u}B_{ul}),$$ and $L$ is the relevant length scale.  The frequency width is calculated from the velocity width as $$\delta\nu = \frac{\nu}{c} \delta v \times \sqrt{\pi/(4\ln2)}.$$ Where does the strange factor in the above equation come from?  For a line profile $$\phi(\nu) = \frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(\nu-\nu_0)^2}{2\sigma^2}},$$ its full-width-at-half-maximum (FWHM) is $$\text{FWHM} = 2\sqrt{2\ln2}\sigma,$$ thus the term $\sqrt{2\pi}\sigma$ can be written as $$\sqrt{2\pi}\sigma \equiv \sqrt{2\pi}\frac{\text{FWHM}}{2\sqrt{2\ln2}} = \text{FWHM}\sqrt{\frac{\pi}{4\ln2}}.$$ From this equation we see that **$\delta v$ is to be interpreted as the FWHM velocity width** (which is commonly reported in radioastronomy observations), while $\delta\nu$ is the $\sqrt{2\pi}\sigma$ term in the denominator of the Gaussian line profile function, which is more for coding convenience.  With this interpretation in mind, $\alpha$ is the peak of the absorption coefficient, and $\tau$ is the peak optical depth.
- `Tr`: $$(I - I_{\text{bg}}) \frac{c^2}{2k_\text{B}\nu^{2}},$$ where $$I= (1-e^{-\tau}) B_\nu(T_\text{ex}) + e^{-\tau} I_{\text{bg}}.$$ Note that the $\tau$ used here is the peak optical depth.
- `fup`: occupation fraction of upper level
- `flow`: occupation fraction of lower level
- `flux_K`: $${\text{Tr}} \times \delta v\times\sqrt{\pi/(4\ln2)}.$$ Why the strange factor again?  Assuming being **optically thin**, $$I_{\nu}= I_{\nu,\text{peak}} e^{-(v-v_0)^2/2\sigma_v^2},$$ thus $$\int I_{\nu} dv = I_{\nu,\text{peak}}\sqrt{2\pi}\sigma_{v} = I_{\nu,\text{peak}}\times\text{FWHM}\times\sqrt{\pi/(4\ln2)}.$$
- `flux_int`: $$(I_\nu - I_{\nu,\text{continuum}}) \times4\pi \times\delta\nu_\text{FWHM}  \times\sqrt{\pi/(4\ln2)},$$ thus it is the spectrally integrated intensity times $4\pi$.
- `flux_Jy`: $$(I_\nu - I_{\nu,\text{continuum}}) \times \Omega_\text{beam} \times10^{23},$$ where $\Omega_\text{beam}$ is the beam solid angle $$\Omega_\text{beam} = \frac{\pi}{4\ln2}(\theta_\text{FWHM})^2.$$
- `beta`: escape probability ($\beta$)
- `Jnu`: $$\frac{A_{ul}y_{u}}{B_{lu}y_{l} - B_{ul}y_{u}} (1-\beta) + J_\text{cont,bg}\beta + J_\text{cont,bulk}$$
- `gup`: statistical weight of upper level
- `glow`: statistical weight of lower level
- `Aul`: Einstein A coefficient
- `Bul`: Einstein B coefficient
- `Blu`: Einstein B coefficient
- `n_crit`: critical density according to Shirley2015, equation (9)
- `n_crit_old`: critical density based on the two-level simplification
- `q`: quality of the computation

---

## Calculation of critical density

By default equation (9) of [Shirley2015](https://iopscience.iop.org/article/10.1086/680342) is used.

A simpler definition that is only valid for two-level system is also used (labeled as "n_crit_old" in the output).  This is described in the text below equation (5) of [Shirley2015](https://iopscience.iop.org/article/10.1086/680342).

---

This code solves essentially the same problem as
[`RADEX`](http://home.strw.leidenuniv.nl/~moldata/radex.html) written by Van der Tak, except that we take a different approach to solve the statistical equilibrium problem.
Given an initial distribution, what we do is to evolve the system towards equilibrium using an ODE solver.

## Usage

To use this code, you first need to compile it using the `makefile` (the executable is named `my_radex`) by running `make` in the command line, then edit the configuration file (`configure.dat`) to meet your needs.  After this execute the command
```
./my_radex configure.dat
```
and you will get the results.

A python wrapper is also included.  To make the wrapper, run `make wrapper` in the command line.
For its usage, see [this Jupyter notebook](https://github.com/fjdu/myRadex/blob/master/example.ipynb).
The wrapper is preliminary; not all the functionalities in the Fortran source code are included in the wrapper (though usually they are not needed).

We use the [LAMDA](http://home.strw.leidenuniv.nl/~moldata/molecules.html) format for the input energy levels and transition rates.

This code has not been thoroughly tested, and it is possible that some input files won't load properly.
Usually one can work around the problem by changing the input file to the format of a file that is known to work.
