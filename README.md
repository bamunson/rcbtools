# rcbtools
An analysis toolkit for RCB stars simulated in MESA

## Installing/Uninstalling rcbtools

To install rcbtools, simply clone or download this repo, then `cd` into it and type the following command:

`pip install .`

To uninstall rcbtools, type the following command:

`pip uninstall rcbtools`

## Using rcbtools

Once installed, you should be able to access the tools like you would import any package:

`import rcbtools`

After importing rcbtools, you can begin to use some of the tools as you wish. Here are some examples to get you started:

```
import rcbtools
import matplotlib.pyplot as plt

p = rcbtools.profile2dict('profile1.data') # Create a dict of values from a MESA profile (column headers are keys)
plt.semilogy(p["mass"],10**p["log_T"]) # Create log plot of temperature v mass

abunds = rcbtools.makeabund('profile1.data') # Create a dict of elements from MESA profile
plt.loglog(p["mass"],abunds["C"]) # Create a log plot of total carbon abundance (sum of all isotopes) v mass

rcbtools.surfabund('profile1.data') # Make surface abundance plot and print surface information

rcbtools.surfabund2('profile1.data',elements=['Li','C','N','O','Ne']) # Same as above, but specify elements

h = rcbtools.profile2dict('history.data') #Import history datafile as python dict

plt.plot(h['log_Teff'],h['log_L']) # Create HR plot from history datafile
plt.xlim(max(h['log_Teff']),min(h['log_Teff'])) # reverse x-axis
```
