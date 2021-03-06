# rcbtools
An analysis toolkit for RCB stars simulated in MESA

## Installing/Uninstalling rcbtools
This toolkit is now in the PyPI repository. This means it can be downloaded via `pip install rcbtools` in the command line.

Alternatively, simply clone or download this repo, then `cd` into it and type the following command:

`pip install .`

If you do not have PyPI installed, you can download this repo and use the following command:

`python setup.py install`

To uninstall rcbtools, type the following command:

`pip uninstall rcbtools`

Otherwise, remove all `rcbtools` directories inside your `$PYTHONPATH`.

## Using rcbtools

Once installed, you should be able to access the tools like you would import any package:

`import rcbtools`

After importing rcbtools, you can begin to use some of the tools as you wish. Here are some examples to get you started:

```
import rcbtools
import matplotlib.pyplot as plt

p = rcbtools.profile2dict('profile1.data') # Create a dict of values from a MESA profile (column headers are keys)
p.keys() # See a list of available keys for the dictionary "p"
plt.semilogy(p["mass"],10**p["logT"]) # Create log plot of temperature v mass

abunds = rcbtools.makeabund('profile1.data') # Create a dict of elements from MESA profile
plt.loglog(p["mass"],abunds["C"]) # Create a log plot of total carbon abundance (sum of all isotopes) v mass

rcbtools.rcbsurf('profile1.data') # Make surface abundance plot and print surface information

h = rcbtools.profile2dict('history.data') #Import history datafile as python dict

plt.plot(h['log_Teff'],h['log_L']) # Create HR plot from history datafile
plt.xlim(max(h['log_Teff']),min(h['log_Teff'])) # reverse x-axis
```
