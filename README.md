# BaryonEffectsEmulator
It is known that N-body simulation cannot trace the effects of baryons at medium and small cosmological scales. 
With this emulator we propose a fast and accurate remedy to this problem. 
For an accurate description of the BaryonEffectsEmulator refer to the research project MSc Thesis in the documentation folder.

## Content
- **Main source** 
  * *BaryonEffectsEmulator.py*: main source code of the emulator.
  * *bee.pkl*: data file on whih the emulator is constructed.
- **Documentation**
  * *MSc Thesis*: research project with detailed documentation of the emulator.
- **Example**
  * *get_boost_example.py*: example code for the generation of the boost.
  * *plot_example.py*: example code for plotting the boost.
  
  
## Installation
- **Prerequisites**: *pandas*, *scipy*, *numpy*.
- The installation process is straightforward, one has simply to download the ZIP file. Consequntely one has to import the functions located in the file *BaryonEffecsEmulator.py* into its own python code.

## Usage

This section explains precedure to follow in order to get the boost P<sub>BCM</sub>/P<sub>DMO</sub>. After the download of the files one can create a new python code and proceed as follow:
```ruby
from BaryonEffectsEmulator import *
```
If the python file *BaryonEffectsEmulator.py* is located in a different folder, then first change the path for example with
```ruby
import sys
sys.path.insert(0, '/path/to/BaryonEffectsEmulator')
```
Once this is done one has to define a baryonic correction model and the redshifts at which the boost has to be calculated, e.g.
```ruby
MyBCM = {'f_b': 0.14, 'logMc': 14, 'mu': 0.4, 'theta_ej':4.0, 'eta_tot': 0.3, 'eta_cga': 0.6}
z=[0.3,0.9,1.5]
```
As last step one should simply call the function *get_boost* and access the results as follows:
```ruby
result=get_boost(z,MyBCM)
k=result['k']
boost0=result['z0']
boost1=result['z1']
boost2=result['z2']
```
