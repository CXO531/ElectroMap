# ElectroMap

ElectroMap - optical and electrophysiological mapping software.

5th Jan 2020 - 
ONLY WORK WITH MATLAB 2019B CURRENTLY

Updated version of software including optical wave similarity and dFdt mapping can be found in new folder 'Source Code - OWS and dFdt update'. .exe version of updated software will be once it has been fully tested. If any issues are found in new version, please contact C.O'Shea@bham.ac.uk.

See below for more details on methodology. 
Setting are available from top menu. Setting dFdt level to 0 will map maximum dFdt. 
Setting OWS threshold to zero will calculate OWS, setting to any other value will calculate regularity index. 

dFdt - https://www.frontiersin.org/articles/10.3389/fphys.2019.01295/full

OWS/RI - https://www.sciencedirect.com/science/article/pii/S0079610719302135



The software is available as 

[1] MATLAB source code. To run this version, MATLAB is required along with the Image Processing, Signal Processing, Statistics and Machine Learning and the Curve Fitting toolbox. Recommended for use with MATLAB 2016b or later. Previous MATLAB releases have not been fully tested, altough can be used to run the software. 

To run the software -

Download/clone ALL files in 'Source Code' folder.
Unzip files
Open and run 'ElectroMap.m'

[2] Standalone .exe (windows) or .dmg (MAC) applications

Installer for standalone versions of the software can be found at https://drive.google.com/open?id=1nJyI07w9WIt5zWcit0aEyIbtg31tANxI. Installation will require internet connection to download MATLAB runtime, and are called -

ElectroMapInstaller_web.exe (Windows)

(MAC)

Alternatively, the MATLAB runtime can be separately downloaded from https://uk.mathworks.com/products/compiler/matlab-runtime.html. Version 9.4 is required for use of ElectroMap. Further detailed instructions are included in the readme.txt file on the google drive.

Two example datasets are also included in the google drive for use in the software - 

'Atria.tif' - Isolated mouse left atria, paced at 3 and 8.33Hz and loaded with voltage sensitive dye. Sampling rate =0.987kHz, pixel size = 71.4micrometers

'GP.tif' - Guinea Pig whole heart paced at 5Hz. Sampling rate =0.5kHz, pixel size = 320micrometers

For any questions, please contact C.O'Shea@bham.ac.uk.
