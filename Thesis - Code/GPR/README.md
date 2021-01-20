# Average Square Difference Function Layer Picking System (ALPS)
- A GPR layer tracing algorithm
- Described in Chapter 2 of my [PhD thesis](../Thesis_AustinLines_FINAL.pdf)

## CONTENTS
[ALPS.mlappinstall](./ALPS.mlappinstall)
> Install file for a MATLAB application. Just open it with MATLAB and it will install, or go to APPS -> Install App and open the file. Once installed, it will be listed under the APPS tab. If you prefer to have the .m files, then extract all the files (use any unzip tool on the .mlappinstall file). The top level file is called "GPR_LayerTracker.mlapp"

[GPR_Example_1.DZT](./GPR_Example_1.DZT) / [GPR_Example_2.DZT](./GPR_Example_2.DZT)
> Example GPR files collectd with a GSSI SIR30 in Greenland that have been preprocessed to test the functionality of the ALPS application.

[GPR_Example_1_layers.mat](./GPR_Example_1_layers.mat) / [GPR_Example_2_layers.mat](./GPR_Example_2_layers.mat)
> Layer files to test the "Import Layers" function of the ALPS application.

[GPR_Example_1.DZG](./GPR_Example_1.DZG) / [GPR_Example_2.DZG](./GPR_Example_2.DZG)
> Native GPS record for the SIR30 of where the data was collected. The ALPS application automatically searches for either a .DZG or .kml file to geolocate the layers when opening a GPR file (.dzt).

[GPR_Example_1.KML](./GPR_Example_1.KML) / [GPR_Example_2.KML](./GPR_Example_2.KML)
> KML file of GPS record of where the data was collected. The ALPS application automatically searches for either a .DZG or .kml file to geolocate the layers when opening a GPR file (.dzt).
