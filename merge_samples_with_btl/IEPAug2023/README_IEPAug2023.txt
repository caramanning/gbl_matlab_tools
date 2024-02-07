Goal of the code

We have sample data from the GC and we need to get the T, S, and density data for these samples from the CTD file.


import_IEP_btl.m:
First use import data function in matlab to load in the bottle file. Can use Generate Script to make code that you can run, or import data to just make a variable.
Save as IEP_btl.mat

-------
import_IEP_SID.m:
Import the sample IDs that correspond to each station and niskin number.


-----
get_TS_data: this is used to match up the data for an individual run with the corresponding T, S, Density
use the daily sample run file and the text to columns function to make a table with the station, niskin and sample ID

copy the sample IDs in at the start of the matlab file and then run the code
export the CSV that is generated, open in excel, and then copy into the daily run file

then there is code at the end of the file to calculate the density for the air-eq files and these can also be copied into excel

