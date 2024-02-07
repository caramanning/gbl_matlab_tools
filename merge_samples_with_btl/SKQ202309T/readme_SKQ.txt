Cara Manning
2024/02/06


goal of the code
-import SKQ202309 CTD Niskin bottle data (.btl) into matlab
-import gas sample IDs into matlab
-find the T/S data for each GC run of samples
-export the sample data back to matlab and add in other CTD data


-CTD files are taken from: C:\Users\ccm21008\OneDrive - University of Connecticut\research\AICC 2023\SKQ202309T\data\SKQ202309T\ctd\proc

-To delete header manually (there are ways to do this in git, but this is is a manual way)
-place cursor at the first point in file you want to keep, CTRL+SHIFT+HOME, then DELETE or BACKSPACE to delete
-remove rows 1-273 (this will vary for each file)

Then import new file into matlab
- Column delimiters: Space, Delimiter Options: Treat multiple delimiters as one
- Range A3:W42
Variable names row: 1
Then ned to rename the variable names because the date appears across multiple columns
-Import selection: generate script (useful if you are going to need to importa a series of similar files)

-Then add code to remove the rows that contain the standard deviation and save the file. Repeat for all CTD files.

