# Matlab_code
The experimental data are given in the files "Ch1603_experimental_data.xlsx"  and "Ch1616_experimental_data.xlsx"   for validation purposes.

 1. For Ch1603 
    1. Ch1603_28_182.m:  parameter estimation for  28 to 182 days. First, run this file for Ch1603 because Ch1603_182_385.m needs some data from the output of Ch1603_28_182.m.
    2. Ch1603_182_385.m:  parameter estimation for  182 to 385 days.
    3. Ch1603_validation_with_experimental_data.m:  To see the final plots, run this file.
    4. Ch1603_experimental_data.xlsx:  Experimental data. The corresponding article is mentioned in the manuscript.

     
 2. For Ch1616
    1. Ch1616_28_182.m:  parameter estimation for  28 to 182 days. First, run this file for Ch1616 because Ch1616_182_413.m needs some data from the output of Ch1616_28_182.m
    2. Ch1616_182_413.m:  parameter estimation for  182 to 413 days.
    3. Ch1616_validation_with_experimental_data.m: To see the final plots, run this file.
    4. Ch1616_experimental_data.xlsx:  Experimental data. The corresponding article is mentioned in the manuscript.

 3. fminsearchbnd.m: Keep this m-file in the same folder. This function minimizes errors in the time of parameter estimation.
