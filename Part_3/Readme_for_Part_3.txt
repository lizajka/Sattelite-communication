With this MATLAB program, "Satellite_Project_Final_Version_Part_3" you are able to place an arbitrary number of satellites around the moon and an arbitrary number of receivers on the moon.

To use this program, follow these steps:

Step_1: there is a CSV file where you can enter the specific information about your satellites. When you run the program, the program reads this file (Lines 14-15 in the main program).
the first caloumn is Semi-major axis in meters, the second column is Ecentricity, then we have three inclinations based on Keplerian parameters (W,I,Omega) in degrees and finally,
the last column is the initial true anomaly of the satellites in radians.

Step_2: You have to set the general data (Lines 5-12 in the main program).

Step_3: You have to choose a specific part of the moon with your desired resolution (Lines 34-59 in the main program) and place your receivers there.

Step_4: Finally, you can run the program and after finishing the runnig process you will see 2 plots (In_View_Total and Actual_Coverage).



