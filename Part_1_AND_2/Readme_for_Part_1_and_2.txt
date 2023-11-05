With this MATLAB program, "Satellite_Project_Final_Version_Part_1_and_2" you can place one receiver on the moon and one satellite around the moon.

To use this program, follow these steps:

Step_1:there is a CSV file where you can enter the specific information about your satellite and receiver. you can also enter some general data.
When you run the program, the program reads this file (Lines 26_27 in the main program).

The following parameters can be changed:
Initial geodetic coordinates of the receiver in degrees.
Initial true anomaly of the satellite in radians.
Keplerian parameters 'a' (in meters) and 'e'.
Satellite orbit inclination based on Keplerian parameters (W,I,Omega) in degrees.
Time step, t_u and delta_t (all in seconds), and the Masking_Angle in degrees.

Step_2:Finally, you can run the program,it will process the input parameters and simulate the satellite's movement around the moon.
The running process will finish after one orbital period, but the program is flexible and you can change ending_time (Line 106).
After finishing the runnig process you will see 2 plots (Pseudorange,Pseudorange_Rate).
