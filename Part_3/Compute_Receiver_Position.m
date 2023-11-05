function rec_pos = Compute_Receiver_Position(Lat_0_R,Long_0_R,time_step,ending_time)

%Define the start time
Time(1)=0;
k=time_step:time_step:ending_time;

%Define the radius of the moon in [m]
moon_radius = 1737500;

%Calculate the receiver's initial position in Cartesian Coordinate System.
X_0_R = moon_radius*cosd(Lat_0_R)*cosd(Long_0_R);
Y_0_R = moon_radius*cosd(Lat_0_R)*sind(Long_0_R);
Z_0_R = moon_radius*sind(Lat_0_R);

%Pre_Allocation
rec_pos = nan(4,length(k)+1);
Phi = nan(1,length(k)+1);

rec_pos(:,1) = [Time(1);X_0_R;Y_0_R;Z_0_R];

%Calculate the moon's angular speed in [r/s]
rotational_period=27.322*24*60*60;
moon_angular_speed=(2*pi)/(rotational_period);

%Calculate the receiver's initial angle (Theta) in Spherical Coordinate System
Theta=acos(Z_0_R/moon_radius);

%Calculate the receiver's initial angle (Phi) in Spherical Coordinate System 
if(Theta==0||Theta==pi)
    
Phi(1)=0;

else
    
Phi(1)=acos(X_0_R/(moon_radius*sin(Theta)));

end

%Calculate step duration for the receiver
d=moon_angular_speed*time_step;

i=1;

for (k=time_step:time_step:ending_time) %#ok
    
 Time(i+1)=Time(i)+time_step; %#ok
    
 %Calculate the receiver's new angle
 Phi(i+1) =Phi(i)+d;
    
 %Calculate the receiver's new position
 rec_pos(:,i+1)=[Time(i+1);moon_radius*sin(Theta)*cos(Phi(i+1));moon_radius*sin(Theta)*sin(Phi(i+1));moon_radius*cos(Theta)];
 
 i=i+1;
end

end