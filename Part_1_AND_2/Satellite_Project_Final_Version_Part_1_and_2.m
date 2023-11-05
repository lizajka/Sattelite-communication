%We assumed that the center of the moon is at the origin of the Cartesian Coordinate System and the receiver is fixed on the moon surface.
%START_TIME is equal to ZERO.
%The radius of the moon is 1737500 [m].
%The moon rotates around the Z-axis at a constant speed.
%The initial orbit is on the X-Y plane.
%The running process will finish after one orbital period, but the program
%is flexible and you can change ending_time (Line 106).
%With these assumptions, you need to give the program the initial geodetic coordinates of the receiver (in degrees) and the initial true anomaly of the satellite (in radians).
%You should also give the program Keplerian parameters a (in meters) and e as well as satellite orbit inclination based on Keplerian parameters (in degrees).
%Finally, you need to give time_step,t_u,delta_t (in seconds) and Masking_Angle (in degrees) to the program.
%Be aware that when the satellite is in sight of the receiver, the receiver_color is green and when it is not in sight of the receiver the receiver_color turns to red.

clear all; %#ok
close all;
clc;

%Define the start time
Time(1)=0;

%Define the radius of the moon in [m]
moon_radius = 1737500;

%Define the speed of light in [m/s]
c=299792458;

%Reading The Inputs
input_parameters = csvread('Input_Parameters.csv',1,0);

Lat_0_R=input_parameters(1,1);
Long_0_R=input_parameters(1,2);
true_anomaly_0=input_parameters(1,3);
a=input_parameters(1,4);
e=input_parameters(1,5);
W=input_parameters(1,6);
I=input_parameters(1,7);
Omega=input_parameters(1,8);
time_step=input_parameters(1,9);
t_u=input_parameters(1,10);
delta_t=input_parameters(1,11);
Masking_Angle=input_parameters(1,12);

%Calculate the receiver's initial position in Cartesian Coordinate System.
X_0_R = moon_radius*cosd(Lat_0_R)*cosd(Long_0_R);
Y_0_R = moon_radius*cosd(Lat_0_R)*sind(Long_0_R);
Z_0_R = moon_radius*sind(Lat_0_R);

%Calculate the satellite's initial position in Cartesian Coordinate System.
sat_ang(1) = true_anomaly_0;
r_t_minus_1 = ((a*(1-e^2))/(1+e*cos(sat_ang(1))));
X_0_S = r_t_minus_1*cos(sat_ang(1));
Y_0_S = r_t_minus_1*sin(sat_ang(1));
Z_0_S = 0;

sat_pos_0=[X_0_S;Y_0_S;Z_0_S];

%Generate a grid of points to plot the moon
[X,Y,Z] = sphere(100);

%Scale the grid by the moon's radius
X = X * moon_radius;
Y = Y * moon_radius;
Z = Z * moon_radius;

%Plot the moon
surf(X,Y,Z, 'FaceColor', [1 1 1], 'FaceAlpha', 1)

%Set the axis limits
axis equal;
axis([-10000000 10000000 -10000000 10000000 -10000000 10000000])

%Set the plot title and axis labels
title('Satellite Orbiting the Moon.')
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

%Hold on to add additional plots
hold on;

%Calculate the moon's angular speed in [r/s] and display it
rotational_period=27.322*24*60*60;
moon_angular_speed=(2*pi)/(rotational_period);
disp(['Moon Angular Speed: ' num2str(moon_angular_speed) ' [r/s]']);

%Define the rotation matrixes based on Keplerian Parameters
rot_W = [cosd(W),sind(W),0;-sind(W),cosd(W),0;0,0,1];

rot_I = [1,0,0;0,cosd(I),sind(I);0,-sind(I),cosd(I)];

rot_Omega = [cosd(Omega),sind(Omega),0;-sind(Omega),cosd(Omega),0;0,0,1];

r_m=rot_Omega*rot_I*rot_W;

%Plot the rotated orbit
r_0_minus_eps=((a*(1-e^2))./(1+e*cosd(0:360)));
Oribit_points = [r_0_minus_eps.*cosd(0:360);r_0_minus_eps.*sind(0:360);zeros(1, 361)];
rotated_orbit_points = r_m*Oribit_points;
plot3(rotated_orbit_points(1,:),rotated_orbit_points(2,:), rotated_orbit_points(3,:),'g')

%Calculate the Period of the Orbit in [s] and display it
G = 6.67428*(10^(-11)); % Gravitational constant in [m^3 kg^-1 s^-2]
M = 7.34767309*(10^22); % Mass of the Moon in [kg]
period_of_the_orbit=(2*pi)*(sqrt((a^3)/(G*M))); 
disp(['Period of The Orbit: ' num2str(period_of_the_orbit) ' [s]']);

%Setting ending_time
ending_time=period_of_the_orbit;

%Define the satellite's intitial position,size and color
satellite_initial_position=r_m*sat_pos_0;
sat_pos = [Time(1);satellite_initial_position(1);satellite_initial_position(2);satellite_initial_position(3)];
sat_size = 50;
sat_color = 'm';

%Plot the satellite's initial position
Satellite_Position = scatter3(sat_pos(2), sat_pos(3), sat_pos(4), sat_size, 'filled', 'MarkerFaceColor', sat_color);

%Define the receiver's initial position and size at t=0;
rec_pos = [Time(1);X_0_R;Y_0_R;Z_0_R];
rec_size=15;

%Calculate the receiver's initial angle (Theta) in Spherical Coordinate System
Theta=acos(Z_0_R/moon_radius);

%Calculate the receiver's initial angle (Phi) in Spherical Coordinate System 
if(Theta==0||Theta==pi)
    
Phi(1)=0;

else
    
Phi(1)=acos(X_0_R/(moon_radius*sin(Theta)));

end

 %Determine wether the satellite is insight of the receiver or not
 insight_threshold_angle = (Masking_Angle+90)*(pi/180);
 vec1 = -rec_pos(2:4,1);
 vec2 =  sat_pos(2:4,1) - rec_pos(2:4,1) ;
 angle_between_sat_and_rec = acos((vec1'*vec2)/(norm(vec1)*norm(vec2)));
    
 if (angle_between_sat_and_rec>insight_threshold_angle)
        
     rec_color = 'g';
        
 else
        
     rec_color = 'r';
        
 end
    
%Plot the receiver's initial position
Receiver_Position = scatter3(rec_pos(2), rec_pos(3), rec_pos(4), rec_size, 'filled', 'MarkerFaceColor', rec_color);

%Calculate the initial distance between the satellite and the receiver and also pseudorange
Distance(:,1)=[Time(1);sqrt((sat_pos(2)-rec_pos(2))^2+(sat_pos(3)-rec_pos(3))^2+(sat_pos(4)-rec_pos(4))^2)];
Pseudorange(:,1)=[Time(1);Distance(2)-c*(t_u-delta_t)];

%Calculate step duration for the receiver
d=moon_angular_speed*time_step;

i=1;
j=1;

%Loop over time steps

for (k=time_step:time_step:ending_time) %#ok
    
    Time(i+1)=Time(i)+time_step; %#ok
    
    %Calculate the receiver's new angle
    Phi(i+1) =Phi(i)+d; %#ok
    
    %Calculate the receiver's new position
    rec_pos(:,i+1)=[Time(i+1);moon_radius*sin(Theta)*cos(Phi(i+1));moon_radius*sin(Theta)*sin(Phi(i+1));moon_radius*cos(Theta)];
    
    %Calculate the satellite's new position in Polar system
    v = sqrt(G*M*(2/r_t_minus_1 - 1/a));
    
    if (i<2)
        
     r_0_minus_eps = ((a*(1-e^2))/(1+e*cos(true_anomaly_0-10*eps)));
     X_0_minus_eps = r_0_minus_eps*cos(true_anomaly_0-10*eps);
     Y_0_minus_eps = r_0_minus_eps*sin(true_anomaly_0-10*eps);
     Z_0_minus_eps = 0;
     
     sat_pos_0_minus_eps=[X_0_minus_eps;Y_0_minus_eps;Z_0_minus_eps];

     vec1 = -sat_pos_0;
     vec2 = sat_pos_0 - sat_pos_0_minus_eps;
     Alpha = acos((vec1'*vec2)/(norm(vec1)*norm(vec2)));
    else
        
    vec1 = -sat_pos(2:4,i);
    vec2 = sat_pos(2:4,i) - sat_pos(2:4,i-1) ;
    Alpha = acos((vec1'*vec2)/(norm(vec1)*norm(vec2)));
    Zaviye=Alpha*(180/pi);
    end
    
    omega(:,i) =[Time(i); (v/r_t_minus_1)*sin(Alpha)];%#ok
    sat_ang(i+1)=(omega(2*i)*time_step)+sat_ang(i); %#ok
    r_t =((a*(1-e^2))/(1+e*cos(sat_ang(i+1))));  
    r_t_minus_1 = r_t;
       
    %Calculate the satellite's new position in Cartesian Coordinate System
    x_s = r_t*cos(sat_ang(i+1));
    y_s = r_t*sin(sat_ang(i+1));
    sat_rotated_position = r_m*[x_s;y_s;0];         
    sat_pos(:,i+1)=[Time(i+1);sat_rotated_position(1);sat_rotated_position(2);sat_rotated_position(3)];
    
    %Determine wether the satellite is in sight of the receiver or not
    insight_threshold_angle = (Masking_Angle+90)*(pi/180);
    vec1 = -rec_pos(2:4,i+1);
    vec2 =  sat_pos(2:4,i+1) - rec_pos(2:4,i+1) ;
    angle_between_sat_and_rec = acos((vec1'*vec2)/(norm(vec1)*norm(vec2)));
    
    if (angle_between_sat_and_rec>insight_threshold_angle)
        
        rec_color = 'g';
        
    else
        
        rec_color = 'r';
        
    end

    %Update the receiver's new position
    set(Receiver_Position,'XData', rec_pos(2,i+1),'YData', rec_pos(3,i+1),'ZData', rec_pos(4,i+1), "MarkerEdgeColor", rec_color, "MarkerFaceColor", rec_color);
    
    %Update the satellite's new position
    set(Satellite_Position, 'XData', sat_pos(2,i+1), 'YData', sat_pos(3,i+1), 'ZData', sat_pos(4,i+1));

    %Refresh the plot
     drawnow;
    
    %Calculate the new distance between the satellite and the receiver and
    %also Pseudorange
     Distance(:,i+1)=[Time(i+1);sqrt((sat_pos(j+5)-rec_pos(j+5))^2+(sat_pos(j+6)-rec_pos(j+6))^2+(sat_pos(j+7)-rec_pos(j+7))^2)]; 
     Pseudorange(:,i+1)=[Time(i+1);Distance(2*(i+1))-c*(t_u-delta_t)];
     
    i=i+1;
    j=j+4;
end

%Calculate Pseudorange_Rate
Pseudorange_Rate=diff(Pseudorange(2,:))/(time_step);

%Plot the results
figure(2);
plot(Pseudorange(1,:),Pseudorange(2,:));
title('Pseudorange over time.');
xlabel('Time (s)');
ylabel('Pseudorange (m)');

figure(3);
plot(Pseudorange_Rate);
title('Pseudorange Rate over time.');
xlabel('Number of Sample.');
ylabel('Pseudorange Rate (m/s)');