function sat_pos = Compute_Satellite_Position(a,e,I,W,Omega,true_anomaly_0,time_step,ending_time)

%Define the start time
Time(1)=0;   
k=time_step:time_step:ending_time;

G = 6.67428*(10^(-11)); % Gravitational constant in [m^3 kg^-1 s^-2]
M = 7.34767309*(10^22); % Mass of the Moon in [kg]

%Define the rotation matrix around X-axis,Y-axis and Z-axis
r_m=rotation_matrix(I,W,Omega);

%Pre_Allocation
sat_pos = nan(4,length(k)+1);
sat_ang = nan(1,length(k)+1);

%Calculate the satellite's initial position in Cartesian Coordinate System
sat_ang(1) = true_anomaly_0;
r_t_minus_1 = ((a*(1-e^2))/(1+e*cos(true_anomaly_0)));
x_s = r_t_minus_1*cos(sat_ang(1));
y_s = r_t_minus_1*sin(sat_ang(1));
sat_pos_0=[x_s;y_s;0];
sat_rotated_position = r_m*sat_pos_0;         
sat_pos(:,1)=[Time(1);sat_rotated_position(1);sat_rotated_position(2);sat_rotated_position(3)];

i=1;

for (k=time_step:time_step:ending_time) %#ok
    
    Time(i+1)=Time(i)+time_step; %#ok
    
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
 
     end
    
    omega(:,i) =[Time(i); (v/r_t_minus_1)*sin(Alpha)];%#ok
    sat_ang(i+1)=(omega(2*i)*time_step)+sat_ang(i); 
    r_t =((a*(1-e^2))/(1+e*cos(sat_ang(i+1))));  
    r_t_minus_1 = r_t;
       
    %Calculate the satellite's new position in Cartesian Coordinate System
    x_s = r_t*cos(sat_ang(i+1));
    y_s = r_t*sin(sat_ang(i+1));
    sat_rotated_position = r_m*[x_s;y_s;0];         
    sat_pos(:,i+1)=[Time(i+1);sat_rotated_position(1);sat_rotated_position(2);sat_rotated_position(3)];
           
    i=i+1;
end

end