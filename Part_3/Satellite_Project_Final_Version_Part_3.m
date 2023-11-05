clear all; %#ok
close all;
clc;

%General data
time_step = 300;
ending_time = 172650;
k=time_step:time_step:ending_time;
Masking_Angle =0;
satelite_antenna_cone_angle =90;
t_transmit=0;
delta_t=0;

%Reading The Inputs
satellites = csvread('satellites_for_the_whole_moon_surface_rectangular_configuration.csv',1,0);

%Computing the position of the satellites
number_of_satellites = size(satellites);
number_of_satellites=number_of_satellites(1);

%Pre_Allocation
all_sat_pos = nan(4,length(k)+1,number_of_satellites);

for i=1:number_of_satellites
    
    all_sat_pos(:,:,i) = Compute_Satellite_Position(satellites(i,1),satellites(i,2),satellites(i,3),satellites(i,4),satellites(i,5),satellites(i,6),time_step,ending_time);
    
end

disp('Satellites Pose computaion is done.');
disp(['Number of the Satellites: ',num2str(number_of_satellites)]);
disp(' ');

%Define the initial position of the receivers
rec_lat_start_point = -90;
rec_lat_end_point =90;
rec_lat_res = 5;

rec_long_start_point = -180;
rec_long_end_point = 180;
rec_long_res =10;

rec_lat_array = rec_lat_start_point:rec_lat_res:rec_lat_end_point;
rec_long_array = rec_long_start_point:rec_long_res:rec_long_end_point;

%Pre_Allocation
receivers = nan(length(rec_lat_array)*length(rec_long_array),2);

%Calculate Cartesian product
for i=1:length(rec_lat_array)
    
    for j=1:length(rec_long_array)
        
        receivers((i-1)*length(rec_long_array)+j,1) = rec_lat_array(i);
        receivers((i-1)*length(rec_long_array)+j,2) = rec_long_array(j);
        
    end
    
end

%Computing the position of the receivers
number_of_receivers = size(receivers);
number_of_receivers=number_of_receivers(1);

%Pre_Allocation
all_rec_pos = nan(4,length(k)+1,number_of_receivers);

for i=1:number_of_receivers
    
    all_rec_pos(:,:,i) = Compute_Receiver_Position(receivers(i,1),receivers(i,2),time_step,ending_time);
    
end

disp('Receivers Pose computaion is done.');
disp(['Number of the Receivers: ',num2str(number_of_receivers)]);

%Calculate in_view_total
in_view_total = zeros(number_of_receivers+1,length(k)+1);
in_view_total(1,:)=0:time_step:ending_time;

for i=1:number_of_receivers
    
    for j=1:number_of_satellites
        
        in_view = evaluate_sat_wrt_rec(all_sat_pos(:,:,j),all_rec_pos(:,:,i),Masking_Angle,satelite_antenna_cone_angle);
        in_view_total(i+1,:) = in_view_total(i+1,:) + in_view(2,:);
        plot(in_view(1,:),in_view(2,:))
        
    end
    
    disp(' ');
    disp(['Evaluating Rec Number ',num2str(i),' is done.']);
    
end

%Drawing Plot 1
for i=2:(number_of_receivers+1) 
    
    figure(1);
    plot(in_view_total(1,:),in_view_total(i,:));
    title('In View Total.');
    xlabel('Time in Seconds.');
    ylabel('Number of satellites in view.');
    hold on;
    
end

%Calculate The Actual Coverage for each receiver
coverage_percent_rec = zeros(1,number_of_receivers);

for i=2:(number_of_receivers+1) 
    
    total_duration_in_coverage = 0;
    total_duration = (length(k))*(time_step);
    
    for j=1:(length(k)+1)
        
        if in_view_total(i,j) > 3
            
            total_duration_in_coverage = total_duration_in_coverage+time_step; 
            
        end
        
    end
    
    coverage_percent_rec(i-1) = (total_duration_in_coverage/total_duration)*100;
    
end

%Calculate The CDF of the Actual Coverage for all receivers
coverage_cdf = zeros(1,100);

for i=1:100
    
    number_of_rec_in_coverage=0;
    
    for j=1:number_of_receivers
        
        if coverage_percent_rec(j)>i
            
            number_of_rec_in_coverage = number_of_rec_in_coverage+1;
            
        end
        
    end
    
    coverage_cdf(i) =number_of_rec_in_coverage;
    
end

%Drawing Plot 2
figure(2);
plot(coverage_cdf);
title('Actual Coverage.');
xlabel('Fraction of ending time in percent.');
ylabel('The total number of receivers that see at least 4 satellites.');