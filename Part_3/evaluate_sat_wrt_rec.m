function in_view = evaluate_sat_wrt_rec(sat_pos,rec_pos,masking_angle,satelite_antenna_cone_angle)

       %Pre_Allocation
       in_view = nan(2,length(sat_pos));
       
       %Setting timestamps
       in_view(1,:) = sat_pos(1,:);
              
       %Define the radius of the moon in [m]
        moon_radius = 1737500;
        
       %Compute wether the satellite is insight of the receiver or not
       
       for (i=1:1:length(sat_pos)) %#ok
             
       insight_threshold_angle= masking_angle+90;
              
       theta=pi-asin((norm(sat_pos(2:4,i))/moon_radius)*sind(satelite_antenna_cone_angle));
       theta=theta*(180/pi);
       
      sat_inview_angle = max(insight_threshold_angle,theta); 
       
                    vec1 = -rec_pos(2:4,i);
                    vec2 =  sat_pos(2:4,i) - rec_pos(2:4,i) ;
                    beta= acos((vec1'*vec2)/(norm(vec1)*norm(vec2)));
                    beta= beta*(180/pi);
                     
                         if (sat_inview_angle<beta && beta <= 180)
                             
                             in_view(2,i) = 1;
                             
                         else
                             
                             in_view(2,i) = 0;
                             
                         end
                         
      end

end