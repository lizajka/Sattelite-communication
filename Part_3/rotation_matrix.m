function  r_m=rotation_matrix(I,W,Omega)

rot_W=[cosd(W),sind(W),0;-sind(W),cosd(W),0;0,0,1];
     
rot_I=[1,0,0;0,cosd(I),sind(I);0,-sind(I),cosd(I)];

rot_Omega=[cosd(Omega),sind(Omega),0;-sind(Omega),cosd(Omega),0;0,0,1];

r_m=rot_Omega*rot_I*rot_W;

end