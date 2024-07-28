%Attempt to upscale(FLOW BASED UPSCALING METHOD)
%PETROLEUM DESIGN PROJECT

%STUDENT NAME - SAN HTOO AUNG
%ID - 10368431
%EMAIL - san.aung@student.manchester.ac.uk

clear;
clc;
Kx = load('PERMX_5 - Copy.GRDECL');        %loading x-perm from GRDECL file
Kx = reshape(Kx,120,60,10);     %Turning the vector into array
Kx = flip(Kx,2);       %Flipping
Kx_z = Kx(:,:,1);  %layer number

p_in = 100;      %Assigning inlet Pressure to 100 bar
p_out = 0;      %Assigning outlet Pressure to 90 bar
mew = 0.05;       %Assigning arbitrary viscosity 

Nx = 60;       %Number of grid blocks in x-direction 
Ny = 120;        %Number of grid blocks in y-direction
Nz = 10;        %Number of layer in z-direction 
Lx = 60;       %length of the reservoir model in x-direction  
Ly = 60;       %length of the reservoir model in y-direction
Lz = 60;       %length of the reservoir model in z-direction 
layer = 1:Nz;   %Layering of z

%GRIDDING
dx = Lx/Nx; % grid size in x-direction for x direction flow
dy = Ly/Ny; % grid size in y-direction for y direction flow
Area = (2*dx)*(2*dy);    %Area of a reservoir
Kf = zeros(60,30);      %targeted upscaling parameters

%% UPSCALING FOR FLOW BASED MODEL [X]

%Upscaling in a way that 2x2 blocks into 1 single block
for I = 1:Ny/2
   for J = 1:Nx/2
       N = 2*2;              % number of unknowns
       A = zeros(N,N);       %N rows and N columns for coefficient matrix
       B = zeros(N,1);       %N rows and 1 column for the right hand side vector
       Perm = Kx(2*I-1:2*I,2*J-1:2*J,10); %rows and columns are 2 grid blocks apart and forming 2x2 grid, change z value in matrix from 1 to 10


%Converting Arrays into 2D matrix
       for i=1:2
           for j=1:2
               m = (i-1)*Nx + j;
           end
       end

%Generating Coefficient Matrix 2x2 blocks
%Counting second row as a first row and first row as a second row 
%top left corner
for i = 2
    j = 1
    m = j+(i-1)*2;  %nth row of this i and j 
    ktop = 0;    %no flow at top
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    kright = harmonic_avg_fun(Perm(i,j),Perm(i,j+1));
    kbot = harmonic_avg_fun(Perm(i,j),Perm(i-1,j));
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m+1,m) = kright;
    A(m-2,m) = kbot;
end

%top right corner 
for i = 2
    j = 2
    m = j+(i-1)*2;  %nth row of this i and j
    ktop = 0; %no flow at top
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j-1));
    kright = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    kbot = harmonic_avg_fun(Perm(i,j),Perm(i-1,j));
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m-1,m) = kleft;
    A(m-2,m) = kbot;
end

%bottom right corner
for i = 1
    j = 2
    m = j+(i-1)*2;  %nth row of this i and j
    ktop = harmonic_avg_fun(Perm(i,j),Perm(i+1,j));
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j-1));
    kright = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    kbot = 0 %no flow at bottom
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m-1,m) = kleft;
    A(m+2,m) = ktop;
end

%bottom left corner
for i = 1
    j = 1
    m = j+(i-1)*2;  %nth row of this i and j
    ktop = harmonic_avg_fun(Perm(i,j),Perm(i+1,j));
    kright = harmonic_avg_fun(Perm(i,j),Perm(i,j+1));
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    kbot = 0 %no flow at bottom
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m+1,m) = kright;
    A(m+2,m) = ktop;
end

%Right hand side matrix B
for i = 1:2
    j = 1
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    m = j+(i-1)*2;
    B(m,1) = -kleft*p_in;
end

%Solving for A*phi = B
A_inv = inv(A);
phi_vec = A_inv*B;

%Convering soluion vector into 2D array
for i = 1:2
    for j = 1:2
         m = j+(i-1)*2;
         phi(i,j)= phi_vec(m);
    end
end
pphi = zeros(2,4);
pphi(1,:) = [100;64.5847;31.0586;0];
pphi(2,:) = [100;67.4755;33.9493;0];
phi=pphi;

%Flow coming in from left hand side and going out from right hand side
Q1 = ((harmonic_avg_fun(Perm(1,1),Perm(1,1))*(phi(1,1)-phi(1,2)))*Area)/(2*dx*dx); %Flow of first row
Q2 = ((harmonic_avg_fun(Perm(2,1),Perm(2,1))*(phi(2,1)-phi(2,2)))*Area)/(2*dx*dx); %Flow of second row
Q = Q1 + Q2; %Total of 2 flows 
Kb = (Q*(2*dx)*dx)/(Area*(p_in-p_out)); %Extracting permeability from total flow
Kf(I,J) = Kb; 
Kf_f = flip(Kf,2);
   end
end

PERM_X = Kf_f(:);

save('fb_kx_ly10.txt','PERM_X','-ASCII'); %producing txt file to create GRDECL

figure (1)
layer = 10; %layer can be changed here from 1 to 10
imagesc(rot90(log10(Kx(:,:,layer))));  
colormap('jet');
caxis([0 log10(20000)]);
axis image 
xlabel('Nx')
ylabel('Ny')
title('Fine Scale Permeability X Direction Layer 10')
colorbar
grid minor

figure (2)
layer = 1; %layer must stay the same, change in layer can be done in line 40
imagesc(rot90(log10(Kf(:,:,layer))));
colormap('jet');
caxis([0 log10(20000)]);
axis image 
xlabel('Nx')
ylabel('Ny')
title('Coarse Scale Permeability X Direction Layer 10')
colorbar
grid minor



       
       
       
       
