%Attempt to upscale(FLOW BASED UPSCALING METHOD)
%PETROLEUM DESIGN PROJECT

%STUDENT NAME - SAN HTOO AUNG
%ID - 10368431
%EMAIL - san.aung@student.manchester.ac.uk

clear;
clc;
Ky = load('PERMY_5 - Copy.GRDECL');        %loading y-perm from GRDECL file
Ky = reshape(Ky,120,60,10);     %Turning the vector into array
Ky = flip(Ky,2);       %Flipping
Ky_z = Ky(:,:,1);  %layer number

p_in = 100;      %Assigning inlet Pressure to 100 bar
p_out = 0;      %Assigning outlet Pressure to 90 bar
mew = 0.05;       %Assigning arbitrary viscosity

Nx = 60;       %Number of grid blocks in x-direction 
Ny = 120;        %Number of grid blocks in y-direction
Nz = 10;        %Number of layer in z-direction 
Lx = 60;       %length of the reservoir model in x-direction  
Ly = 60;        %length of the reservoir model in y-direction
layer = 1:Nz;   %Layering of z

%GRIDDING
dx = Lx/Nx; % grid size in x-direction for x direction flow
dy = Ly/Ny; % grid size in y-direction for y direction flow
Area = (2*dx)*(2*dy);    %Area of a reservoir
Kf = zeros(60,30);      %targeted upscaling parameters

%% UPSCALING FOR FLOW BASED MODEL [Y]

%Upscaling in a way that 2x2 blocks into 1 single block
for I = 1:Ny/2
   for J = 1:Nx/2
       N = 2*2;              % number of unknowns
       A = zeros(N,N);       %N rows and N columns for coefficient matrix
       B = zeros(N,1);       %N rows and 1 column for the right hand side vector
       Perm = Ky(2*I-1:2*I,2*J-1:2*J,1); %rows and columns are 2 grid blocks apart and forming 2x2 grid, change z value in matrix from 1 to 10


%Converting Arrays into 2D matrix
 for i=1:2
           for j=1:2
               m = (i-1)*Nx + j;
           end
       end
       
%Generating Coefficient Matrix 2x2 blocks
%Counting second row as a first row and first row as a second row 
%bottom left corner
for i = 1
    j = 1
        m = j+(i-1)*2;  %nth row of this i and j
        ktop = harmonic_avg_fun(Perm(i,j),Perm(i+1,j));
        kleft = 0 %no flow at left
        kright = harmonic_avg_fun(Perm(i,j),Perm(i,j+1));
        kbot = harmonic_avg_fun(Perm(i,j),Perm(i,j));
        A(m,m) = -(ktop+kright+kleft+kbot);
        A(m+1,m) = kright;
        A(m+2,m) = ktop;
end

%bottom right corner
for i = 1
    j = 2
    m = j+(i-1)*2;  %nth row of this i and j
    ktop = harmonic_avg_fun(Perm(i,j),Perm(i+1,j));
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j-1));
    kright = 0; %no flow at right
    kbot = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m+2,m) = ktop;
    A(m-1,m) = kleft;
end

%top right corner
for i = 2
    j = 2
    m = j+(i-1)*2;  %nth row of this i and j
    ktop = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    kleft = harmonic_avg_fun(Perm(i,j),Perm(i,j-1));
    kright = 0 %no flow at right 
    kbot = harmonic_avg_fun(Perm(i,j),Perm(i-1,j));
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m-2,m) = kbot;
    A(m-1,m) = kleft;
end

%top left corner
for i = 2
    j = 1
    m = j+(i-1)*2;  %nth row of this i and j
    ktop = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    kleft = 0 %no flow at left
    kright = harmonic_avg_fun(Perm(i,j),Perm(i,j+1));
    kbot = harmonic_avg_fun(Perm(i,j),Perm(i-1,j));
    A(m,m) = -(ktop+kright+kleft+kbot);
    A(m+1,m) = kright;
    A(m-2,m) = kbot;
end

%Right hand side matrix B
for j = 1:2
    i = 1
    ktop = harmonic_avg_fun(Perm(i,j),Perm(i,j));
    m = j+(i-1)*2;
    B(m,1) = -ktop*p_in;
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

pphi(2,:)=phi(1,:);
pphi(3,:)=phi(2,:);
pphi(1,:)=100;
pphi(4,:)=0;
phi=pphi;

%Flow coming in from top and going out from bottom
Q1 = ((harmonic_avg_fun(Perm(2,1),Perm(2,1))*(phi(1,1)-phi(2,1)))*Area)/(2*dy*dy); %Flow of first row
Q2 = ((harmonic_avg_fun(Perm(2,2),Perm(2,2))*(phi(1,2)-phi(2,2)))*Area)/(2*dy*dy); %Flow of second row
Q = Q1 + Q2; %Total of 2 flows 
Kb = (Q*(2*dy)*dy)/((4*dx*dy)*(p_in-p_out)); %Extracting permeability from total flow
Kf(I,J) = Kb;
Kff_y = flip(Kf,2);
   end
end

PERM_Y = Kff_y(:);
save('fb_ky_ly1.txt','PERM_Y','-ASCII'); %producing txt file to create GRDECL

figure (1)
layer = 1; %layer can be changed here from 1 to 10
imagesc(rot90(log10(Ky(:,:,layer))));  
colormap('jet');
caxis([0 log10(20000)]);
axis image 
xlabel('Nx')
ylabel('Ny')
title('Fine Scale Permeability Y Direction Layer 1')
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
title('Coarse Scale Permeability y Direction Layer 1')
colorbar
grid minor

