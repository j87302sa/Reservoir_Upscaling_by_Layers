%Attempt to upscale(ANALYTICAL METHOD)
%PETROLEUM DESIGN PROJECT

%STUDENT NAME - SAN HTOO AUNG
%ID - 10368431
%EMAIL - san.aung@student.manchester.ac.uk

clear;
clc;

%%POROSITY UPSCALING%%
po = load('PORO_5 - Copy.GRDECL');      %loading porosity from GRDECL file
po = reshape(po,120,60,10);     %Turning the vector into array
po = flip(po,2);        %Flipping

Nx = 120;       %Number of grid blocks in x-direction 
Ny = 60;        %Number of grid blocks in y-direction
Nz = 10;        %Number of grid blocks in z-direction 
Lx = 60;       %length of the reservoir model in x-direction  
Ly = 60;        %length of the reservoir model in y-direction
layer = 1:Nz;   %Layering of z

%% UPSCALING POROSITY USING ARITHEMATIC AVERAGING METHOD

%Upscaling of porosity in all 10 layers
for i = 1:1:60
    for j = 1:1:30
        for layer = 1:1:Nz
            I = (2*i)-1;
            J = (2*j)-1;
            K = layer;
            upscl_po(i,j,K)  = (po(I,J,K)+po(I,J+1,K)+po(I+1,J,K)+po(I+1,J+1,K))/4;
        end
    end
end

upscl_po_new = flip(upscl_po,2);
PORO = upscl_po_new(:);  %Assigning all rows and columns into one column to include as GRDECL in RE studio
save('upscl_po_NEW.txt','PORO','-ASCII'); %producing txt file to create GRDECL

%% Figures
k = input('Layer Number:'); %to produce each layer 

figure(11)
layer = k;
imagesc((rot90(po(:,:,layer))));
colormap('jet');
axis image 
xlabel('Nx')
ylabel('Ny')
title('FINE SCALE POROSITY LAYER 10')
colorbar
grid minor

figure(12)
layer = k;
imagesc((rot90(upscl_po(:,:,layer))));
colormap('jet');
axis image
xlabel('Nx')
ylabel('Ny')
title('COARSE SCALE POROSITY LAYER 10')
colorbar
grid minor



            
            
