%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: For_loop_SI-DHM_demodulation                                          %
%                                                                              %
% The script cperforms the demodulation and composition of SR images from      %
% SI-DHM hologram witha random phase shift between them                        %
%                                                                              %                                       
% Authors: Sofia Obando, Raul Castaneda, Carlos Trujillo, Rene Restrepo,       %
%          Ana Doblas.                                                         %
% Applied Optics Group EAFIT univeristy, Colombia                              % 
% and Optica Imaging Research Laboratory University of MAssachusetts Dartmouth % 
%                                                                              %
% Email: racastaneq@eafit.edu.co; adoblas@umassd.edu                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%% Input parameters to record the holograms
lambda = 0.532;   % Source's wavelength in microns
dxy = 3.75;       % Pitch in X and Y direction in microns
region = 2;  % Region for filtering (quadrant follows the Cartession convention)
step = 0.025;  % set the step for the For Loop searc, in the range = 0 to 2Pi

% Read recorded holograms
path_1 = 'C:\Users\racastaneq\Documents\MEGA\MEGAsync\RACQ\Universities\05 EAFIT\Research projects\2023\SI-SOFI\GitHub versions\'; % Directory for read holograms
name1 =  'holo1.tif'; 
name2 =  'holo2.tif'; 

h1 = function_heuristic.holo_read(strcat(path_1, name1)); % Read first hologram
h2 = function_heuristic.holo_read(strcat(path_1, name2)); % Read second hologram

% Display the first hologram
figure,imagesc(h1),colormap(gray),colorbar,title('Hologram h1'),daspect([1 1 1])

% output
out = function_heuristic.SIDHM(h1, h2, lambda, dxy);
figure, imagesc(abs(out).^0.2), colormap gray, title('FT SIM')
out = fftshift(ifft2(fftshift(out)));    
figure, imagesc(angle(out)), colormap gray, title('Phase SI'),daspect([1 1 1])






