%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: functions_evaluation                                                  %
%                                                                              %
% The script contains all implemented function for SI_DHM_heuristic            %
%                                                                              %                                       
% Authors: Raul Castaneda, Sofia Obando, Carlos Trujillo, Rene Restrepo,       %
%           Ana Doblas.                                                        %
% Applied Optics Group EAFIT univeristy                                        % 
%                                                                              %
% Email: racastaneq@eafit.edu.co; adoblas@umassd.edu                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef function_heuristic
    methods(Static)

        function [holo,M,N,m,n] = holo_read(filename)
            % function for loading holograms
            % INPUTUS..
            % filename:XXXXXXXXXXXXXXXXXXXXXXXXXX
            holo = double(imread(filename));
            holo = holo(:,:,1);
            [N,M] = size(holo);
            [n,m] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
        end


        function [holo_filter,holo_FT,fx_max,fy_max] = spatialFilter_SIDHM(holo,M,N,visual,factor)
            % function to performance the spatial filter over a SI-DHM
            % INPUTUS..
            % holo:XXXXXXXXXXXXXXXXXXXXXXXXXX
            % M,N:xxxxxxxxxxxxxxxxxxxxxx
            % visual:xxxxxxxxxxxxxxxxxxxxxxxxxx
            % factor:xxxxxxxxxxxxxxxxxxxxxxxxxxx
            fft_holo = fftshift(fft2(fftshift(holo(:,:,1)))); 
            figure,imagesc(log(abs(fft_holo).^2)),colormap(gray),title('FT Hologram'),daspect([1 1 1]) 
            
            [pointsY,pointsX] = ginput(2);
            x3 = round(pointsX(1));
            x4 = round(pointsX(2));
            y3 = round(pointsY(1));
            y4 = round(pointsY(2));
           
            mask = zeros(M,N);
            mask(x3:x4,y3:y4) = 1;
            fft_holo_I = fft_holo .* mask;
            %figure,imagesc((abs( fft_holo_I).^0.1)),colormap(gray),title('FT Hologram'),daspect([1 1 1]) 
            
            % max values first peak
            maxValue_1 = max(max(abs(fft_holo_I)));
            [fy_max_1 fx_max_1] = find(abs(fft_holo_I) == maxValue_1);
            mask(fy_max_1 - 20:fy_max_1 + 20,fx_max_1 - 20:fx_max_1 + 20)=0;
            fx_max_L = fx_max_1;
            fy_max_L = fy_max_1;
            fft_holo_I = fft_holo_I .* mask;
            %figure,imagesc(log(abs(fft_holo_I).^2)),colormap(gray),title('FT FIrst peak'),daspect([1 1 1])
            
            maxValue_1 = max(max(abs(fft_holo_I)));
            [fy_max_1 fx_max_1] = find(abs(fft_holo_I) == maxValue_1);
            fx_max_D = fx_max_1(1);
            fy_max_D = fy_max_1(1);
        
            fy_max = [fx_max_L,fx_max_D];
            fx_max = [fy_max_L,fy_max_D];
        
            %find the centers between both peaks 
            middlePoint_X = (fx_max_D + fx_max_L) / 2;
            middlePoint_Y = (fy_max_D + fy_max_L) / 2;
            
                    
            mask = zeros(M,N);
            mask(x3:x4,y3:y4) = 1;
            
            fft_filter_holo = fft_holo .* mask;
            holo_filter_1 = fftshift(ifft2(fftshift(fft_filter_holo)));
        
            fft_holo_2 = fftshift(fft2(fftshift(holo(:,:,2)))); 
            %fft_filter_holo_2 = fft_holo_2 .* filter;
            fft_filter_holo_2 = fft_holo_2 .* mask;
            holo_filter_2 = fftshift(ifft2(fftshift(fft_filter_holo_2)));
        
        
            holo_filter(:,:,1) = holo_filter_1;
            holo_filter(:,:,2) = holo_filter_2;
        
            holo_FT(:,:,1) = fft_filter_holo;
            holo_FT(:,:,2) = fft_filter_holo_2;
        
            if visual == 'Yes'
                figure,imagesc(log(abs(fft_filter_holo).^2)),colormap(gray),title('FT Filter Hologram'),daspect([1 1 1]) 
                figure,imagesc((abs(holo_filter_1).^2)),colormap(gray),title('Filter Hologram'),daspect([1 1 1]) 
            end
        end


        function [cf] = costFunction_SIDHM(theta, FTHolo, fx_max,fy_max)
            cf = 0;
            [M,N] = size(FTHolo);
            [Dtemp] = function_heuristic.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,1);
            cf =  abs(Dplus(fx_max(1),fy_max(1))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end


        function [cf] = costFunction_SIDHMII(theta, FTHolo, fx_max,fy_max)
            cf = 0;
            [M,N] = size(FTHolo);
            [Dtemp] = function_heuristic.demComp2SIDHM(theta,FTHolo);
            Dplus = Dtemp(:,:,2);
            cf =  abs(Dplus(fx_max(2),fy_max(2))) / (abs(Dplus(fx_max(1),fy_max(1))) + abs(Dplus(fx_max(2),fy_max(2))));
        end


        function [D] = demComp2SIDHM(theta, H)
            [X, Y, no] = size(H); 
            D = zeros(X,Y,no);
            M = 1/2*[exp(1i*0) exp(-1i*0);exp(1i*theta) exp(-1i*theta)];
            Minv = pinv(M);
            
            D(:,:,1) = Minv(1,1).*H(:,:,1) + Minv(1,2).*H(:,:,2);
            D(:,:,2) = Minv(2,1).*H(:,:,1) + Minv(2,2).*H(:,:,2);
        end


        function  ref = phase_rec(filename, dx, dy, lambda, region,  save)
            % Main function to perform phase reconstruction
            holo = filename;
    
            % Get the size of the hologram
            [N,M] = size(holo);
            % Create a meshgrid for the hologram
            [m,n] = meshgrid(-M/2:M/2-1,-N/2:N/2-1);
            
            % Calculate the Fourier Transform of the hologram and shift the zero-frequency component to the center
            ft_holo = fftshift(fft2(fftshift(holo)));
            
            
            % Initialize a filter with zeros
            filter = zeros(N,M);
            
            if region==1
                filter(1:round(N/2-(N*0.1)),round(M/2+(M*0.1)):M) = 1; % 1nd quadrant
            elseif region==2
                filter(1:round(N/2-N*0.1),1:round(M/2 -M*0.05)) = 1;  % 2nd quadrant
            elseif region==3
                filter(round(N/2+(N*0.1)):N,1:round(M/2-(M*0.1))) = 1; % 3nd quadrant
            else
                filter(round(N/2+(N*0.1)):N,round(M/2+(M*0.1)):M) = 1; % 4nd quadrant
            end
    
            % Apply the filter to the Fourier Transform of the hologram
            ft_filtered_holo = ft_holo .* filter;
            
            filtered_spect = log(abs(ft_filtered_holo).^2);
            % Find the maximum value in the filtered spectrum
            [~,idx] = max(filtered_spect(:));
            
            % Define wavenumber
            k = 2 * pi / lambda;
            
            % Calculate the center frequencies for fx and fy
            fx_0 = M/2;
            fy_0 = N/2;
            
            % Get the maximum values of fx and fy
            [fy_max,fx_max] = ind2sub([N,M],idx);
            
            % Define the step size for the search
            step = 0.9;
            
            % Initialize variables for the search
            j = 0;
            
            % Calculate the Inverse Fourier Transform of the filtered hologram
            holo_rec = fftshift(ifft2(fftshift(ft_filtered_holo)));
            
            % Define the search range (G)
            G = 3;
            % Initialize flag for the search loop
            fin = 0;
            
            % Set initial values for fx and fy
            fx = fx_max;
            fy = fy_max;
            
            % Initialize temporary search range
            G_temp = G;
            
            % Loop to find the optimal fx and fy values
            while fin == 0
              i = 0;
              j = j + 1;
              
              % Initialize the maximum sum (for thresholding)
              suma_maxima = 0;
              
              % Nested loops for searching in the range of fx and fy
              for fy_tmp = fy-step*G_temp:step:fy+step*G_temp
                for  fx_tmp = fx-step*G_temp:step:fx+step*G_temp
                  i = i + 1;
                  
                  % Calculate the metric for the current fx and fy
                  [suma] = function_heuristic.metric(holo_rec, fx_0, fy_0, fx_tmp, fy_tmp, lambda, M, N , dx, dy, m, n, k);
                  
                  % Update maximum sum and corresponding fx and fy if 
                  % current sum is greater than the previous maximum
                  if (suma > suma_maxima)
                    x_max_out = fx_tmp;
                    y_max_out = fy_tmp;
                    suma_maxima=suma;
                  end         
                end
              end
              
              % Update the temporary search range
              G_temp = G_temp - 1;
              
              % Check if the optimal values are found, set the flag to exit the loop
              if x_max_out == fx && y_max_out == fy
                fin = 1;        
              end
              
              % Update fx and fy for the next iteration
              fx = x_max_out;
              fy = y_max_out; 
            end

            % Calculate the angles for the compensation wave
            theta_x = asin((fx_0 - x_max_out) * lambda / (M * dx));
            theta_y = asin((fy_0 - y_max_out) * lambda / (N * dy));
            
            % Calculate the reference wave
            ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));
            
            % Apply the reference wave to the hologram reconstruction
            holo_rec2 = holo_rec .* (ref);
            Ampl = abs(holo_rec2);
            phase = angle(holo_rec2);

        end

        
        function suma = metric(holo_rec, fx_0, fy_0, fx_tmp, fy_tmp, lambda, M, N , dx, dy, m, n, k)
            % Function to calculate the metric for the current fx and fy
            % Calculate the angles for the compensation wave
            theta_x = asin((fx_0 - fx_tmp) * lambda / (M * dx));
            theta_y = asin((fy_0 - fy_tmp) * lambda / (N * dy));
    
            % Calculate the reference wave
            ref = exp(1i * k * (sin(theta_x) * m * dx + sin(theta_y) * n * dy));
    
            % Apply the reference wave to the hologram reconstruction
            holo_rec2 = holo_rec .* ref;
                
            % Calculate the phase of the hologram reconstruction
            phase = angle(holo_rec2);
        
            % Normalize the phase and convert it to uint8
            phase = mat2gray(phase);
            phase = uint8(phase * 255);
        
            % Threshold the phase image
            BW = imbinarize(phase, 0.1);
        
            % Calculate the sum of all elements in the resulting binary image
            suma = sum(sum(BW));
        end

        function output = SIDHM(h1, h2, lambda, dxy)
            % Sometimes a shift in the phase is necessary to avercome some phase jumps.
            % If the viasualization of the phase demodulation has any phase jump, change
            % the value of these two variables (only 1 or -1 values are possible)
            PhaseOrientation_Plus = -1;
            PhaseOrientation_Minus = 1;

            % Read dimensions of the holograms
            [N, M] = size(h1);

            % Create a meshgrid for the hologram
            [m, n] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
            
            holo(:,:,1) = h1;
            holo(:,:,2) = h2;

            % Lines to implement the spatial filter using a circular mask
            [holo_filter,holo_FT,fx_max,fy_max] = function_heuristic.spatialFilter_SIDHM(holo,M,N,'Not',5);

            % minimization particleswarm
            % Blind demodulation to recover the two shifted object spectrum
            % individually without prior knowledge of the phase in the SI
            % pattern Gplus
            test_theta = randi([0 360])*pi/180;
            lb = 0; 
            ub = 360*pi/180; 

            options = optimoptions('particleswarm', 'Display', 'off', 'SwarmSize', 4, 'MaxIterations', 100, 'FunctionTolerance', 1e-6);
            cf_particleswarm = function_heuristic.costFunction_SIDHM(test_theta,holo_FT,fx_max,fy_max);
            [theta, cf_particleswarm] = particleswarm(@(params) function_heuristic.costFunction_SIDHM(params,holo_FT,fx_max,fy_max), 1, lb, ub, options);
            [Gdemod] =  function_heuristic.demComp2SIDHM(theta, holo_FT);
            Gplus_demod = Gdemod(:,:,1);

            gplus = fftshift(ifft2(fftshift(Gplus_demod)));
            gplus = gplus * PhaseOrientation_Plus;
            ref_gplus = function_heuristic.phase_rec(gplus, dxy, dxy, lambda, 2, false);
            Gplus = gplus .* ref_gplus;
            phasePlus = angle(Gplus);
            figure,imagesc(phasePlus),colormap(gray),colorbar,title('Phase gplus'),daspect([1 1 1])

            % minimization particleswarm
            % Blind demodulation to recover the two shifted object spectrum
            % individually without prior knowledge of the phase in the SI
            % pattern Gplus
            cf_particleswarm = function_heuristic.costFunction_SIDHMII(test_theta,holo_FT,fx_max,fy_max);
            [theta, cf_particleswarm] = particleswarm(@(params) function_heuristic.costFunction_SIDHMII(params,holo_FT,fx_max,fy_max), 1, lb, ub, options);
            [Gdemod] =  function_heuristic.demComp2SIDHM(theta, holo_FT);
            Gminus_demod = Gdemod(:,:,2);
            
            gminus = fftshift(ifft2(fftshift(Gminus_demod)));
            gminus = gminus * PhaseOrientation_Minus;
            ref_minus = function_heuristic.phase_rec(gminus, dxy, dxy, lambda, 2, false);
            Gminus = gminus .*ref_minus;
            phaseMinus = angle(Gminus);
            figure,imagesc(phaseMinus),colormap(gray),colorbar,title('Phase gminus'),daspect([1 1 1])

            % Combine and display the results
            Ft_Gminus = fftshift(fft2(fftshift(Gminus)));
            Ft_Gminus_norm = (Ft_Gminus - min(min(Ft_Gminus)))/(max(max(Ft_Gminus))+ min(min(Ft_Gminus)));
            
            Ft_Gplus = fftshift(fft2(fftshift(Gplus)));
            Ft_Gplus_norm = (Ft_Gplus - min(min(Ft_Gplus)))/(max(max(Ft_Gplus))+ min(min(Ft_Gplus)));
            
            Gsim = Ft_Gminus_norm + Ft_Gplus_norm;
            output = Gsim;

        end

    end
end
