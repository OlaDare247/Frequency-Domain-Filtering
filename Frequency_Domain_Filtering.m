
% Exercise 1: Try the following: in the above process skip the padding
% step and perform filtering. Observe any differences.
% Load the Image
A = imread('lena_gray_256.png');
figure, imagesc(A), axis image, colormap gray, colorbar;
title('Original Image');
% Compute Fourier Transfrom without padding
FA = fftshift(fft2(A));
figure, imagesc(log(1+abs(FA))), axis image, colormap gray, colorbar;
title('Magnitude Spectrum (Without Padding)');
% Create Ideal Low Pass Filter (Same Size as A)
[x, y] = meshgrid(-128:127, -128:127); % 256x256 meshgrid
z = sqrt(x.^2 + y.^2);
ILPF = (z < 15);
figure, imagesc(ILPF), axis image, colormap gray, colorbar;
title('Ideal Low_Pass Filter (Cutoff = 15)');
% Apply Low-Pass Filter in Frequency Domain
FFA = FA .* ILPF;
figure, imagesc(log(1+abs(FFA))), axis image, colormap gray, colorbar;
title('Filtered Spectrum (Without Padding)');
% Apply Inverse Fourier Transform
IFFA = ifft2(ifftshift(FFA));
IFFA = real(IFFA); % We are considering the real part only
figure, imagesc(IFFA), axis image, colormap gray, colorbar;
title('Filtered Image Without Padding)');
% Exercise 2: Please apply the following processes.
% 1. Open the test image house.
% 2. Generate ideal lowpass filters with radii 30, 50, 100.
% 3. Apply the filter and show results.
% 4. Make observations on the amount of blurring and presence of ringing artifacts.
% Let us clear the workspace and close all figures
clear; clc; close all;
% 1. We will load in the test image 'house'
A = imread('house.png'); 
A = im2double(A);
figure, imagesc(A), axis image, colormap gray, colorbar;
title('Original House Image');
% Compute Fourier Transform
FA = fftshift(fft2(A)); % Compute centered DFT
figure, imagesc(log(1+abs(FA))), axis image, colormap gray, colorbar;
title('Magnitude Spectrum of Original Image');
% Get image size
[M, N] = size(A);
% 3. Generate Ideal Low-Pass Filters with Radii 20, 50 and 100
[x, y] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
% Define the frequency domain distance from center
z = sqrt(x.^2 + y.^2);
% Create filters
ILPF_30 = (z < 30);
ILPF_50 = (z < 50);
ILPF_100 = (z < 100);
% Display the filters
figure;
subplot(1,3,1), imagesc(ILPF_30), axis image, colormap gray, colorbar, title('Ideal LPF (Radius = 30');
subplot(1,3,2), imagesc(ILPF_50), axis image, colormap gray, colorbar, title('Ideal LPF (Radius = 50');
subplot(1,3,3), imagesc(ILPF_100), axis image, colormap gray, colorbar, title('Ideal LPF (Radius = 100');
% 4. Apply the filters in the frequency domain
FA_30 = FA .* ILPF_30;
FA_50 = FA .* ILPF_50;
FA_100 = FA .* ILPF_100;
% Compute inverse DFT to obtain filtered images
IFFA_30 = real(ifft2(ifftshift(FA_30)));
IFFA_50 = real(ifft2(ifftshift(FA_50)));
IFFA_100 = real(ifft2(ifftshift(FA_100)));
% Display filtered images
figure;
subplot(1,3,1), imagesc(IFFA_30), axis image, colormap gray, colorbar, title('Filtered Image (R=30)');
subplot(1,3,2), imagesc(IFFA_50), axis image, colormap gray, colorbar, title('Filtered Image (R=50)');

subplot(1,3,3), imagesc(IFFA_100), axis image, colormap gray, colorbar, title('Filtered Image (R=100)');
% Exercise 3 (10 pt). Repeat above process for D0 = 30, 50, 100. What do you observe?
% Clear workspace and close figures
clear; clc; close all;
% 1. Load the test image 'Lena'
A = imread('lena_gray_256.png');
% Normalize the pixel values
A = im2double(A); 
% Make a subplot
figure, imagesc(A), axis image, colormap gray, colorbar;
title('Original Image');
% 2. Compute Fourier Transform
FA = fftshift(fft2(A)); % Compute centered DFT
figure, imagesc(log(1+abs(FA))), axis image, colormap gray, colorbar;
title('Magnitude Spectrum of Original Image');
% Get image size
[M,N] = size(A);
% 3. Generate Butterworth Low-Pass Filters with D0 = 30, 50, 100
[x, y] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
n = 1; % Order of Butterworth filter
D0_values = [30, 50, 100]; % Different cutoff frequencies
BLPFs = cell(1, length(D0_values));
for i = 1:length(D0_values)
 
D0 = D0_values(i);
 
BLPFs{i} = 1./ (1+(sqrt(x.^2 + y.^2) / D0).^(2 * n));
end
% Display the filters
figure;
for i = 1:length(D0_values)
 
subplot(1,3,1), imagesc(BLPFs{i}), axis image, colormap gray, colorbar;
 
title(['Butterworth LPF (D0 = ', num2str(D0_values(i)), ')']);
end
% Apply the filters in the frequency domain
FA_filtered = cell(1, length(D0_values));
IFFA_filtered = cell(1, length(D0_values));
for i = 1:length(D0_values)
 
FA_filtered{i} = FA .* BLPFs{i};
 
IFFA_filtered{i} = real(ifft2(ifftshift(FA_filtered{i}))); % Inverse FFT
end
% Display filtered images
figure;
for i = 1:length(D0_values)
 
subplot(1,3,i), imagesc(IFFA_filtered{i}), axis image, colormap gray, colorbar;
 
title(['Filtered Image (D0 = ', num2str(D0_values(i)), ')']);
end
% Exercise 4 (10 pt). Apply the above bandpass filter for multiple center frequencies and bandwidths to extract different characteristics from 
the image.
% clear workspace and close figures
clear; clc; close all;
% 1. load the test image 'lena'
A = imread('lena_gray_256.png');
A = im2double(A), axis image, colormap gray, colorbar;
title('Original Image');
% 2. Compute Fourier Transform
FA = fftshift(fft2(A)); % Compute centered DFT
figure, imagesc(log(1+abs(FA))), axis image, colormap gray, colorbar;
title('Magnitude Spectrum of Original Image');
% Get image size
[M,N] = size(A);
% Generate Gaussain Band-Pass Filters For Multiple D0 and W values
[x,y] = meshgrid(-M/2:M/2-1, -N/2:N/2-1);
D = sqrt(x.^2 + y.^2);
% Define multiple center frequncies and bandwidths
D0_values = [20, 50, 100]; % Center frequencies
W_values = [10, 30, 50]; % Bandwidths

