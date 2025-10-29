clc
clear
close all;

%% First part (Execution time of DFT and FFT)
L = 8192;
x = rand(1, L);        % Xi[n] random signal

% Measure execution time for DFT {x[n]}
tic;
X_dft = DFT(x);
dft_time = toc;
fprintf('Execution Time for DFT = %.4f seconds \n', dft_time);

% Measure execution time for FFT {x[n]}
tic;
X_fft = fft(x);        % Calculate fast fourier transform directly using built-in function fft()
fft_time = toc;
fprintf('Execution Time for FFT = %.4f seconds\n', fft_time);

%% DFT implementation Function
function X = DFT(x)
    N = length(x);
    X = zeros(1, N); 
    for k = 0:N-1
        for n = 0:N-1
            X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * n * k / N);
        end
    end
end