function yn = hpf(yo)
% High-frequency Pass Filter at the same domain of y
%
% 05-Aug-2013 18:08:10
% hanligong@gmail.com

N = length(yo);
%Correction function
f_corr = -2*abs((1:N)-(N/2+1))/N+1; % ramp filter
%Fourier transform
f = fft(yo);
%Enhance high-frequency
f = f.*f_corr;
%Inverse Fourier transform
yn = ifft(f);
yn = real(yn);

% N = length(yo);
% yo(end+1:end+N) = 0; N = 2*N;
% %Correction function
% f_corr = -2*abs((1:N)-(N/2+1))/N+1; % ramp filter
% %Fourier transform
% f = fft(yo);
% %Enhance high-frequency
% f = f.*f_corr;
% %Inverse Fourier transform
% yn = ifft(f);
% yn = real(yn);
% yn = yn(1:N/2);
end
