
% Read the file
data = dlmread('data_6302.txt');
z = data(:,1);
w = data(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT method1 to find Frequency
[fhz,fft_spectrum] = do_fft(z(:),w(:));
[amplitude1, idx] = max(fft_spectrum);
frequency = fhz(idx);
period = 1/frequency;

figure(2)
plot(fhz, fft_spectrum)
hold on
xlabel('Frequency')
ylabel('Amplitude')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % FFT method2 to find Frequency
% L = length(z);
% fs = 1/mean(diff(z));                                 % Sampling Frequency
% fn = fs/2;                                              % Nyquist Frequency
% yc = w - mean(w);                                       % Sampling Interval 
% FTv = fft(yc)/L;                                        % Fourier Transform
% Fv = linspace(0, 1, fix(L/2)+1)*fn;                     % Frequency Vector
% Iv = 1:length(Fv);                                      % Index Vector
% [max_FTv, maxidx] = max(FTv(Iv));
% Fv_at_max_FTv = Fv(maxidx);
% wavelength1=1/Fv_at_max_FTv;
% 
% figure(3)
% plot(Fv, abs(FTv(Iv))*2)
% grid
% hold on
% xlabel('Frequency')
% ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zero crossing method -alternative for clean periodic signals-
zct = find_zc(z,w,0);
wavelength2 = mean(diff(zct));
amplitude2 = 0.5*(max(w) - min(w));

figure(4)
plot(z,w,'g',zct,zeros(1,numel(zct)),'*r','markersize',25)
xlabel('Spanwise length')
ylabel('Velocity')
legend('sinosodial shape','wavelength points')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zct = find_zc(x,w,threshold)

% positive slope "zero" crossing detection, using linear interpolation
w = w - threshold;
zci = @(data) find(diff(sign(data))>0);  %define function: returns indices of +ZCs
ix=zci(w);                      %find indices of + zero crossings of x
ZeroX = @(x0,y0,x1,y1) x0 - (y0.*(x0 - x1))./(y0 - y1); % Interpolated x value for Zero-Crossing 
zct = ZeroX(x(ix),w(ix),x(ix+1),w(ix+1));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq_vector,fft_spectrum] = do_fft(time,data)

dt = mean(diff(time));
Fs = 1/dt;
nfft = length(data); % maximise freq resolution => nfft equals signal length

% window : hanning
window = hanning(nfft);
cor_coef = length(window)/sum(window);

% fft scaling
% fft_spectrum = abs(fft(data))/nfft;
fft_spectrum = abs(fft(data.*window))*2*cor_coef/nfft;

% one sidded fft spectrum  % Select first half 
    if rem(nfft,2)    % nfft odd
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)';
    end
fft_spectrum = fft_spectrum(select,:);
freq_vector = (select - 1)*Fs/nfft;

end
