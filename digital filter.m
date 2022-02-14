close all;
clear;
clc;

frequency_band_filters = ['0 - 170','170 - 310','310 - 600','600 - 1000','1k - 3k','3k - 6k','6k - 12k','12k - 14k','14k - 16k'];
frequency_bands = [0 170 310 600 1000 3000 6000 12000 14000 16000];
bandgains = [];

% Getting audio file name from the user
% audiofile = input('Please enter the file name:  ','s');
[y ,Fs] = audioread('sample_mono.wav');
N = length(y);
t = linspace(0, N/Fs, N);
Fm = Fs/2;

%Prompting gains to be used in each frequency band
prompt = ('Please enter the gain in the bands ');
for i = 1:9
    fprintf('Please enter the gain in the bands ');
    fprintf('%d - %d',frequency_bands(i),frequency_bands(i+1));
    x = input('Hz: ');
    bandgains = [bandgains,x];


end
% for i = 1:9
%     chr = convertStringsToChars(frequency_bands(i));
%     current =  input([prompt chr ' Hz in dB: ']);
%     current_in_watt = 10^(current/20);
%     bandgains = [bandgains current_in_watt];
%  end

% for i = 1:9
%     prompt(frequency_bands(i));
% end

% bandgains = [1 0 0 0 0 0 0 0 0];

%Getting and validating type of filter FIR or IIR
filter_type = input([newline 'Choose the type of filter you want (type 1 or 2)' newline '1. FIR' newline '2. IIR' newline]);
 while(filter_type ~= 2 && filter_type ~= 1)
     filter_type = input(['Choose a valid type of filter (type 1 or 2)' newline]);
 end

if filter_type == 1
    filter_order  = 40;
else
    filter_order = 4;
end

%Checking if user wants to change output sample frequency or not, and if
%yes asking the user for that frequency
% new_sample_rate = input(['Do you want a specific sample rate for the output? or leave it as it is?' newline '1. Add a new sample rate' newline '2. Leave it as it is' newline]);
% while(new_sample_rate ~= 2 && new_sample_rate ~= 1)
%     new_sample_rate = input(['Please choose a valid option' newline '1. Add a new sample rate' newline '2. Leave it as it is' newline]);
% end
new_Fs = 1*Fs;
% if(new_sample_rate == 1)
%     factor = input('Please enter the factor to be multiplied by the original sample rate: ');
%     new_Fs = factor*Fs;
% end
yt = 0;
%IIR
if filter_type == 2
    for i = 1:9
        if i == 1
            Wc = frequency_bands(i+1);
            Wd = Wc/Fm; %normalised form
            if (Wd > 0 && Wd < 1 )
                [num ,den] = butter(filter_order,Wd,'low');
            end
        else
            Wc = frequency_bands(i:i+1);
            Wd = Wc/Fm; %normalised form
            if ((Wd(1)>0) && (Wd(1)<1) && (Wd(2)>0) && (Wd(2)<1))
                [num,den] = butter(filter_order,Wd,'bandpass');
            end
        end
        y2=filter(num,den,y);
        figure; subplot(3,2,1);
        stepz(num,den); %plotting step response
        subplot(3,2,3);
        impz(num,den); %plotting impulse response
        [H, f] = freqz(num,den,512,Fs);
        mag = 20*log10(abs(H)); %getting magnitude of H
        phase = angle(H)*180/pi; %getting phase of H
        subplot(3,2,2);
        plot(f*180/pi,mag); %plotting magnitude of H
        xlabel('Frequency (Hz)');
        ylabel('|H(f)|');
        title('Magnitude');
        subplot(3,2,4);
        plot(f*180/pi,phase); %plotting phase of H
        xlabel('Frequency (Hz)');
        ylabel('<H(f) (degrees)');
        title('Phase');
        subplot(3,2,5);
        zplane(num,den); %plotting zeros and poles of the system
        title('Zeros and Poles');
        subplot(3,2,6);
        plot(t,y2); %plotting the signal after filteration in time domain
        title('Gain');
        xlabel('time (s)');
        ylabel('Amplitude');
        yt = yt + (y2*bandgains(i));
    end

%FIR
else
    for i = 1:9
        if i == 1
            Wc = frequency_bands(i+1);
            Wd = Wc/Fm; %normalised form
            if (Wd > 0 && Wd < 1 )
                [num, den] = fir1(filter_order,Wd,'low');
            end
        else
            Wc = frequency_bands(i:i+1);
            Wd = Wc/Fm; %normalised form
            if ((Wd(1)>0) && (Wd(1)<1) && (Wd(2)>0) && (Wd(2)<1))
                [num,den]=fir1(filter_order,Wd,'bandpass');
            end
        end
        y2=filter(num,den,y);
        figure; subplot(3,2,1);
        stepz(num,den); %plotting step response
        subplot(3,2,3);
        impz(num,den); %plotting impulse response
        [H, f] = freqz(num,1,512,Fs);
        mag = 20*log10(abs(H)); %getting magnitude of H
        phase = angle(H)*180/pi; %getting phase of H
        subplot(3,2,2);
        plot(f*180/pi,mag); %plotting magnitude of H
        xlabel('Frequency (Hz)');
        ylabel('|H(f)|');
        title('Magnitude');
        subplot(3,2,4);
        plot(f*180/pi,phase); %plotting phase of H
        xlabel('Frequency (Hz)');
        ylabel('<H(f) (degrees)');
        title('Phase');
        subplot(3,2,5);
        zplane(num,den); %plotting zeros and poles of the system
        title('Zeros and Poles');
        subplot(3,2,6);
        plot(t,y2); %plotting the signal after filteration in time domain
        title('Gain');
        xlabel('time (s)');
        ylabel('Amplitude');
        yt = yt + (y2*bandgains(i));
    end
end

Y = fftshift(fft(y));
Ymag = abs(Y); %magnitude of the original signal
Yphase = angle(Y); %phase of the original signal
Yt = fftshift(fft(yt));
Ytmag = abs(Yt); %magnitude of the filtered signal
Ytphase = angle(Yt); %phase of the filtered signal
Fvec = linspace(-Fs/2,Fs/2,N);
figure;subplot(2,2,1);
plot(Fvec,20*log10(Ymag),Fvec,20*log10(Ytmag)); %plotting magnitudes of both signals
xlabel('Frequency (Hz)');
ylabel('Magnitude response (dB)');
legend('Original Signal','Filtered Signal');
subplot(2,2,2);
plot(Fvec,Yphase*180/pi,Fvec,Ytphase*180/pi); %plotting phases of both signals
ylabel('Phase response');
xlabel('Frequency (Hz)');
legend('Original Signal','Filtered Signal');
subplot(2,2,[3 4]);
plot(t,y,t,yt); %plotting both signals in time domain
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal','Filtered Signal');

sound(yt,new_Fs);