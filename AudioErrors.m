clear all, clc, close all;

[y1, Fs1] = audioread("./WavFiles/16.wav");
[y2, Fs2] = audioread("./WavFiles/45.wav");
[y3, Fs3] = audioread("./WavFiles/26.wav");

s1 = y1(:, 1);
s2 = y2(:, 1);
s3 = y3(:, 1);

crt = 0.0001;
n = 3;

disp(size(y1));
disp(size(y2));
disp(size(y3));

len = size(y1);
t = linspace(0, 5, len(1));

figure(1);
subplot(3, 2, 1);
plot(t, y1(:, 1));
title("Audio Signal 1");
xlabel('---> Time(s)');   ylabel('---> s_1(t)');
grid on;


subplot(3, 2, 3);
plot(t, y2(:, 1));
title("Audio signal 2");
xlabel('---> Time(s)');   ylabel('---> s_2(t)');
grid on;


subplot(3, 2, 5);
plot(t, y3(:, 1));
title("Audio signal 3");
xlabel('---> Time(s)');   ylabel('---> s_3(t)');
grid on;



s1 = transpose(lowpass(s1, 10000, Fs1));
s2 = transpose(lowpass(s2, 15000, Fs2));
s3 = transpose(lowpass(s3, 15000, Fs3));


fshift1 = (-length(s1)/2:length(s1)/2-1)*(Fs1/length(s1));
yshift1 = fftshift(fft(s1))/length(s1);
subplot(3, 2, 2);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 1');
xlabel('---> Frequency(Hz)');   ylabel('---> S_1(f)');
% xlim([-100, 100]);      yticks(0:0.25:1);
grid on;

fshift2 = (-length(s2)/2:length(s2)/2-1)*(Fs2/length(s2));
yshift2 = fftshift(fft(s2))/length(s2);
subplot(3, 2, 4);
plot(fshift2, yshift2);
title('Frequency Spectrum of signal 2');
xlabel('---> Frequency(Hz)');   ylabel('---> S_2(f)');
% xlim([-100, 100]);      yticks(0:0.25:1);
grid on;

fshift3 = (-length(s3)/2:length(s3)/2-1)*(Fs3/length(s3));
yshift3 = fftshift(fft(s3))/length(s3);
subplot(3, 2, 6);
plot(fshift3, yshift3);
title('Frequency Spectrum of signal 3');
xlabel('---> Frequency(Hz)');   ylabel('---> S_3(f)');
% xlim([-100, 100]);      yticks(0:0.25:1);
grid on;


figure(2);

signals = [s1; s2; s3];

[transmission, tTime] = transmitterCommutator(crt, signals, n, t);
transmission = awgn(transmission, 30);

subplot(2, 1, 1);
plot(tTime, transmission);
xlabel('---> t(s)');     ylabel('---> S(t)');
title(sprintf('Time Domain Multiplexed Signal'));
grid on;

audiowrite('TransmissionOutput.wav', transmission, 3/crt);

fshift1 = (-length(transmission)/2:length(transmission)/2-1)*(Fs1/length(transmission));
yshift1 = fftshift(fft(transmission))/length(fftshift(abs(fft(transmission))));
subplot(2, 1, 2);
plot(fshift1, yshift1);
% xlim([-15000, 15000]);
xlabel('---> Frequency(Hz)');     ylabel('---> S(f)');
title(sprintf('Frequency Spectrum of Multiplexed Signal'));
grid on;


[decoded, timestamps] = demultiplexer(crt, transmission, n, tTime);

sizes = size(timestamps);

endInd = sizes(2); 
disp("EndInd = " + endInd);
figure(3);
titles = ["Audio signal 1", "Audio signal 2", "Audio Signal 3"];
k = 1;

for i=1:n
    subplot(3, 2, k);
    k = k + 1;
    plot(timestamps(i, 1:endInd-1), decoded(i, 1:endInd-1));
    xlabel("---> Time(secs)");           ylabel("---> s(t)");
    title(titles(i));
    grid on;
    
    n1 = length(timestamps(i, :));
    fshift1 = (-n1/2:n1/2-1)*(1/crt/n1);
    yshift1 = fftshift(fft(decoded(i, :)));
    yshift1 = yshift1/length(yshift1);

    subplot(3, 2, k);
    k = k + 1;
    plot(fshift1, yshift1);
    grid on;
    title("Frequency spectrum of " + titles(i));
    xlabel('---> Frequency(Hz)');   ylabel('---> S(f)');
end
sgtitle('Demultiplexed Signals');


audiowrite('output1.wav', decoded(1, 1:endInd-1), Fs1/5);
audiowrite('output2.wav', decoded(2, 1:endInd-1), Fs2/5);
audiowrite('output3.wav', decoded(3, 1:endInd-1), Fs3/5);





function [trans, transtime] = transmitterCommutator(crt, signals, noOfSig, time)
    trans = [];                  %Transmission signal
    transtime = [];              %Time for transmission
    j = 1;                       %Current signal index
    delT = crt/noOfSig;          %Switching time
    T = 0;                       %T keeps track of time to switch
    k = 1;
    for i=1:length(time)
        if(time(i) >= T)
            trans(k) = signals(j ,i);
            transtime(k) = T;
            T = T + delT;
            k = k + 1;
            j = j + 1;
            if(j > noOfSig)
                j = 1;
            end
        end
    end
end

function [decoded, timestamps] = demultiplexer(crt, transmission, noOfSig, time)
    
    for j = 1:noOfSig
        i = 1;
        while (j+noOfSig*(i-1)) <= length(transmission)
            decoded(j, i) = transmission(j+(noOfSig*(i-1)));
            timestamps(j, i) = time(j+(noOfSig*(i-1)));
            i = i + 1;
        end
    end
end
