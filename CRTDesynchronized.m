clear all, clc, close all;

fs = 10000;
t = 0:1/fs:5;

%PART 1: DEFINING INPUT SIGNALS ==========================================>
%Signal 1: (Single tone periodic)
fc1 = 10;
s1 = 2*cos(2*pi*fc1*t);

%Signal 2: 
fc2 = 10;
s2 = sawtooth(2*pi*fc2*t,0.5);

%Signal 3:
fc3 = 10;
s3 = square(2*pi*fc3*t);

%Signal 4: (Multitone tone periodic)
fc41 = 1;
fc42 = 10;
s4 = 3*sin(2*pi*fc41*t) + 1*sin(2*pi*fc42*t);

%Signal 5: (Rectangular pulse [Infinite Bandwidth])
s5 = rectangularPulse(t);

n = 5;                            %Number of Signals
signals = [s1; s2; s3; s4; s5];   %Two dimensional array of signals
crt = 0.001;                      %Commutator rotation time

%PART 2: PLOTTING SIGNALS AND SPECTRUM ===================================>
figure(1);
%<--------Signal 1--------->
subplot(5,2,1);
plot(t, s1);
xlabel('---> t(s)');     ylabel('---> s_1(t)');
title(sprintf('Signal 1: 2cos(2%sf_{c1}t) for f_{c1} = %.2f', '\pi', fc1));
grid on;

fshift1 = (-length(s1)/2:length(s1)/2-1)*(fs/length(s1));
yshift1 = fftshift(fft(s1))/length(s1);
subplot(5, 2, 2);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 1');
xlabel('---> Frequency(Hz)');   ylabel('---> S_1(f)');
xlim([-100, 100]);      yticks(0:0.25:1);
grid on;

%<--------Signal 2--------->
subplot(5,2,3);
plot(t, s2);
xlabel('---> t(s)');     ylabel('---> s_2(t)');
title(sprintf('Signal 2: Triangular Wave'));
grid on;

fshift1 = (-length(s2)/2:length(s2)/2-1)*(fs/length(s2));
yshift1 = fftshift(abs(fft(s2)))/length(s2);
subplot(5, 2, 4);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 2');
xlabel('---> Frequency(Hz)');   ylabel('---> S_2(f)');
xlim([-100, 100]);              yticks(0:0.2:0.4);
grid on;

%<--------Signal 3--------->
subplot(5,2,5);
plot(t, s3);
xlabel('---> t(s)');     ylabel('---> s_3(t)');
title(sprintf('Signal 3: Square Wave'));
grid on;

fshift1 = (-length(s3)/2:length(s3)/2-1)*(fs/length(s3)); 
yshift1 = fftshift(abs(fft(s3)))/length(s3);
subplot(5, 2, 6);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 3');
xlabel('---> Frequency(Hz)');       ylabel('---> S_3(f)');
xlim([-1000, 1000]);                yticks(0:0.2:0.7);
grid on;

%<--------Signal 4--------->
subplot(5,2,7);
plot(t, s4);
xlabel('---> t(s)');     ylabel('---> s_4(t)');
title(sprintf('Signal 4: 3sin(2%sf_{c1}t) + sin(2%sf_{c2}t) for f_{c1} = %.2f & f_{c2} = %.2f', '\pi', '\pi', fc41, fc42));
grid on;

fshift1 = (-length(s4)/2:length(s4)/2-1)*(fs/length(s4));
yshift1 = fftshift(abs(fft(s4)))/length(s4);
subplot(5, 2, 8);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 4');
xlabel('---> Frequency(Hz)');   ylabel('---> S_4(f)');
xlim([-100, 100]);              yticks(0:0.5:1.5);
grid on;

%<--------Signal 5--------->
subplot(5,2,9);
plot(t, s5);
xlabel('---> t(s)');     ylabel('---> s_5(t)');
title(sprintf('Signal 5: Rectangular Pulse'));
grid on;

fshift1 = (-length(s5)/2:length(s5)/2-1)*(fs/length(s5));
yshift1 = fftshift(abs(fft(s5)))/length(s5);
subplot(5, 2, 10);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 5');
xlabel('---> Frequency(Hz)');   ylabel('---> S_5(f)');
xlim([-200, 200]);                yticks(0:0.1:0.5);
grid on;

sgtitle('Signals To Transmit');

%PART 3: BANDLIMITING THE SIGNALS=========================================>
s1 = lowpass(s1, 100, fs);
s2 = lowpass(s2, 100, fs);
s3 = lowpass(s3, 400, fs);
s4 = lowpass(s4, 100, fs);
s5 = lowpass(s5, 100, fs);

%PART 4: PLOTTING THE BANDLIMITED SIGNALS=================================>
figure(2);
subplot(5,2,1);
plot(t, s1);
xlabel('---> t(s)');     ylabel('---> s_1(t)');
title(sprintf('Signal 1: 2cos(2%sf_{c1}t) for f_{c1} = %.2f', '\pi', fc1));
grid on;

fshift1 = (-length(s1)/2:length(s1)/2-1)*(fs/length(s1));
yshift1 = fftshift(fft(s1))/length(s1);
subplot(5, 2, 2);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 1');
xlabel('---> Frequency(Hz)');   ylabel('---> S_1(f)');
xlim([-100, 100]);      yticks(0:0.25:1);
grid on;

%<--------Signal 2--------->
subplot(5,2,3);
plot(t, s2);
xlabel('---> t(s)');     ylabel('---> s_2(t)');
title(sprintf('Signal 2: Triangular Wave'));
grid on;

fshift1 = (-length(s2)/2:length(s2)/2-1)*(fs/length(s2));
yshift1 = fftshift(abs(fft(s2)))/length(s2);
subplot(5, 2, 4);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 2');
xlabel('---> Frequency(Hz)');   ylabel('---> S_2(f)');
xlim([-100, 100]);              yticks(0:0.2:0.4);
grid on;

%<--------Signal 3--------->
subplot(5,2,5);
plot(t, s3);
xlabel('---> t(s)');     ylabel('---> s_3(t)');
title(sprintf('Signal 3: Square Wave'));
grid on;

fshift1 = (-length(s3)/2:length(s3)/2-1)*(fs/length(s3)); 
yshift1 = fftshift(abs(fft(s3)))/length(s3);
subplot(5, 2, 6);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 3');
xlabel('---> Frequency(Hz)');       ylabel('---> S_3(f)');
xlim([-1000, 1000]);                yticks(0:0.2:0.7);
grid on;

%<--------Signal 4--------->
subplot(5,2,7);
plot(t, s4);
xlabel('---> t(s)');     ylabel('---> s_4(t)');
title(sprintf('Signal 4: 3sin(2%sf_{c1}t) + sin(2%sf_{c2}t) for f_{c1} = %.2f & f_{c2} = %.2f', '\pi', '\pi', fc41, fc42));
grid on;

fshift1 = (-length(s4)/2:length(s4)/2-1)*(fs/length(s4));
yshift1 = fftshift(abs(fft(s4)))/length(s4);
subplot(5, 2, 8);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 4');
xlabel('---> Frequency(Hz)');   ylabel('---> S_4(f)');
xlim([-100, 100]);              yticks(0:0.5:1.5);
grid on;

%<--------Signal 5--------->
subplot(5,2,9);
plot(t, s5);
xlabel('---> t(s)');     ylabel('---> s_5(t)');
title(sprintf('Signal 5: Rectangular Pulse'));
grid on;

fshift1 = (-length(s5)/2:length(s5)/2-1)*(fs/length(s5));
yshift1 = fftshift(abs(fft(s5)))/length(s5);
subplot(5, 2, 10);
plot(fshift1, yshift1);
title('Frequency Spectrum of signal 5');
xlabel('---> Frequency(Hz)');   ylabel('---> S_5(f)');
xlim([-200, 200]);              yticks(0:0.1:0.5);
grid on;

sgtitle('Bandlimited Signals');

%<========SIGNAL TRANSMISSION===============>
[transmission, tTime] = transmitterCommutator(crt, signals, n, t);

figure(3);
subplot(2,1,1);
plot(tTime, transmission);
xlabel('---> t(s)');     ylabel('---> S(t)');
title(sprintf('Time Domain Multiplexed Signal'));
grid on;

fshift1 = (-length(transmission)/2:length(transmission)/2-1)*(fs/length(transmission));
yshift1 = fftshift(abs(fft(transmission)))/length(fftshift(abs(fft(transmission))));
subplot(2,1,2);
plot(fshift1, yshift1);
% xlim([-15000, 15000]);
xlabel('---> Frequency(Hz)');     ylabel('---> S(f)');
title(sprintf('Frequency Spectrum of Multiplexed Signal'));
grid on;

%PART 6: UTILITY PLOT=====================================================>
figure(4);
plot(tTime, transmission,'LineWidth', 0.5);
hold on;
plot(t, s1,'LineWidth',1.5);
hold on;
plot(t, s2,'LineWidth',1.5);
hold on;
plot(t, s3,'LineWidth',1.5);
hold on;
plot(t, s4,'LineWidth',1.5);
hold on;
plot(t, s5,'LineWidth',1.5);
hold off;
xlabel('---> t(s)');     ylabel('---> S(t)');
legend('TDM Signal','Signal 1','Signal 2','Signal 3','Signal 4', 'Signal 5');
title(sprintf('Visualizing Time Multiplexed Signals'));
grid on;

%[decoded, timestamps] = demultiplexer(crt, transmission, n, tTime);
nWindow = 2;
nCols = length(transmission);
if (mod(nCols, nWindow) ~= 0)
    nCols = nCols + (nWindow - mod(nCols, nWindow));
end

decoded = zeros(n, nCols);
timestamps = zeros(n, nCols);
currentP = 1;
currentS = 1;
currentC = 1;
while (currentP <= length(transmission))
    decoded(currentS, currentC:currentC+1) = transmission(currentP:currentP+1);
    timestamps(currentS, currentC:currentC+1) = tTime(currentP:currentP+1);
    currentP = currentP + nWindow;
    currentS = currentS + 1;
    if (currentS > n)
        currentS = 1;
        currentC = currentC + nWindow;
    end
end

[r, c] = size(decoded);
while (all(decoded(:, c) == 0))
    decoded = decoded(:, 1:c-1);
    timestamps = timestamps(:, 1:c-1);
    [r, c] = size(decoded);
end

sizes = size(timestamps);

endInd = sizes(2); 
disp("EndInd = " + endInd);
figure(5);
titles = ["Cosine signal", "Triangular signal", "Square wave", "Sine combination", "Rectangular pulse"];
k = 1;
for i=1:n
    subplot(5, 2, k);
    k = k + 1;
    plot(timestamps(i, 1:endInd-1), decoded(i, 1:endInd-1));
    xlabel("---> Time(secs)");           ylabel("---> s(t)");
    title(titles(i));
    grid on;
    
    n1 = length(timestamps(i, :));
    fshift1 = (-n1/2:n1/2-1)*(1/crt/n1);
    yshift1 = fftshift(abs(fft(decoded(i, :))));
    yshift1 = yshift1/length(yshift1);

    subplot(5, 2, k);
    k = k + 1;
    plot(fshift1, yshift1);
    grid on;
    title("Frequency spectrum of " + titles(i));
    xlabel('---> Frequency(Hz)');   ylabel('---> S(f)');
end
sgtitle('Demultiplexed Signals');

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
