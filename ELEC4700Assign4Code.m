% Rudi Hidvary
close all 
clear all
clc

G1 = 1;
G2 = 0.5;
G3 = 0.1;
G4 = 10;
G5 = 0.001

C = 0.25;

L = 0.2

alpha = 100;

Vin = -10:1:10;

GMatrix = [
    1 0 0 0 0 0 0;
    0 1 -1 0 0 0 0;
    G1 (-(G1+G2)) 0 0 0 -1 0;
    0 0 G3 0 0 -1 0;
    0 0 G3 0 0 0 -1;
    0 0 0 G4 (-(G4+G5)) 0 0;
    0 0 0 1 0 0 -alpha]

CMatrix = [
    0 0 0 0 0 0 0;
    0 0 0 0 0 -L 0;
    C -C 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;]


% DC sweep from -10 to 10 Volts 
V3 = [];
V0 = [];
n = 1;
w = 0;
for Vin = -10:1:10
    A = ((CMatrix.*(1i.*0) + GMatrix));% CMatrix is multipled by zero as there is no frequency effects
    FVector = [Vin; 0; 0; 0; 0; 0; 0];
    V = A\FVector;
    V3(n) = V(3); % Input voltage 
    V0(n) = V(5); % Output Voltage 
    n = n+1;
end

n = linspace(-10,10,n-1);

% % Plotting the outputs 
% figure(1)
% subplot(2,1,1)
% plot(n,V3,'k')
% xlabel('Voltage IN')
% ylabel('V3')
% title('V3 Output Voltage')
% grid on 
% 
% 
% subplot(2,1,2)
% plot(n,V0,'b')
% xlabel('Voltage IN')
% ylabel('V0')
% title('V0 Output Voltage')
% grid on 


V0 = [];
n = 1;
Vin = 1; % DC bias of 1 Volt
freq = 10000; % maximum frequency
for w = 0:1:freq
    A = ((CMatrix.*(1i.*w) + GMatrix));
    FVector = [Vin; 0; 0; 0; 0; 0; 0];
    V = A\FVector;
    V3(n) = V(3);
    V0(n) = V(5);
    n = n+1;
end


% Plotting the outputs 
figure(2)
semilogx(0:1:freq,V0/Vin,'k')
xlabel('Frequency')
ylabel('V0')
title('V0 as Function of Frequency')
grid on 

figure(3)
semilogx(0:1:freq,angle(V0/Vin),'k')
xlabel('Frequency')
ylabel('Phase')
title('V0 Phase')
grid on 


figure(4)
semilogx(0:1:freq,20*log10(abs(V0/Vin)),'b')
xlabel('Frequency')
ylabel('V0 in dB')
title('V0 Output Voltage')
grid on

r = normrnd(C,0.05,20000,1);
figure(5)
hist(r,200)
title('C perturbation Values')
xlabel('Value of C')
ylabel('Bin Value')
grid on



V0 = [];
n = 1;
Vin = 1; % DC bias of 1 Volt
freq = pi; % maximum frequency
CVal = C;
for w = 1:1:20000
    
    C = r(w);
    
    CMatrix = [
    0 0 0 0 0 0 0;
    0 0 0 0 0 -L 0;
    C -C 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;];

    A = ((CMatrix.*(1i.*freq) + GMatrix));
    FVector = [Vin; 0; 0; 0; 0; 0; 0];
    V = A\FVector;
    V3(n) = V(3);
    V0(n) = V(5);
    n = n+1;
end

figure(6)
hist(abs(V0/Vin),200)
title('Gain perturbation Values')
xlabel('Value Gain')
ylabel('Bin Value')
grid on

% Transient Simulation 
V1 = [];
V3 = [];
V0 = [];
time = 1000;
Vin = [zeros(1,30) ones(1,970)]; % DC bias of 1 Volt
freq = 0; % maximum frequency
for t = 1:time
    A = ((CMatrix.*(1i.*freq) + GMatrix));
    FVector = [Vin(t); 0; 0; 0; 0; 0; 0]
    V = A\FVector
    V1(t) = V(1);
    V3(t) = V(3);
    V0(t) = V(5);
end

t = 1:time;

figure(7)
subplot(3,1,1)
plot(t,V1)
title('Input Signal')
xlabel('time')
ylabel('Input Voltage')
grid on

subplot(3,1,2)
plot(t,V3)
title('V3')
xlabel('time')
ylabel('V3 Voltage')
grid on

subplot(3,1,3)
plot(t,V0)
title('V0')
xlabel('time')
ylabel('V0 Voltage')
grid on

figure(8)
subplot(2,1,1)
plot(abs(fft(V1)))
title('V1 fft')
xlabel('frequency')
ylabel('magnitude')
grid on

subplot(2,1,2)
plot(fftshift(V1))
title('V1 fftshift')
xlabel('frequency')
ylabel('magnitude')
grid on

figure(9)
subplot(2,1,1)
plot(abs(fft(V0)))
title('V0 fft')
xlabel('frequency')
ylabel('magnitude')
grid on


subplot(2,1,2)
plot(fftshift(V0))
title('V0 fftshift')
xlabel('frequency')
ylabel('magnitude')
grid on

% Sinusoidal input

V1 = [];
V3 = [];
V0 = [];
time = 1000;
f = 1/0.03;
Vin = sin(2*pi*f*(t/1000)); % DC bias of 1 Volt
freq = 0; % maximum frequency
for t = 1:time
    A = ((CMatrix.*(1i.*freq) + GMatrix));
    FVector = [Vin(t); 0; 0; 0; 0; 0; 0];
    V = A\FVector;
    V1(t) = V(1);
    V3(t) = V(3);
    V0(t) = V(5);
end

t = 1:time;

figure(10)
subplot(3,1,1)
plot(t,V1)
title('Input Signal')
xlabel('time')
ylabel('Input Voltage')
grid on

subplot(3,1,2)
plot(t,V3)
title('V3')
xlabel('time')
ylabel('V3 Voltage')
grid on

subplot(3,1,3)
plot(t,V0)
title('V0')
xlabel('time')
ylabel('V0 Voltage')
grid on

figure(11)
subplot(2,1,1)
plot(abs(fft(V1)))
title('V1 fft')
xlabel('frequency')
ylabel('magnitude')
grid on

subplot(2,1,2)
plot(fftshift(V1))
title('V1 fftshift')
xlabel('frequency')
ylabel('magnitude')
grid on

figure(12)
subplot(2,1,1)
plot(abs(fft(V0)))
title('V0 fft')
xlabel('frequency')
ylabel('magnitude')
grid on


subplot(2,1,2)
plot(fftshift(V0))
title('V0 fftshift')
xlabel('frequency')
ylabel('magnitude')
grid on


% Gaussian Pulse

V1 = [];
V3 = [];
V0 = [];
time = 1000;
stdDev = 30;
delay = 60;
x = 1:1000;
A = sqrt(1/(stdDev*sqrt(pi)));
B = exp((-((x-delay).^2))./(2*(stdDev^2)));
Vin = (A*7.2939).*B; % Gaussian Pulse
% 
freq = 0; % maximum frequency
for t = 1:time
    A = ((CMatrix.*(1i.*freq) + GMatrix));
    FVector = [Vin(t); 0; 0; 0; 0; 0; 0];
    V = A\FVector;
    V1(t) = V(1);
    V3(t) = V(3);
    V0(t) = V(5);
end

t = 1:time;

figure(13)
subplot(3,1,1)
plot(t,V1)
title('Input Signal')
xlabel('time')
ylabel('Input Voltage')
grid on

subplot(3,1,2)
plot(t,V3)
title('V3')
xlabel('time')
ylabel('V3 Voltage')
grid on

subplot(3,1,3)
plot(t,V0)
title('V0')
xlabel('time')
ylabel('V0 Voltage')
grid on

figure(14)
subplot(2,1,1)
plot(abs(fft(V1)))
title('V1 fft')
xlabel('frequency')
ylabel('magnitude')
grid on

subplot(2,1,2)
plot(fftshift(V1))
title('V1 fftshift')
xlabel('frequency')
ylabel('magnitude')
grid on

figure(15)
subplot(2,1,1)
plot(abs(fft(V0)))
title('V0 fft')
xlabel('frequency')
ylabel('magnitude')
grid on


subplot(2,1,2)
plot(fftshift(V0))
title('V0 fftshift')
xlabel('frequency')
ylabel('magnitude')
grid on

