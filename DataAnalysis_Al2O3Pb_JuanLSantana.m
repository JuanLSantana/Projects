clc
close all
%%
close all
I1mbar = table2array(I_V_1mbar(:,2));
V1mbar = table2array(I_V_1mbar(:,1));
Cond1mbar = derivada(I1mbar',V1mbar');
Cond1mbar = Cond1mbar';
figure(1)
grid on, hold on,
plot(V1mbar,I1mbar, '. r', 'Markersize', 4) 
%plot(V1mbar, I1mbar, '- k')
%axis ([-1e-2 1e-2 -4e-6 4e-6])
hold off
numrow = numel(Cond1mbar);
OscCond = ones(numrow/2-3,1);
for k = 1:numrow/2-3
    OscCond(k,1) =  Cond1mbar(numrow/2+3+k,1);
end
Cond = smooth(OscCond);
figure(2)
grid on, hold on
plot(V1mbar,Cond1mbar, '. r', 'Markersize', 4) 
%plot(V1mbar,Cond,'- k')
xlabel('Voltage')
ylabel('\sigma')
Vcut = -V1mbar(numrow/2+4:numrow,1);

%Phonon-Electron coupling demo
figure(3)
grid on, hold on
plot(Vcut,OscCond, '. r', 'Markersize', 4) 
plot(Vcut,Cond,'- k')
xlabel('Voltage')
ylabel('\sigma')

%% Comparación de métodos de derivación


dV = diff(V1mbar);
Cond2 = ones(numrow-1,1);
for k = 1:numrow-1
    Cond2(k,1) = (I1mbar(k+1,1)-I1mbar(k,1))./dV(k);
end

SCond2 = smooth(Cond2);
figure(4)
grid on, hold on
plot(V1mbar(1:numrow-1,1),Cond2, '. r', 'Markersize', 4) 
%plot(V1mbar(1:numrow-1,1),SCond2,'- k')
xlabel('Voltage')
ylabel('\sigma')



%% Corriente BCS teórica

%Conductancia a Temp Ambiente

I = table2array(IVcte(:,2));
V = table2array(IVcte(:,1));

Condaux = I./V;
Cond = sum(Condaux)./numel(Condaux);

Condder = sum(derivada(I',V'))./numel(derivada(I',V'));
n = 5;
Iserie = ones(n,1);
gap = 0.0015;
%q = 1.6e-19;
q =1;
k = 8.62e-5;
T = 3;
syms v

for m = 1:n
    Iserie(m,1) = 2.*Condder.*(-1).^(m+1).*besselk(1,m.*q.*gap./(k.*T)).*sinh(m.*q.*v./(k.*T));
end

%%
close all
syms x
F = [2.*Condder.*gap.*(-1).^(1+1).*besselk(1,1.*gap./(k.*T)).*sinh(1.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(2+1).*besselk(1,2.*gap./(k.*T)).*sinh(2.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(3+1).*besselk(1,3.*gap./(k.*T)).*sinh(3.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(4+1).*besselk(1,4.*gap./(k.*T)).*sinh(4.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(5+1).*besselk(1,5.*gap./(k.*T)).*sinh(5.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(6+1).*besselk(1,6.*gap./(k.*T)).*sinh(6.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(7+1).*besselk(1,7.*gap./(k.*T)).*sinh(7.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(8+1).*besselk(1,8.*gap./(k.*T)).*sinh(8.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(9+1).*besselk(1,9.*gap./(k.*T)).*sinh(9.*q.*v./(k.*T)), ...
    2.*Condder.*gap.*(-1).^(10+1).*besselk(1,10.*gap./(k.*T)).*sinh(10.*q.*v./(k.*T))];
Isum = sum(F);
Iteo = matlabFunction(Isum);

Vteo = 0:0.00001:0.0015;

figure(6)
grid on, hold on,
plot(Vteo,Iteo(Vteo))























