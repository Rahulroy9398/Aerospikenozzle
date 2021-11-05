AR = 3.81;
gamma=1.23;
eta_b = 0.8;
% input tube inner diameter (the units of this parameter determine the units of the dimensional results)
t_diam = 1.164;
A_t = pi*((t_diam/2)^2); % throat area
A_e = AR*(A_t);
r_e = sqrt(A_e/pi);
M_e = AR2Mach(AR,gamma);
nu_e = Mach2Prandtl(M_e,gamma);
N = 10000000;
M_vals = linspace(1,M_e,N);
AR_vals = Mach2AR(M_vals,gamma);
nu_vals = Mach2Prandtl(M_vals,gamma);
mu_vals = Mach2Mangle(M_vals);
alpha_vals = nu_e-nu_vals+mu_vals;
% non-dimensional values
l_nondim_vals = (1-sqrt(1-(AR_vals.*(1-(eta_b.^2)).*M_vals.*(sin(alpha_vals)./AR))))./sin(alpha_vals);
r_nondim_vals = 1-(l_nondim_vals.*sin(alpha_vals));
x_nondim_vals = l_nondim_vals.*cos(alpha_vals);
y_nondim_vals = l_nondim_vals.*sin(alpha_vals);
Length_nondim = max(x_nondim_vals)-min(x_nondim_vals);
% dimensional values
l_vals = l_nondim_vals.*r_e;
r_vals = r_nondim_vals.*r_e;
x_vals = x_nondim_vals.*r_e;
y_vals = y_nondim_vals.*r_e;
Length = Length_nondim.*r_e;
%Plotting
figure
plot(r_nondim_vals,x_nondim_vals);
xlabel('r/r_e')
ylabel('x/r_e')
figure
plot(x_nondim_vals,y_nondim_vals,0,0,'o');
xlabel('x/r_e')
ylabel('y/r_e')
fprintf('Exit Mach number = %g \n',M_e)
fprintf('Length = %g [in] \n',Length)
fprintf('Cowl Seperation = %g [in] \n \n',min(l_vals))
fprintf( 'Flow Turn Angle = %g [deg] \n',nu_e*180/pi)
% Create a text file containing coordinates for input in CAD
n = 500;
m = N/n;
p = length(x_vals);
x = x_vals(1:m:p);
y = y_vals(1:m:p);
z = zeros(1,n);
A = [x;y;z];
fileID = fopen('Angelino_Aerospike.txt','w');
fprintf(fileID,'%6.10f %12.10f %12.10f\n',A);
fclose(fileID);
function M = Prandtl2Mach(nu,gamma)
M0u = 6; % first initial guess (upper bound)
M0l = 1.25; %second initial guess (lower bound)
maxits = 50; % maximum number of iterations
tolerance = 0.000001; % tolerance
i = 0; % initialize i
f1 = sqrt((gamma+1)/(gamma-1)); % constant function of gamma
f2 = sqrt((gamma-1)/(gamma+1)); % constant function of gamma
f = f1*atan(f2*sqrt((M0u^2)-1))-atan(sqrt((M0u^2)-1))-nu; % Prandtl-Meyer function: f(M) - nu = 0
j = abs(f); % initialize j which checks how close to zero our current guess getsus
while j>tolerance && i<=maxits
f = f1*atan(f2*sqrt((M0l^2)-1))-atan(sqrt((M0l^2)-1))-nu;
k = f1*atan(f2*sqrt((M0u^2)-1))-atan(sqrt((M0u^2)-1))-nu;
y = (M0l*k-M0u*f)/(k-f);
M0l = M0u;
M0u = y;
f = f1*atan(f2*sqrt((M0u^2)-1))-atan(sqrt((M0u^2)-1))-nu;
j = abs(f);
i = i+1;
end
M = M0u;
end
%Converting Mach number to Prandtl-Meyer angle
function [nu] = Mach2Prandtl(M,gamma)
f1 = (gamma+1)/(gamma-1);
f2 = 1/f1;
nu = sqrt(f1).*atan(sqrt((f2.*((M.^2)-1))))-atan(sqrt((M.^2)-1));
end
% Getting the Mach angle from Mach number
function [mu] = Mach2Mangle(M)
mu = asin(1./M);
end
% Mach number to area ratio (AR)
function [AR] = Mach2AR(M,gamma)
f1 = 2/(gamma+1);
f2 = (gamma-1)/2;
f3 = (gamma+1)/(gamma-1);
AR = sqrt((1./(M.^2)).*((f1.*(1+(f2.*(M.^2)))).^f3));
end
% Area ratio to Mach number
function [M] = AR2Mach(AR,gamma)
f1 = 2/(gamma+1);
f2 = (gamma-1)/2;
f3 = (gamma+1)/(gamma-1);
M0 = 4;
M = M0;
Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^f3)-(AR^2);
Tolerance = 0.0001;
maxits = 200;
J = abs(Function);
i = 1;
while J>Tolerance && i<=maxits
Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^(f3))-(AR^2);
dFunction = (f3/(M^2))*((f1*(1+(f2*(M^2))))^(f3-1))*2*f2*f1*M + (-2/(M^3))*((f1*(1+f2*(M^2)))^(f3));
y = M-(Function/dFunction);
M = y;
Function = (1/(M^2))*((f1*(1+(f2*(M^2))))^(f3))-(AR^2);
J = abs(Function);
i = i+1;
end
end