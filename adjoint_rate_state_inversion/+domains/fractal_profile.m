% self-similar fractal fault profile
% Eric M. Dunham
% original: November 13, 2012
% modified: August 31, 2023 by Vidar Stiernström

% fault length L, grid spacing h=L/N
% periodicity is assumed for profile, so y(L/2)=y(-L/2)

% NOTE: fault will have ODD number N+1 grid points located at nodes
% x = [-L/2:h:L/2] (length L, grid spacing h=L/N)
% Fourier method will give N points with replica length L
% so an extra point, at same (x,y) as first point, is appended at end

% write_files = false; % write files in FDMAP format
% endian = 'n'; % endian of binary files
function [x,y] = fractal_profile(xlim, m, alpha, cutoff, shift_profile, x_shift, seed)
default_arg('alpha', 10^(-2)) % amplitude-to-wavelength ratio of roughness
default_arg('cutoff', true)
default_arg('shift_profile', true); % shift profile in y direction so that yend=0
default_arg('x_shift', 0);  % shift profile in x direction by this amount)
default_arg('seed', 4);

L = xlim(2)-xlim(1);
N = m-1;
h = L/N;

if cutoff
  Lmin = 20*h; % minimum wavelength (nominally 20*h)
end

x = [xlim(1):h:xlim(2)]; % coordinate at nodes

% wavenumber  
kN = pi/h; % Nyquist wavenumber
k = zeros(1,N);
k(1:N/2+1) = 2*[0:N/2]/N;
k(N:-1:N/2+2) = -k(2:N/2);
k = k*kN;

% white noise, unit normal distribution
s = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(s);
y = randn([1,N]); % Gaussian white noise

% scale so PSD has unit amplitude
y = y*sqrt(N/L);

% FFT
Y = fft(y)*h;

% calculate PSD and check for unit amplitude
PSDy = abs(Y).^2/L;
disp(['PSD = ' num2str(mean(PSDy))])

% multiply Y by square-root of desired PSD
PSDy_exact = (2*pi)^3*alpha^2*k.^(-3);
Y = Y.*sqrt(PSDy_exact);

% remove k=0 component
Y(1) = 0;

% add short wavelength (high wavenumber) cutoff
if cutoff
  kmax = 2*pi/Lmin;
  I = find(abs(k)>kmax); 
  Y(I) = 0;
end

% inverse FFT, explointing conjugate symmetry of real-valued function
y = ifft(Y,'symmetric')/h;
  
% check alpha
alpha_check = sqrt(h*sum(y.^2)/L)/L;
ratio = alpha/alpha_check;
disp(['input alpha=',num2str(alpha),' calculated alpha=',num2str(alpha_check)])
disp(['ratio=',num2str(ratio)])

% repeat first point, shift if needed
y = [y y(1)];
x = x-x_shift;

yend = y(1); 
if shift_profile 
  y = y-yend;
end

return
