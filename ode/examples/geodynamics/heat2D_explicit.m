%heat2D_explicit.m
%
% Solves the 2D heat equation with an explicit FD scheme

clear

% Physical parameters
L   = 150e3   % Width of lithosphere [m]
H   = 100e3   % Height of lithosphere [m]
Tbot = 1573.15;  % [K]  %1300;   % Temperature of bottom lithosphere [°C]
Tsurf = 273.15;   %[K] %0;     % Temperature of country rock [°C]
Tplume = 1773.15;  % [K] %1500; % Temperature of plume [°C]
kappa = 1e-6;   % Thermal diffusity of rock [m2/s]
Wplume = 25e3; % Width of plume
day = 3600*24; % # seconds per day
year = 365.25*day; % # seconds per year

% Numerical parameters
nx = 101;   % # gridpoints in x-direction
nz = 51;    % # gridpoints in z-direction
nt = 500;   % Number of timesteps to compute
dx = L/(nx-1) % spacing of grid in x-direction
dz = H/(nz-1) % spacing of grind in z-direction
[x2d,z2d] = meshgrid(-L/2:dx:L/2, -H:dz:0); % create grid

% compute stable timestep
dt = min([dx,dz])^2/kappa/4

% Setup initial linear temperature profile
T = abs(z2d./H)*Tbot;   % vector ./ componentwise division


% Imping plume beneath lithosphere
ind = find(abs(x2d(1,:)) <= Wplume/2);
T(1,ind) = Tplume;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START OF DISCR. SCHEME %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = 0;
for n=1:nt
    % Compute new temperature
    Tnew = zeros(nz,nx);
    % short cuts -- I prefer kappa*dt( .../dx^2 + .../dz^2)
    %sx = kappa*dt/dx^2;
    %sz = kappa*dt/dz^2;
     
    % Set BOUNDARY CONDITIONS
    % const temperature conditions on top and bottom i.e.
    Tnew(1,:) = T(1 ,: );
    %% TODO
    Tnew(nz,:) = Tsurf;
    
    for i=2:nz-1
        % employ zero flux on right and left-hand side of domain i.e.
        % e.g. on LEFT side of domain: T(i,2)-T(i,0)/2dx = 0 ==> T(i,0) =
        % T(i,2). This will be inserted in the discretisation scheme (4)       
        Tnew(i,1) =  T(i,1) + kappa*dt*((2*T(i,2) - 2*T(i,1))/dx^2 + (T(i+1,1) - 2*T(i,1) + T(i-1,1))/dz^2); %%LEFT
        Tnew(i,nx) = T(i,nx) + kappa*dt*((2*T(i,nx-1) - 2*T(i,nx))/dx^2 + (T(i+1,nx) - 2*T(i,nx) + T(i-1,nx))/dz^2); 
    end
    
    
    for j=2:nx-1
        for i=2:nz-1  
            %% TODO. Discretisation scheme (4):
            Tnew(i,j) = T(i,j) + kappa*dt*( (T(i,j+1) - 2*T(i,j) + ...
                                             T(i,j-1))/dx^2 + (T(i+1,j) - ...
                                                              2*T(i,j) + T(i-1,j))/dz^2);
        end
    end
   
    
    T = Tnew;
    time = time+dt;
    
    scale = 1; % 1.e-03
    
    % Plot solution every 10. timesteps
    if (mod(n,10)==0)
       figure(1), clf
       pcolor(x2d*scale,z2d*scale,Tnew); shading interp, colorbar
       hold on
       contour(x2d*scale,z2d*scale,Tnew,[100:100:1500],'k');
       xlabel('x [km]')
       ylabel('z [km]')
       zlabel('Temperature [°C]')
       title(['Temperature evolution after ',num2str(time/year/1e6),' Myrs'])
       drawnow
    end
end                                     
    