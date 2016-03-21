Grid.Nx=60;  Grid.hx=20*.3048;          % Dimension in x-direction
Grid.Ny=220; Grid.hy=10*.3048;          % Dimension in y-direction
Grid.Nz=1;   Grid.hz=2*.3048;           % Dimension in z-direction
N=Grid.Nx*Grid.Ny*Grid.Nz;              % Number of grid celles
Grid.V=Grid.hx*Grid.hy*Grid.hz;         % Volume of each cells
Fluid.vw=3e-4; Fluid.vo=3e-3;           % Viscosities
Fluid.swc=0.2; Fluid.sor=0.2;           % Irreducible saturations
St = 5;                                 % Maximum saturation time step
Pt = 100;                               % Pressure time step
ND = 200;                              % Number of days in simulation
Q=zeros(Grid.Nx,Grid.Ny,Grid.Nz);               % Source term for injection
IR=795*(Grid.Nx*Grid.Ny*Grid.Nz/(60*220*85));   % and production. Total
Q(1,1,:)=IR; Q(Grid.Nx,Grid.Ny,:)=-IR; Q=Q(:);  % rate scaled to one layer

load spe_perm.dat;
disp(['spe_perm size: ', num2str(size(spe_perm(:)))]);
Perm=reshape(spe_perm',60,220,85,3);
Kx=Perm(1:Grid.Nx,1:Grid.Ny,1:Grid.Nz,1);
Ky=Perm(1:Grid.Nx,1:Grid.Ny,1:Grid.Nz,2);
Kz=Perm(1:Grid.Nx,1:Grid.Ny,1:Grid.Nz,3);

Grid.K=ones(3,Grid.Nx,Grid.Ny,Grid.Nz);
for iz=1:Grid.Nz
	for iy=1:Grid.Ny
		for ix=1:Grid.Nx
			Grid.K(1,ix,iy,iz)=Kx(ix,iy,iz);
			Grid.K(2,ix,iy,iz)=Ky(ix,iy,iz);
			Grid.K(3,ix,iy,iz)=Kz(ix,iy,iz);
		end
	end
end
figure;surf(log10(reshape(Grid.K(1,:,:,:),Grid.Nx,Grid.Ny)));view(2);drawnow;
figure;surf(log10(reshape(Grid.K(2,:,:,:),Grid.Nx,Grid.Ny)));view(2);drawnow;
figure;surf(log10(reshape(Grid.K(3,:,:,:),Grid.Nx,Grid.Ny)));view(2);drawnow;

load spe_phi.dat;
Phi=reshape(spe_phi',60,220,85);
Por=Phi(1:Grid.Nx,1:Grid.Ny,1:Grid.Nz);
Grid.por=max(Por(:),1e-3);
figure;
surf(reshape(Grid.por,Grid.Nx,Grid.Ny));view(2);drawnow;


figure;
S=Fluid.swc*ones(N,1);                  % Initial saturation
Pc=[0; 1]; Tt=0;                        % For production curves
for tp=1:ND/Pt;
        [P,V]=Pres(Grid,S,Fluid,Q);     % Pressure solver
        for ts=1:Pt/St;
                S=NewtRaph(Grid,S,Fluid,V,Q,St);                % Implicit saturation solver
                subplot('position' ,[0.05 .1 .4 .8]);           % Make left subplot
                pcolor(reshape(S,Grid.Nx,Grid.Ny,Grid.Nz)');    % Plot saturation
                shading flat; caxis([Fluid.swc 1-Fluid.sor]);   %   
                [Mw,Mo]=RelPerm(S(N),Fluid); Mt=Mw+Mo;          % Mobilities in well-block
                Tt=[Tt,(tp-1)*Pt+ts*St];                        % Compute simulation time
                Pc=[Pc,[Mw/Mt; Mo/Mt]];                         % Append production data
                subplot('position' ,[0.55 .1 .4 .8]);           % Make right subplot
                plot(Tt,Pc(1,:),Tt,Pc (2,:));                   % Plot production data
                axis([0,ND,-0.05,1.05]);                        % Set correct axis
                legend('Water cut','Oil cut');                  % Set legend
                drawnow;                                        % Force update of plot
        end
        fname=sprintf('output/saturation_%i',tp);
        fid=fopen(fname,'w');
        if (fid==-1) disp(['Cannot open file ', fname]); return; end 
        fwrite(fid,S,'single');
        fclose(fid); 
end
