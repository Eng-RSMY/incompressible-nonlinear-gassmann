Grid.Nx=60;  Grid.hx=20*.3;          % Dimension in x-direction
Grid.Ny=60; Grid.hy=10*.3;          % Dimension in y-direction
Grid.Nz=80;  Grid.hz=2*.3;           % Dimension in z-direction
N=Grid.Nx*Grid.Ny*Grid.Nz;              % Number of grid celles
Grid.V=Grid.hx*Grid.hy*Grid.hz;         % Volume of each cells
Fluid.vw=3e-4; Fluid.vo=3e-3;           % Viscosities
Fluid.swc=0.0; Fluid.sor=0.0;           % Irreducible saturations
St = 5;                                 % Maximum saturation time step
Pt = 100;                               % Pressure time step
ND = 4000;                              % Number of days in simulation
Q=zeros(Grid.Nx,Grid.Ny,Grid.Nz);               % Source term for injection
IR=20; %795*(Grid.Nx*Grid.Ny*Grid.Nz/(60*220*85));   % and production. Total
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
%save_bin(0,'perm1',reshape(Grid.K(1,:,:,:),Grid.Nx*Grid.Ny*Grid.Nz,1));
%save_bin(0,'perm2',reshape(Grid.K(2,:,:,:),Grid.Nx*Grid.Ny*Grid.Nz,1));
%save_bin(0,'perm3',reshape(Grid.K(3,:,:,:),Grid.Nx*Grid.Ny*Grid.Nz,1));

load spe_phi.dat;
Phi=reshape(spe_phi',60,220,85);
Por=Phi(1:Grid.Nx,1:Grid.Ny,1:Grid.Nz);
Grid.por=max(Por(:),1e-2);
%save_bin(0,'por',Grid.por);

fname1=sprintf('%s/spe_606080_2/perm1_1_orig',getenv('SCRATCH'));
fname2=sprintf('%s/spe_606080_2/perm2_1_orig',getenv('SCRATCH'));
fname3=sprintf('%s/spe_606080_2/perm3_1_orig',getenv('SCRATCH'));
fname4=sprintf('%s/spe_606080_2/por_1_orig',getenv('SCRATCH'));
fid1=fopen(fname1,'w')
fid2=fopen(fname2,'w')
fid3=fopen(fname3,'w')
fid4=fopen(fname4,'w')
for i = 1:Grid.Ny
    fwrite(fid1, Grid.K(1,:,i,:), 'single');
    fwrite(fid2, Grid.K(2,:,i,:), 'single');
    fwrite(fid3, Grid.K(3,:,i,:), 'single');
    fwrite(fid4, Por(:,i,:), 'single');
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);

return;

I=ones(N,1);
S=Fluid.swc*I;                          % Initial saturation
Pc=[0; 1]; Tt=0;                        % For production curves


rhoGr = 2.670;        % Density of grains, g/cm^3
rhoW  = 1.100;        % Density of water, g/cm^3
rhoHc = 0.600;        % Density of hydrocarbon, g/cm^3

KGr = 38.0;           % Bulk modulus of grains, GPa
KW  = 3.5;            % Bulk modulus of water, GPa
KHc = 0.8;            % Bulk modulus of hydrocarbon, GPa

KC = 15.0;            % bulk modulus of wet clay, GPa


C = 0;               % Clay content
Pe = 0.4;            % Effective pressure in kBar
%Pe = 9.8 * 2500 * 1000 / 1e+5 / 1e+3; % Effective pressure in kBar

Vp = 5.77*I - 6.94*Grid.por - 1.73*sqrt(C)*I + 0.446*(Pe - exp(-16.7*Pe))*I; % Vp according to Eberhart-Phillips et al, km/s
Vs = 3.70*I - 4.94*Grid.por - 1.57*sqrt(C)*I + 0.361*(Pe - exp(-16.7*Pe))*I; % Vs according to Eberhart-Phillips et al, km/s
save_bin(0,'vp',Vp);
save_bin(0,'vs',Vs);

bulk_modulus = @(rho, vp, vs) rho .* (vp.^2 - 4.0/3.0*vs.^2);
shear_modulus = @(rho, vs) rho .* vs.^2;

rhoB  = (I - Grid.por).*rhoGr + Grid.por.*rhoW;
K = bulk_modulus(rhoB,Vp,Vs);
G = shear_modulus(rhoB,Vs);
save_bin(0,'Kinit',K);
save_bin(0,'Ginit',G);

K0 = 0.5 * (C*KC + (1-C)*KGr + 1.0/(C/KC + (1-C)/KGr));

Kfl = I ./ (S./KW + (I - S)./KHc);

Kstar=(K.*(Grid.por*K0./KW + I - Grid.por) - K0*I) ./ (Grid.por*K0./KW + K/K0 - I - Grid.por);
save_bin(0,'Kstar',Kstar);

Ip_initial = rhoB .* Vp; % initial impedance
        
rhoFl = S*rhoW + (I - S)*rhoHc;
rhoB  = (I - Grid.por).*rhoGr + Grid.por.*rhoFl;
save_bin(0,'rhoB',rhoB);

for tp=1:ND/Pt;
        [P,V]=Pres(Grid,S,Fluid,Q);                             % Pressure solver
        for ts=1:Pt/St;
                S=NewtRaph(Grid,S,Fluid,V,Q,St);                % Implicit saturation solver
                %subplot('position' ,[0.05 .1 .4 .8]);           % Make left subplot
                %Splane=S(1:Grid.Nx*Grid.Ny,1);
                %pcolor(reshape(Splane,Grid.Nx,Grid.Ny)');       % Plot saturation
                %shading flat; caxis([Fluid.swc 1-Fluid.sor]);   %   
                [Mw,Mo]=RelPerm(S(N),Fluid); Mt=Mw+Mo;          % Mobilities in well-block
                Tt=[Tt,(tp-1)*Pt+ts*St];                        % Compute simulation time
                Pc=[Pc,[Mw/Mt; Mo/Mt]];                         % Append production data
                %subplot('position' ,[0.55 .1 .4 .8]);           % Make right subplot
                %plot(Tt,Pc(1,:),Tt,Pc (2,:));                   % Plot production data
                %axis([0,ND,-0.05,1.05]);                        % Set correct axis
                %legend('Water cut','Oil cut');                  % Set legend
                %drawnow;                                        % Force update of plot
        end 

        Kfl = I ./ (S./KW + (I - S)./KHc);
        rhoFl = S*rhoW + (I - S)*rhoHc;

        Ksat = Kstar + (I-Kstar./K0).^2 ./ (Grid.por./Kfl + (I-Grid.por)./K0 - Kstar./K0.^2);

        rhoB  = (I - Grid.por).*rhoGr + Grid.por.*rhoFl;
        vp = sqrt((Ksat + 4.0/3.0*G)./rhoB);
        vs = sqrt(G./rhoB);

        save_bin(tp,'pressure',P);
        save_bin(tp,'saturation',S);
        save_bin(tp,'Ksat',Ksat);
        save_bin(tp,'rhoB',rhoB);
        save_bin(tp,'vp',vp);
        save_bin(tp,'vs',vs);

        Vp_diff = Vp - vp; % Vp is initial velocity computed according to Eberhart-Phillips et al.
        Ip = vp .* rhoB; % Impedance

        rho_above = 2.2;
        vp_above = 4.2;
        I_above = vp_above * rho_above * I;
        RefAmp = (I_above - Ip) ./ (I_above + Ip);
        RefAmp_initial = (I_above - Ip_initial) ./ (I_above + Ip_initial);
        RefAmp_diff = RefAmp_initial - RefAmp;

        save_bin(tp,'VPdiff',Vp_diff);
        save_bin(tp,'Ip',Ip);
        save_bin(tp,'RefAmp',RefAmp);
        save_bin(tp,'RefAmpDiff',RefAmp_diff);
end
