% This script is limited for educational use. As an assignment the code is
% partially complete. 

% Do NOT share this code without the author's permission
%
% Assignment: % As part of this assignment you need to fix the terms of the
% A0 matrix. All non-zero elements in A0 are correct. Only some of the zero
% elements are correct. Study the script and fix matrix A0.

% Purpose: Linear Stability Analysis of 3D Prandtl Slope Flow
% Author: Inanc Senocak
% Date: 10/21/2019
% Please report bugs to senocak@pitt.edu
%
% This script is a modification of the matlab script of Jerome
% Hoepffner that computes the temporal stability of 2D mixing layers
% http://www.lmm.jussieu.fr/~hoepffner/mixingstab.pdf 
%
% The script makes use of the Matlab differentiation matrix suite of
% Weideman and Reddy. The theory is explained in the following work: 
% J. A. C. Weidemann and S. C. Reddy, A MATLAB Differentiation Matrix
% Suite, ACM Transactions on Mathematical Software, 26, (2000): 465-519

format compact; clear all; clc; close all;

flow.type = -1;       % -1 for katabatic, 1 for anabatic
flow.PI_s = 18.5;     % stratification perturbation number
flow.slope =70;       % slope angle in degrees
flow.prandtl = 0.71;  % Prandtl number

%%%% parameters of the numerical method
%%%%
params.n=150;          % number of collocation nodes
params.h=100.0;        % domain height in z direction

params.kx=[0.24 ];  
params.ky=[0.0];  

params.nkx=1;        
params.nky=1;        

%%%% dimensionless parameters
sign = flow.type;        % -1 for katabatic, 1 for anabatic
PI_s = flow.PI_s;        % stratification perturbation number
slope = flow.slope;      % slope angle
alpha = pi* slope /180;  % slope in radian
Pr = flow.prandtl;       % Prandtl number

N=params.n;          % number of collocation nodes
NKX=params.nkx;      % number of wavenumbers in x
NKY=params.nky;      % number of wavenumbers in y
H=params.h;          % domain height in z direction
kx_min=min(params.kx); kx_max=max(params.kx);
ky_min=min(params.ky); ky_max=max(params.ky);
kx=linspace(kx_min,kx_max,NKX);     % wavenumber space in x direction
ky=linspace(ky_min,ky_max,NKY);     % wavenumber space in y direction
[X, Y] = meshgrid(kx,ky);

[z,DM] = chebdif(N,2); 

aa=0.999999;
xx=H*(1+z)./(1-aa*z);
x=xx;
dzdx =  (1/H)*(1+aa)./((aa*x/H+1).^2);
d2zdx2 = (1/H^2)*aa*(1+aa)./((aa*x/H+1).^3);

M = zeros(N,N);
M2 = zeros(N,N);
M = diag(dzdx);
M2 = diag(d2zdx2);
MM = M*M;

D1=M*DM(:,:,1); 
D2=M2*D1 + MM*DM(:,:,2);

Coe1 = Pr * sin(alpha) / PI_s;
Coe2 = sin (alpha) / PI_s;

%%%% Prandtl base flow
u  = sign * (sqrt(2.0)*sin(x/sqrt(2.0)).*exp(-x/sqrt(2.0)));                                % along slope flow velocity
up = sign * (cos(x/sqrt(2.0)).*exp(-x/sqrt(2.0)) - sin(x/sqrt(2.0)).*exp(-x/sqrt(2.0)));    % first derivative of along slope flow velocity
b  = sign * (sqrt(2.0)*cos(x/sqrt(2.0)).*exp(-x/sqrt(2.0)));                                % base profile for buoyancy 
bp = sign * (-sin(x/sqrt(2.0)).*exp(-x/sqrt(2.0)) -  cos(x/sqrt(2.0)).*exp(-x/sqrt(2.0)));  % first derivative of buoyancy profile

%%% coordinate selectors
uu=1:N; vv=uu+N; ww= vv+N; pp=ww+N; bb=pp+N;

%%%% dynamic matrices: 
Z=zeros(N,N); I=eye(N);

E=blkdiag(-i*I,  -i*I,  -i*I,  Z, -i*I); % temporal operator (u, v, w, p, b)

% As part of this assignment you need to fix the terms of the A0 matrix only.
% None zero elements are correct.

A0=[Coe1*D2           Z,                            Z,       Z,               Coe1*I; ... %u   % [u, v, w, p, b]
          Z,          Z,                            Z,       Z,                    Z; ... %v
          Z,          Z,                      Coe1*D2,     -D1,                    Z; ... %w
          Z,          Z,                          -D1,       Z,                    Z; ... %p
    -Coe2*I,          Z,                            Z,       Z,                   Z];     %b  

A1_kx=[-i*diag(u),            Z,            Z,   -i*I,             Z;  ...%u
                Z,   -i*diag(u),            Z,      Z,             Z;  ...%v
                Z,            Z,   -i*diag(u),      Z,             Z;  ...%w
             -i*I,            Z,            Z,      Z,             Z;  ...%p
                Z,            Z,            Z,      Z,   -i*diag(u)];     %b 

A1_ky=[Z,     Z,    Z,     Z,   Z; ... %u
       Z,     Z,    Z,  -i*I,   Z; ... %v
       Z,     Z,    Z,     Z,   Z; ... %w
       Z,  -i*I,    Z,     Z,   Z; ... %p
       Z,     Z,    Z,     Z,   Z];    %b  DONE!


A2=blkdiag(-I*Coe1,  -I*Coe1,  -I*Coe1,  Z,  -I*Coe2); % (u, v, w, p, b) 

%this next step is for imposing the boundary conditions
bot=[1 N]; % boundaries indices
II=eye(5*N); 
DD=blkdiag(D1,D1,D1,D1,D1);
bbot=[bot;bot+N;bot+2*N;bot+4*N];
bbot1=[bot;bot+N;bot+2*N;bot+3*N;bot+4*N];
bbot2=[bot;bot+N;bot+2*N];

%initialize omega
om_i=zeros(NKY,NKX);
om_r=zeros(NKY,NKX);

for ii = 1 : NKX
    k_x = kx(ii);
    for jj = 1 : NKY
        k_y =ky(jj);
        k2 = (k_x^2+k_y^2);
        % temporal stability analysis: om is the eigenvalue
        AA0 = A0 + k_x*A1_kx + k_y*A1_ky + k2*A2;
        
        %%%% impose boundary conditions
        E(bbot,:) = 0;
        AA0(bbot1,:) = DD(bbot1,:); % Neumann boundary conditions 
        AA0(bbot2,:) = II(bbot2,:); % Dirichlet boundary conditions
        %%%% solve for eigenvalues
        % S: diagonal matrix with eigenvalues
        % U: full matrix whose columns are the corresponding right eigenvectors
        
        [U,S]=eig(AA0, E); 
        %toc;

        %%%% remove small and large eigenmodes
        S=diag(S); rem=abs(S)>100|abs(S)<1e-8; S(rem)=[]; U(:,rem)=[];
        %U=U(1:5*N,:); % remove ku and kv part of eigenvectors

        [eimag,is]=sort(imag(S));
        xs=U(:,is);
        es=S(is);
        om_i(jj,ii)=imag(es(end))
        om_r(jj,ii)=real(es(end))
        %om_i(jj,ii)=max(imag(S));
        
    end
end
%remove negative parts for plotting purposes.
om_i(om_i  < 0.0) = NaN;
% figure;
% contourf(X,Y,om_i,20);
% colorbar;
% xlabel('$k_x$'); ylabel('$k_y$')

set(gca,'FontSize',16)

[row col] = find( om_i == max( max(om_i) ) );
row
col
max(max(om_i))

%solve the eigenvalue problem for the max growth state
%display('ky = ')
%k_y = ky(row)
%display('kx = ')
%k_x = kx(col)
%k2 = (k_x^2+k_y^2);
% temporal stability analysis: om is the eigenvalue
%AA0 = A0 + k_x*A1_kx + k_y*A1_ky + k2*A2;

%%%% impose boundary conditions
%E(bbot,:) = 0;
%E(btop,:) = 0;
%AA0(bbot1,:) = DD(bbot1,:); % Neumann boundary conditions 
%AA0(bbot2,:) = II(bbot2,:); % Dirichlet boundary conditions
%AA0(btop1,:) = DD(btop1,:); % Neumann boundary conditions 
%AA0(btop2,:) = II(btop2,:); % Dirichlet boundary conditions
%%%% solve for eigenvalues
% S: diagonal matrix with eigenvalues
% U: full matrix whose columns are the corresponding right eigenvectors

%[U,S]=eig(AA0, E); 
%toc;

%%%% remove small and large eigenmodes
%S=diag(S); rem=abs(S)>100|abs(S)<1e-8; S(rem)=[]; U(:,rem)=[];
%U=U(1:5*N,:); % remove ku and kv part of eigenvectors

%[eimag,is]=sort(imag(S));
%xs=U(:,is);
%es=S(is);
%display('imag(phase speed)')
%imag(es(end))
%display('real(phase speed)')
%real(es(end))

% zmax = 20;

% figure;
% subplot(1,5,1)
% plot( imag(xs(uu,end)), x, 'k--'); hold on
% plot( real(xs(uu,end)), x, 'k');
% xlabel('$u^{\prime}$'); ylabel('$z_n$');
% ylim([0,zmax])
% set(gca,'fontsize',16)
% 
% subplot(1,5,2)
% plot( imag(xs(vv,end)), x, 'k--'); hold on
% plot( real(xs(vv,end)), x, 'k');
% xlabel('$v^{\prime}$'); 
% ylim([0,zmax])
% set(gca,'fontsize',16)
% 
% subplot(1,5,3)
% plot( imag(xs(ww,end)), x, 'k--'); hold on
% plot( real(xs(ww,end)), x, 'k');
% xlabel('$w^{\prime}$'); 
% ylim([0,zmax])
% set(gca,'fontsize',16)
% 
% subplot(1,5,4)
% plot( imag(xs(pp,end)), x, 'k--'); hold on
% plot( real(xs(pp,end)), x, 'k');
% xlabel('$p^{\prime}$'); 
% ylim([0,zmax])
% set(gca,'fontsize',16)
% 
% subplot(1,5,5)
% plot( imag(xs(bb,end)), x, 'k--'); hold on
% plot( real(xs(bb,end)), x, 'k');
% xlabel('$b^{\prime}$'); 
% ylim([0,zmax])
% set(gca,'fontsize',16)
% 
% figure(2);
% plot( imag(U(vv,end)), x, 'k--'); hold on
% plot( real(U(vv,end)), x, 'k');
% 
% figure(3);
% plot( imag(U(ww,end)), x, 'k--'); hold on
% plot( real(U(ww,end)), x, 'k');
% 
% figure(4);
% plot( imag(U(pp,end)), x, 'k--'); hold on
% plot( real(U(pp,end)), x, 'k');
% 
% figure(5);
% plot( imag(U(bb,end)), x, 'k--'); hold on
% plot( real(U(bb,end)), x, 'k');
%%%% plotting the eigenvalues/eigenvectors

% xc = [real(es(1)) real(es(end))];
% yc = [0 0];
% line(xc,yc,'Color','red'); hold on;
% ylim([-4 1]); xlim('auto');
% 
% figure(2);
% y=x;
%  for ind=1:length(es);
%     %%%% plot one eigenmode
%     h=plot(real(S(ind)),imag(S(ind)),'*'); hold on
%  
% % 
% %   %%%%  plotting command for eigenmodes and callback function 
%    tt=['figure(3);aa=' num2str(ind) '; plot(real(U(uu,aa)),y,''b'',imag(U(uu,aa)),y,''b--'',real(U(vv,aa)),y,''r'',imag(U(vv,aa)),y,''r--'',real(U(ww,aa)),y,''k'',imag(U(ww,aa)),y,''k--'',real(U(pp,aa)),y,''c'',imag(U(pp,aa)),y,''c--'',real(U(bb,aa)),y,''g'',imag(U(bb,aa)),y,''g--''); legend(''ureal'',''uimag'',''vreal'',''vimag'',''wreal'',''wimag'',''preal'',''pimag'',''breal'',''bimag'');xlabel(''amplitude'');ylabel(''y'');'];
%    set(h,'buttondownfcn',tt);
%  end
%  title('Click on eigenvalues to see the eigenmodes')
%  %xlimvec=[-3, 3]; ylimvec=[-3,3];
%  xlim('auto'); ylim([-6 6]); grid on; hold off
