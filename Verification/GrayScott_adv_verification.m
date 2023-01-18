%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gray-Scott model with consideration of advection   
% Reference: 1. https://itp.uni-frankfurt.de/~gros/StudentProjects/
%            Projects_2020/projekt_schulz_kaefer/#header
%       
%            2. https://mrob.com/pub/comp/xmorphia/ogl/index.html
% 
% compare with verification_2d.inp  
%
% Shuolun Wang, 2022 @ND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Uall] = GrayScott_adv_verification()
clear % clear the workspace variables
clearvars % clear all globle variables
clc % clear the command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Problem data
%
% total time
FinalTime = 5000.0;
%
counter = 1;
% time increment
dtime = 1;
%
% time increment for plotting the temperature
dtOut = 1;
%
% length of ``bar''
%
Length = 50;
%
%
%
%
% conductivity(possibly a function of position and temperature)
%
Du = @(x,t) 4.0e-2; % default 
Dv = @(x,t) 1.0e-2; % default 
%
%
% feeding rate of U 
F_const = 0.028;
%
% kill rate of V
k_const = 0.06; 
%
%
% 
%
%  rate constant for ode 
ph = 10.0;
%ph = 1.0;
%
% Threshold 
Tr = 0.1; 
% 
%
% strength constant for electric field coupling 
sigma = 0.05;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the mesh you want
%
% number of elements
%
nElem = 300;
%
% nodes per element
%
nNode = 2;
%
% number of integration points per element
%
nInt = 1;
%
% the total number of nodes in this mesh
%
nNodes = nElem*nNode - (nElem - 1);
%
% generate the mesh
%
coord = zeros(nNodes,1);
for i=2:nNodes
    coord(i) = coord(i-1) + Length/(nNodes-1);
end
%
% here is the element to nodal connectivity in the form
%  node on left  = connect(elem,1)
%  next node moving right = connect(elem,2)
%  ....
%  node on the right = connect(elem,nNode)
%
connect = zeros(nElem,nNode);
for i=1:nNode
    connect(1,i) = i;
end
for i=2:nElem
    for j=1:nNode
        connect(i,j) = connect(i-1,j) + nNode - 1;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intial condition on the Uall
%
 %Uall = zeros(2*nNodes,1);
 Uall = zeros(3*nNodes,1);

% Uall(1) = 0.5;
% Uall(2) = 0.1;


for i = 1:1:nNodes
    Uall(3*i-2) = 1.0;  
    Uall(3*i-1) = 0.0; 
    Uall(3*i) = 0.0; 
end

% initial condition on first five nodes in abaqus 
Uall(1) = 0.5;
Uall(4) = 0.5;
Uall(7) = 0.5;
Uall(10) = 0.5;
Uall(13) = 0.5;

Uall(2) = 0.25;
Uall(5) = 0.25;
Uall(8) = 0.25;
Uall(11) = 0.25;
Uall(14) = 0.25;


% inital condition on phase variable 
c_old = 1e-5*ones(nNodes,1);



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Boundary conditions
%
% flux boundary condition
%
flux = 0.0;      % magnitude AND direction of flux BC
elemFlux = nElem; % element with the flux BC
nodeFlux = nNode; % local node number to apply the flux BC
%
% temperature boundary condition (possible function of time)
%
%Tbc = @(t) 100.0*(1.0 - cos(-t/30));
%Tbc = @(t) 300+50*sin(t/30);
%Tbc = @(t) 500;

Bc_Phi_negative = @(t)    -1.0;
Bc_Phi_positive = @(t)     1.0;



node_Phi_negative = 3*nNodes; % global node number to apply the temperature BC
node_Phi_positive = 3; % global node number to apply the temperature BC



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin the loop over time
%
timeOut = 0.0;
for time=0:dtime:FinalTime
    
    
    % obtain the temperature field from the last converged step
    %
    %Told = T;
    %Dofsold = Dofs;
    Uallold = Uall;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % begin the nonlinear iteration loop
    %
    iter = 0;
    while(iter<=100)
        
        iter = iter + 1;
        
        K = zeros(3*nNodes,3*nNodes);
        R = zeros(3*nNodes,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over elements
        %
        for elem=1:nElem
            
            % get the coordinates for the nodes on this element
            %
            nodeCoords = zeros(nNode,1);
            for i=1:nNode
                nodeCoords(i,1) = coord(connect(elem,i));
            end
            
            
            % get the nodal temperatures for this element
            %
            Te = zeros(nNode,1);
            Ve = zeros(nNode,1);
            Phie = zeros(nNode,1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TeOld = zeros(nNode,1);
            VeOld = zeros(nNode,1);
            PhieOld = zeros(nNode,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            for i=1:nNode
                row = 3*(connect(elem,i)-1)+1;
                Te(i,1) = Uall(row,1);
                Ve(i,1) = Uall(row+1,1);
                Phie(i,1) = Uall(row+2,1);
                TeOld(i,1) = Uallold(row,1);
                VeOld(i,1) = Uallold(row+1,1);
                PhieOld(i,1) = Uallold(row+2,1);
            end
            
            
            % compute the element residual and tangent
            %
%            [ReU,ReV,KUU,KVV,KUV,KVU] = element(nInt,nNode,nodeCoords,...
%                                                Te,TeOld,Ve,VeOld,...
%                                                Du,Dv,F_const,k_const,dtime);


         [ReU,ReV,RePhi,...
          KUU,KVV,KUV,KVU,...
          KUPhi,KVPhi,KPhiU,...
          KPhiV,KPhiPhi] = element(nInt,nNode,nodeCoords,...
                                   Te,TeOld,Ve,VeOld,...
                                   Phie,PhieOld,...
                                   Du,Dv,F_const,k_const,sigma,dtime);                                 
        
                                              
            % check for any flux boundary conditions on the right side
            %
            %         Rbc = zeros(nNode,1);
            %         if(elem==elemFlux)
            %             Rbc(nodeFlux,1) = Area(Length)*flux;
            %         end
            
            
            % assemble this elements tangent into the gloabl
            %
            for i=1:nNode
                for j=1:nNode
                    row = 3*(connect(elem,i)-1)+1;
                    col = 3*(connect(elem,j)-1)+1;
                    
                    K(row,col) = K(row,col) + KUU(i,j);
                    K(row+1,col+1) = K(row+1,col+1) + KVV(i,j);
                    K(row+2,col+2) = K(row+2,col+2) + KPhiPhi(i,j);
                    
                    K(row,col+1) = K(row,col+1) + KUV(i,j);
                    K(row+1,col) = K(row+1,col) + KVU(i,j);
                    
                    K(row,col+2) = K(row,col+2) + KUPhi(i,j);
                    K(row+1,col+2) = K(row+1,col+2) + KVPhi(i,j);
                    K(row+2,col) = K(row+2,col) + KPhiU(i,j);
                    K(row+2,col+1) = K(row+2,col+1) + KPhiV(i,j);
                    
                end
            end
            
            % assemble this elements residual into the global
            %
            for i=1:nNode
                row = 3*(connect(elem,i) - 1) + 1;
                R(row) = R(row) + ReU(i);
                R(row + 1) = R(row + 1) + ReV(i);
                R(row + 2) = R(row + 2) + RePhi(i);                
            end
            
            
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % impose temperature boundary condition
        %
        % enforce the temperature BC
        %
        Uall(node_Phi_negative) = Bc_Phi_negative(time);
        Uall(node_Phi_positive) = Bc_Phi_positive(time);
        
        %
        % modify the tangent and residual
        %
        K(node_Phi_negative,node_Phi_negative) = (1.e6)*(trace(K)/nNodes);
        R(node_Phi_negative,1) = 0.0;
        
        K(node_Phi_positive,node_Phi_positive) = (1.e6)*(trace(K)/nNodes);
        R(node_Phi_positive,1) = 0.0;
        
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for convergence on the residual
        %
        if(norm(R)<1.e-10)
               %  fprintf('Converged on the residual in %i iterations \n', iter);
            break;
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Global system solve for the temperature feild
        %
        % compute the temperature field correction
        %
        Delta = K\R;
        %
        % compute the updated degree of freedom
        %
        Uall = Uall + Delta;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for convergence on the corrections
        %
        if(norm(Delta)<1.e-10)
            % fprintf('Converged on the corrections in %i iterations \n', iter);
            break;
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    % end the nonlinear iteration loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Postprocessing to plot the activator and supressor field
    %
    if(time>=timeOut)
        
        
        [c_new] = postprocess(Length,coord,Uall,time,dtime,nNodes,...
                              c_old,Tr,ph,counter);
        timeOut = timeOut + dtOut;
        counter = counter + 1;
        c_old = c_new;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
% end the loop over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



function [c_new] = postprocess(Length,coord,Uall,time,dtime,nNodes,...
                               c_old,Tr,ph,counter)

% load abaqus data for virification 
% 2d abaqus 
load('NT11.txt')
load('NT12.txt')   
load('NT13.txt')                           
load('SDV1.txt')     
% 3d abaqus  
load('NT11_3d.txt')
load('NT12_3d.txt')   
load('NT13_3d.txt')                           
load('SDV1_3d.txt')   
             
                           
T = zeros(nNodes,1);
V = zeros(nNodes,1);
Phi = zeros(nNodes,1);

global T_mat
global V_mat
global Phi_mat
global c_mat

% extract the nodal solution from the global vector 
%
for i = 1:1:nNodes
    
    T(i) = Uall(3*i-2);
    T_mat(counter,i) = Uall(3*i-2);
    V(i) = Uall(3*i-1);
    V_mat(counter,i) = Uall(3*i-1);
    Phi(i) = Uall(3*i);
    Phi_mat(counter,i) = Uall(3*i);    
    
end


% time integration of phase variable c
%

for i = 1:1:nNodes
    
    % compare u with threshold 
    if  V(i) <= Tr
        au = 0.49; 
    elseif V(i) > Tr
        au = 0.49 - 2.5*(V(i) - Tr);
    end
    
    % call RTsafe solver
    % 
    x1 = 0;
    x2 = 1;
    xacc = 1.0e-5;
    
    arg = [c_old(i),dtime,ph,au];
    maxit = 500;
    gfunc = @Cfunction;
    
    c_new(i) = rtsafe(x1,x2,xacc,arg,maxit,gfunc); 
    c_mat(counter,i) = c_new(i);
    

end


% here I want to make the Abaqus data less dense 
%
pt_2d = 5; 

NT11_less(:,1) = NT11(1:pt_2d:end,1);
NT11_less(:,2) = NT11(1:pt_2d:end,2);
NT12_less(:,1) = NT12(1:pt_2d:end,1);
NT12_less(:,2) = NT12(1:pt_2d:end,2);
NT13_less(:,1) = NT13(1:pt_2d:end,1);
NT13_less(:,2) = NT13(1:pt_2d:end,2);
SDV1_less(:,1) = SDV1(1:pt_2d:end,1);
SDV1_less(:,2) = SDV1(1:pt_2d:end,2);

pt_3d = 9; 

NT11_3d_less(:,1) = NT11_3d(1:pt_3d:end,1);
NT11_3d_less(:,2) = NT11_3d(1:pt_3d:end,2);
NT12_3d_less(:,1) = NT12_3d(1:pt_3d:end,1);
NT12_3d_less(:,2) = NT12_3d(1:pt_3d:end,2);
NT13_3d_less(:,1) = NT13_3d(1:pt_3d:end,1);
NT13_3d_less(:,2) = NT13_3d(1:pt_3d:end,2);
SDV1_3d_less(:,1) = SDV1_3d(1:pt_3d:end,1);
SDV1_3d_less(:,2) = SDV1_3d(1:pt_3d:end,2);


figure(1) 
h1 = plot(coord/Length,T,'r-','LineWidth',2);
hold on
plot(NT11_less(:,1)/Length,NT11_less(:,2),'ro');
hold on
plot(NT11_3d_less(:,1)/Length,NT11_3d_less(:,2),'rv');
hold on

h2 = plot(coord/Length,V,'k-','LineWidth',2);
hold on 
plot(NT12_less(:,1)/Length,NT12_less(:,2),'ko');
hold on 
plot(NT12_3d_less(:,1)/Length,NT12_3d_less(:,2),'kv');
hold on 


h3 = plot(coord/Length,Phi,'m-','LineWidth',2);
hold on 
plot(NT13_less(:,1)/Length,NT13_less(:,2),'mo');
hold on 
plot(NT13_3d_less(:,1)/Length,NT13_3d_less(:,2),'mv');
hold on 

h4 = plot(coord/Length,c_new,'b-','LineWidth',2);
hold on
plot(SDV1_less(:,1)/Length,SDV1_less(:,2),'bo');
hold on 
plot(SDV1_3d_less(:,1)/Length,SDV1_3d_less(:,2),'bv');
hold off;
xlim([0 1]);
%ylim([0 1]);
xlabel('$x/l$','Interpreter','latex');
ylabel('Normalized value','Interpreter','latex');
title(['time=',num2str(time)]);
%h = legend([h1,h2,h3,h4],'activator $u$','suppressor $v$','electric potential $\Phi$','phase field $c$');
%set(h,'Interpreter','latex','Box','off');
%set(h,'location','northout');
set(gca,'FontSize',23)






% figure(2)
% %contourf(T_mat,100,'LineStyle','none')
% contourf(V_mat,100,'LineStyle','none')
% %contourf(c_mat,100,'LineStyle','none')
% colormap bone
% h=colorbar('SouthOutside');
% xlabel('Position')
% ylabel('Time')
% title('activator $u$','interpreter','latex')
% set(gca,'FontSize',23)

drawnow
end


function [ReU,ReV,RePhi,...
          KUU,KVV,KUV,KVU,...
          KUPhi,KVPhi,KPhiU,...
          KPhiV,KPhiPhi] = element(nInt,nNode,nodeCoords,...
                                   Te,TeOld,Ve,VeOld,...
                                   Phie,PhieOld,...
                                   Du,Dv,F_const,k_const,sigma,dtime)                                            
                                 

% obtain gauss points and weights
%
if(nInt==1)
    [xi,w] = GaussInt1Pt();
elseif(nInt==2)
    [xi,w] = GaussInt2Pt();
elseif(nInt==3)
    [xi,w] = GaussInt3Pt();
elseif(nInt==4)
    [xi,w] = GaussInt4Pt();
elseif(nInt==5)
    [xi,w] = GaussInt5Pt();
else
    error('nInt is not programmed');
end


% obtain nodal coordinates
%
if(nNode==2)
    x1 = nodeCoords(1,1);
    x2 = nodeCoords(2,1);
elseif(nNode==3)
    x1 = nodeCoords(1,1);
    x2 = nodeCoords(2,1);
    x3 = nodeCoords(3,1);
elseif(nNode==4)
    x1 = nodeCoords(1,1);
    x2 = nodeCoords(2,1);
    x3 = nodeCoords(3,1);
    x4 = nodeCoords(4,1);
else
    error('nNode is not programmed');
end


% init
%
% Re = zeros(2*nNode,1);
% Ke = zeros(2*nNode,2*nNode);

ReU = zeros(nNode,1);
ReV = zeros(nNode,1);
RePhi = zeros(nNode,1);

KUU = zeros(nNode,nNode);
KVV = zeros(nNode,nNode);
KUV = zeros(nNode,nNode);
KVU = zeros(nNode,nNode);
KUPhi = zeros(nNode,nNode);
KVPhi = zeros(nNode,nNode);
KPhiU = zeros(nNode,nNode);
KPhiV = zeros(nNode,nNode);
KPhiPhi = zeros(nNode,nNode);



%
% loop over integration points
%
for intPt=1:nInt
    
    % compute shape functions and derivatives
    %
    if(nNode==2)
        [N,B,Jac] = shapeLinear(x1,x2,xi(intPt));
    elseif(nNode==3)
        [N,B,Jac] = shapeQuadratic(x1,x2,x3,xi(intPt));
    elseif(nNode==4)
        [N,B,Jac] = shapeCubic(x1,x2,x3,x4,xi(intPt));
    else
        error('nNode is not programmed');
    end
    
    % current location of this integ point
    %
    x = N*nodeCoords;
    
    
    % compute the temperature and dTdX at this integ point based on nodal
    % values of the temperature
    %
    Told = N*TeOld;
    T = N*Te;
    dTdX = B*Te;
    dTdtime = (T - Told)/dtime;
    
    Vold = N*VeOld;
    V = N*Ve;
    dVdX = B*Ve;
    dVdtime = (V - Vold)/dtime;
    
    Phiold = N*PhieOld;
    Phi = N*Phie;
    dPhidX = B*Phie;
    
    
    
    % calculate some lengthy quantities 
    %
    alpha = - T*V^2 + F_const*(1.0 - T); 
    beta  = T*V^2 - (F_const + k_const)*V;
    
    % update the element residual
    %

    ReU = ReU + Jac*w(intPt)*...
        (...
          transpose(N)*dTdtime ...
        + transpose(B)*Du(x,T)*dTdX...
        - transpose(N)*alpha...
        );    
    
    ReV = ReV + Jac*w(intPt)*...
        (...
         transpose(N)*dVdtime ...
        -sigma*transpose(N)*dPhidX*dVdX...
        +transpose(B)*Dv(x,V)*dVdX...
        -transpose(N)*beta...
        );    
    
    RePhi = RePhi + Jac*w(intPt)*...
        (...
         transpose(B)*dPhidX ...
        );         
    
    
    
    % update the element tangent
    %
    
    KUU = KUU + Jac*w(intPt)*...
        (...
        - (1.0/dtime)*transpose(N)*N...
        - transpose(B)*Du(x,T)*B...
        - (V^2 + F_const)*transpose(N)*N...
        );    
    
    
   
    KVV = KVV + Jac*w(intPt)*...
        (...
        - (1.0/dtime)*transpose(N)*N...
        + transpose(N)*sigma*dPhidX*B...
        - Dv(x,V)*transpose(B)*B...
        + (2*T*V - F_const-k_const)*transpose(N)*N...        
        );    
    
    
    KPhiPhi = KPhiPhi + Jac*w(intPt)*...
        (...
        -transpose(B)*B...       
        ); 
    
    
    
    KUV = KUV + Jac*w(intPt)*...
        (...
        - 2*T*V*transpose(N)*N...      
        );
    
    KVU= KVU + Jac*w(intPt)*...
        (...
        V^2*transpose(N)*N...      
        );    
    
    KVPhi = KVPhi + Jac*w(intPt)*...
        (...
        sigma*transpose(N)*B*dVdX...      
        );    
    
    
    
    
    
end % loop over integration points

return;
end

function [xi,w] = GaussInt1Pt()
% Gauss integration locations and weights for 1pt integration
xi = 0.0;
%
w = 2.0;
return;
end

function [xi,w] = GaussInt2Pt()
% Gauss integration locations and weights for 2pt integration
xi(1) = -sqrt(1.0/3.0);
xi(2) =  sqrt(1.0/3.0);
%
w(1) = 1.0;
w(2) = 1.0;
return;
end

function [xi,w] = GaussInt3Pt()
% Gauss integration locations and weights for 3pt integration
xi(1) = -0.7745966692;
xi(2) =  0.0;
xi(3) =  0.7745966692;
%
w(1) = 0.5555555556;
w(2) = 0.8888888889;
w(3) = 0.5555555556;
return;
end

function [xi,w] = GaussInt4Pt()
% Gauss integration locations and weights for 4pt integration
xi(1) = -0.8611363116;
xi(2) = -0.3399810436;
xi(3) =  0.3399810436;
xi(4) =  0.8611363116;
%
w(1) = 0.3478548451;
w(2) = 0.6521451549;
w(3) = 0.6521451549;
w(4) = 0.3478548451;
return;
end

function [xi,w] = GaussInt5Pt()
% Gauss integration locations and weights for 5pt integration
xi(1) = -0.9061798459;
xi(2) = -0.5384693101;
xi(3) =  0.0;
xi(4) =  0.5384693101;
xi(5) =  0.9061798459;
%
w(1) = 0.2369268851;
w(2) = 0.4786286705;
w(3) = 0.5688888889;
w(4) = 0.4786286705;
w(5) = 0.2369268851;
return;
end

function [N,B,Jac] = shapeLinear(x1,x2,xi)
% shape functions and derivatives for a 2-node 1D element

% element length
%
Le = x2 - x1;

% the shape function matrix
%
N = (1.0/2.0)*[1.0-xi 1.0+xi];

% derivatives of shape functions
%
B = (1.0/Le)*[-1.0 1.0];

% the mapping jacobian
%
Jac = Le/2.0;

return;
end

function [N,B,Jac] = shapeQuadratic(x1,x2,x3,xi)
% shape functions and derivatives for a 3-node 1D element

% the shape function matrix
%
N1 = (1/2)*xi*(xi - 1);
N2 = (1 + xi)*(1 - xi);
N3 = (1/2)*xi*(xi + 1);
N = [N1 N2 N3];


% derivatives of shape functions in the xi coordinate
%
dN1dXi = xi - 1/2;
dN2dXi = -2*xi;
dN3dXi = 1/2 + xi;
%
% the mapping jacobian
%
Jac = dN1dXi*x1 + dN2dXi*x2 + dN3dXi*x3;
%
% derivatives of shape functions in the x coordinate
%
B = [dN1dXi dN2dXi dN3dXi]/Jac;

return;
end

function [N,B,Jac] = shapeCubic(x1,x2,x3,x4,xi)
% shape functions and derivatives for a 4-node 1D element

% the shape function matrix
%
N1 = (-9/16)*(xi+1/3)*(xi-1/3)*(xi-1);
N2 = (27/16)*(xi+1)*(xi-1/3)*(xi-1);
N3 = (-27/16)*(xi+1)*(xi+1/3)*(xi-1);
N4 = (9/16)*(xi+1)*(xi+1/3)*(xi-1/3);
N = [N1 N2 N3 N4];

% derivatives of shape functions in the xi coordinate
%
dN1dXi = (9*xi)/8 - (27*xi^2)/16 + 1/16;
dN2dXi = (81*xi^2)/16 - (9*xi)/8 - 27/16;
dN3dXi = 27/16 - (81*xi^2)/16 - (9*xi)/8;
dN4dXi = (27*xi^2)/16 + (9*xi)/8 - 1/16;
%
% the mapping jacobian
%
Jac = dN1dXi*x1 + dN2dXi*x2 + dN3dXi*x3 + dN4dXi*x4;
%
% derivatives of the shape function in the x coordinate
%
B = [dN1dXi dN2dXi dN3dXi dN4dXi]/Jac;

return;
end

function [x] = rtsafe(x1,x2,xacc,arg,maxit,fhandle)


[fl,df] = fhandle(x1,arg);
[fh,df] = fhandle(x2,arg);

% if (fl < 0)
%     xl = x1;
%     xh = x2;
% else
%     xh = x1;
%     xl = x2;
%     temp = fl;
%     fl = fh;
%     fh = temp;
% end

if (fl == 0)
    x = x1;
    return;
elseif (fh == 0)
    x = x2;
    return;
elseif (fl < 0)
    xl = x1;
    xh = x2;
else
    xh = x1;
    xl = x2;
end

x = 0.5 * (x1 + x2);
dxold = abs(x2 - x1);
dx = dxold;

[f,df] = fhandle(x,arg);

for j=1:maxit
    test1 = ((x - xh)*df - f) * ((x - xl)*df - f);
    test2 = abs(2*f) - abs(dxold * df);
    if (test1 > 0 || test2 > 0)
        dxold = dx;
        dx = 0.5 * (xh - xl);
        x = xl + dx;
        if (xl == x)
            return;
        end
    else
        dxold = dx;
        dx = f/df;
        temp = x;
        x = x - dx;
        if (temp == x)
            return;
        end
    end
    if (abs(dx) < xacc)
        return;
    end
    
    [f,df] = fhandle(x,arg);
    
    if (f < 0)
        xl = x;
    else
        xh = x;
    end
    

end
end

function [F,DF] = Cfunction(x,arg)

c_t    = arg(1);
dt     = arg(2);
ph     = arg(3);
au     = arg(4);


F = c_t - x + dt*ph*(au*x^2 - au*x - x^3 + x^2);

     
DF = -1 + dt*ph*(2*au*x - au - 3*x^2 + 2*x); 


end
