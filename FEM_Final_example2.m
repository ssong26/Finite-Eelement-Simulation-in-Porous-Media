% =====================07-Dec-2017 Final Project============================
% Siyuan_Song Final Project for EN234
% 07-Dec-Example 1
%
% Program Description------------------------------------------------------
%
% 1. This Program is intended to calculate the acoustic field in porous media
% for different kinds of structure.
%
% 2. The following variables will be calculated in the current program
%
% ==========================================================================
% =========================The Main Function================================
%
% This example is to investigate the propagation of sound
%
function FEM_Final()
clear;clc;
%
% ----------------------Part 1 Basic parameters----------------------------
%
% ---Material Parameters
density = 1.23;                 % density
viscosity = 1.95*10^(-5);       % viscosity
kata = 0.024;                   % specific heat capacity
gama = 1.4;                     % specific heat ratio
p0 = 1.013*10^(5);              % ambient air pressure
cp = 1004;                      % specifica heat capacity of air at constant pressure
pr = viscosity*cp/kata;         % Prantle number
c0 = 343;                       % sound speed in air
air_impedance = density*c0;     % air impedance
freq = 10000;
coeff = freq*2*pi/c0;
p_background = 1;
%
% ---Structural Parameters
material_length = 0.1;         % the length of the porous material in x direction
material_height = 0.05;         % the height of the porous material in y direction
material_radius = material_height/4;
%
%----------------------------Part 2 Generating Grid------------------------
%
% ----Basic parameters
number_length = 200;            % number of node along radial direction
number_height = 50;            % number of node along circumferential direction
nnode = number_length * number_height; % number of node in the structure
coord_temp = zeros(nnode,2);        % coordinates of the node
x_temp = linspace(0,material_length,number_length);
y_temp = linspace(0,material_height,number_height);
node_count = 0;
% ---Generate the node and the basis connect
for i=1:number_length
    for j=1:number_height
        if (x_temp(i)-material_length/2)^2+(y_temp(j)-material_height/2)^2>material_radius^2
            node_count = node_count + 1;
            coord_temp(node_count,1)= x_temp(i);
            coord_temp(node_count,2)= y_temp(j);
        end
    end
end
nnode = node_count;
coord = coord_temp(1:node_count,:);
connect_basis=delaunay(coord(:,1),coord(:,2));%generate the connectivity array
connect_basis_number=size(connect_basis,1);
connect = [];
connect_temp_number = 0;
for i = 1:connect_basis_number
    point=connect_basis(i,:);%record the three point of the element
    centroids_x=1/3*(coord(point(1),1)+coord(point(2),1)+coord(point(3),1));%record the centroids of the element
    centroids_y=1/3*(coord(point(1),2)+coord(point(2),2)+coord(point(3),2));
    if ((centroids_x-material_length/2)^2+(centroids_y-material_height/2)^2)>material_radius^2%the centroids of the element is not inside the circle r=1
        connect_temp_number=connect_temp_number+1;
        connect(connect_temp_number,:)=point;
    end
end
figure;
triplot(connect,coord(:,1),coord(:,2),'r');
connect_number = connect_temp_number;
%
% ---------------------Part 3 Global Stiffness Matrix-----------------------
%
Stif=zeros(nnode,nnode);
for lmn=1:connect_number    % Loop over all the elements
    %
    %   Set up the stiffness for the current element
    %
    a = connect(lmn,1);
    b = connect(lmn,2);
    c = connect(lmn,3);
    coord_element = zeros(3,2);
    coord_element(1,:)=coord(a,:);
    coord_element(2,:)=coord(b,:);
    coord_element(3,:)=coord(c,:);
    n_nodes=3;n_points=3;
    kel = elstif(coord_element,coeff,n_nodes,n_points);
    %
    %   Add the current element stiffness to the global stiffness
    %
    for i = 1 : 3
        for j = 1 : 3
            Stif(connect(lmn,i),connect(lmn,j)) = Stif(connect(lmn,i),connect(lmn,j)) + kel(i,j);
        end
    end
end
figure
spy(Stif)
p = symrcm(Stif); 
figure
spy(Stif(p,p))
%
%================
%
% ==========================================================================
% =========================Boundary Condition 1=============================
%
% this part is to set u_r=1 when r=a
fixnodes=zeros(number_height,2);
for j=1:number_height
    fixnodes(j,1)=j;
    fixnodes(j,2)=p_background;
end
%
% ==================== Assemble Boundary condition 1=======================
%
% Modify the global stiffness and residual to include constraints
%
resid=zeros(nnode,1);
for i=1:number_height
    for j=1:nnode
        Stif(i,j)=0;
    end
    resid(i) = fixnodes(i,2);
    Stif(i,i)=1;
end
% ================== Solve the FEM equations ==============================
%
pressure = Stif\resid;
%
% ========================Compare the result================================
%
% the results of FEM
pressure_matrix = zeros(number_length,number_height);
x_matrix = zeros(number_length,number_height);
y_matrix = zeros(number_length,number_height);
x_temp = linspace(0,material_length,number_length);
y_temp = linspace(0,material_height,number_height);
number_temp = 0;
for i=1:number_length
    for j=1:number_height
        x_matrix(i,j) = x_temp(i);
        y_matrix(i,j) = y_temp(j);
        if (x_temp(i)-material_length/2)^2+(y_temp(j)-material_height/2)^2>material_radius^2
            number_temp = number_temp + 1 ;
            pressure_matrix(i,j) = pressure(number_temp);
        end
        if (x_temp(i)-material_length/2)^2+(y_temp(j)-material_height/2)^2<=material_radius^2
            node_count = node_count + 1;
            pressure_matrix(i,j) = NaN;
        end
    end
end
figure;
pcolor(x_matrix,y_matrix,pressure_matrix);
colorbar;
%dlmwrite('exp2_results.txt', results);
%===================================end====================================
end
%=================Some Functions may be used in the Main function==========
%
% integrationpoints function for 2D, we can obtain the position and weight
function xi = abq_UEL_2D_integrationpoints(n_points, n_nodes)
xi=zeros(n_nodes,3);        % xi=(x,y,w)
if n_nodes==3               % nodes of the elements
    if n_points==3          % number of integration points
        xi(1,1)=0.5;xi(1,2)=0.5;xi(1,3)=1/6; % integration point 1
        xi(2,1)=0.0;xi(2,2)=0.5;xi(2,3)=1/6; % integration point 2
        xi(3,1)=0.5;xi(3,2)=0.0;xi(3,3)=1/6; % integration point 3
    end
end
end
% shapefunction for 2D, we can obtain the N and dNdx.
function f=abq_UEL_2D_shapefunctions(xi,n_points,n_nodes)
f=zeros(n_nodes,3);
if n_nodes==3
    if n_points==3
        f(1,1)=xi(1);f(1,2)=1;f(1,3)=0;
        f(2,1)=xi(2);f(2,2)=0;f(2,3)=1;
        f(3,1)=1-xi(1)-xi(2);f(3,2)=-1;f(3,3)=-1;
    end
end
end
% stiffness function, we can obtain the stiffness matrix==
function kel = elstif(coord,coeff,n_nodes,n_points)
xi = abq_UEL_2D_integrationpoints(n_points, n_nodes);
kel = 0;
for i = 1 : 3 % integration over all points
    f = abq_UEL_2D_shapefunctions(xi(i,1:2),n_points,n_nodes);
    %
    dNdxi = zeros(3,2);dNdxi(:,1) = (f(:,2))';dNdxi(:,2) = (f(:,3))';
    %
    dxdxi = coord'*dNdxi;
    dxidx = inv(dxdxi);
    %
    dNdx = dNdxi*dxidx;
    determinate = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1);
    % ---C matrix
    if n_points ==3
        if n_nodes==3
            B_wave=zeros(3,2);
            B_wave(:,1) = dNdx(:,1);B_wave(:,2) = dNdx(:,2);
        end
    end
    C = B_wave*B_wave'*xi(i,3)*determinate;
    % ---T matrixv
    if n_points == 3
        if n_nodes == 3
            B_line=zeros(3,1);
            B_line(:,1)=f(:,1);
        end
    end
    T = B_line*B_line'*xi(i,3)*determinate;
    kel = kel + 0.5*(C-coeff^2*T);
end
end
% 
%==================================end=====================================