function phi2stl(phi, varargin)
% Developed by Panagiotis Vogiatzis (panagiotis.vogiatzis@stonybrook.edu)
% Advisor: Prof. Shikui Chen (shikui.chen@stonybrook.edu)
% Computational Modeling Analysis and Design Optimization Lab (CMADO)
% Department of Mechanical Engineering
% State University of New York at Stony Brook
% First version: November 06, 2015
%%% Disclaimer: This code is provided for educational purposes and is not
%%% guaranteed to be free of errors. The authors are not liable for any
%%% problems caused by the use of this code.

% phi2stl(phi) takes the level-set function Phi and converts it into an
% stl file using binary mode. The code contains also post-processing tools.
% After the post-processing, Phi is saved as 'phi.mat'
%   -phi: a matrix with the values of the 3D or 4D level-set function phi
%   -2D (3D) Design
%       A 2D (3D) design is described by a 3D (4D) level-set function phi.
%       phi is a i x j (i x j x k) matrix, where i and j (i, j and k)
%       define the grid of the 2D (3D) design.
%   -i x j x k is the initial grid.
%   -If the design has been discretized in n elements in one direction, the
%    corresponding grid will have n+1 nodes.
%   -The original resolution is set to be 1/(max([i, j, k])-1). (res = 1)
%   -The code transforms, initially, the larger side of the design
%    automatically to 1mm. (r = 1)
%   -Visualization is off (vis = 0)
%   -The default thickness for 2D designs is set to be 20%.

% phi2stl(phi, r) includes the scale factor r as an input
%   r will scale the dimension of the larger side to r mm. The dimension of
%   the other 1 or 2 direction/s will be resized accordingly.
%   i.e. r = 10 will make the largest dimension equal to 10mm.

% phi2stl(phi, r, res) includes resolution res. With default r = 1, if
%   the initial grid is 100x100x100, the resolution is set to be 1/100mm.
%   If res = 0.5, the resolution will be half of the initial (1/50mm).
%   This parameter increases/decreases the accuracy of the design. High
%   resolution should be carefully selected since it will significantly
%   increase the computational cost. res = 1 keeps the initial resolution.

% phi2stl(phi, r, res, vis) manages also the visualization. vis = 0/1
%   disables/enables the visualization. Disabling the visualization will
%   accelerate the code.
%
% phi2stl(phi, r, res, vis, t) Sets thickness for 2D designs.
%   Percentage of the thickness compared to the initial x dimension.
%   i.e. t = 50 makes the thickness 50% of the initial x dimension.
[r, res, vis, t] = inputs(varargin{:}); % Investigate which mode is used.

% % % % %PARAMETERS TO BE CHANGED% % % % %
% Flip parameters
%   The design can be flipped in one or more directions. These parameters
%   can be used in order to correctly mirror the design later.
fl_x = 0; % fl_x = 0 will not flip the design in the x direction.
fl_y = 0;
fl_z = 0; % fl_z = 1 will flip the design in the z direction.

% Symmetry
symx = 1; % symx = 1 will not apply symmetry in x direction.
symy = 1; % symy = 2 will apply symmetry in y direction.
symz = 1; % In 2D designs, symz will not have any effect on the design.

% Periodicity
%   Defining the number that the design is repeated in each direction.
%   If a 4x5x1 structure is needed: per_x = 4, per_y = 5 and per_z = 1.
per_x = 1;
per_y = 1;
per_z = 1; % In 2D designs, instead of using per_z, increase t parameter.

% Orientation
%   Orientation and flip parameters will define the final design placement.
%   1, 2 and 3 refers to initial x, y and z directions.
%   Each of 1, 2 and 3 have to be used once.
%   The orientation of the design can be changed with the following:
x_new = 1; % x_new=2 will move the initial y direction to x direction.
y_new = 2; % Default: x_new = 1, y_new = 2 and z_new = 3.
z_new = 3;


% % % % %EXAMPLES% % % % %
% % Michell beam (online download is required)
% load('michell_beam.mat'); r = 110; vis = 1; t = 10;

% % 3D Microgripper (online download is required)
% load('microgripper.mat'); r=50; vis = 1; fl_z = 1; symz = 2;

% % 2D Example
% [g{1}, g{2}] = ndgrid(0:0.02:1, 0:0.02:1); vis = 1;
% phi = sqrt((g{1}-0.5).^2 + (g{2}-0.5).^2) - 0.25; per_x = 3; per_y = 2;

% Dimensions of initial phi
[i, j, k] = size(phi);
mx = max([i, j, k]);
d(1) = (i-1)/(mx-1);
di=d(1)*r;
d(2) = (j-1)/(mx-1);
d(3) = (k-1)/(mx-1);
dx = max([d(1), d(2), d(3)])/(mx-1);
[g{1}, g{2}, g{3}] = ndgrid(0:dx:d(1), 0:dx:d(2), 0:dx:d(3));
% since the largest dimension is scaled to 1, set max(phi) equal to 0.2.
phi = 0.2*phi/max(phi(:));
% % if vis==1
% %     figure;
% %     if ndims(phi)==2
% %         surf(g{1}, g{2}, phi); title('Initial Phi'); axis equal;
% %         figure;
% %         contourf(g{1}, g{2}, phi, [0 0]); axis equal;
% %     else
% %         show(g{1}, g{2}, g{3}, phi, d(1), d(2), d(3))
% %     end
% %     title('Initial Design')
% %     saveas(gcf,'Initial Design.jpeg');
% % end

% Reshaping Phi if res<1
if res<1
    dx = dx/res;
    [g, phi] = resolution(dx, g, phi, d(1), d(2), d(3));
    [i, j, k] = size(phi);
    mx = max([i, j, k]);
    phi = r*phi;
end

% Flip the design in x, y and z direction
if fl_x==1
    phi = flip(phi,1);
end
if fl_y==1
    phi = flip(phi,2);
end
if fl_z==1
    phi = flip(phi,3);
end

% Symmetry in x, y and z direction
phiy = flip(phi(:,1:j-1,:),2);
phixy = [phi, (symy-1)*phiy;
    (symx-1)*flip(phi(1:i-1,:,:),1),(symx-1)*(symy-1)*flip(phiy(1:i-1,:,:),1)];
phixyz = flip(phixy,3);
phisym(:,:,:) = phixy;
if ndims(phi)==3
    phisym(:,:,k:2*k-1) = (symz-1)*phixyz;
end
phi=phisym(1:(mx-1)*d(1)*symx+1,1:(mx-1)*d(2)*symy+1,1:(mx-1)*d(3)*symz+1);

% Periodicity in x, y and z direction
d(1) = per_x*r*d(1)*symx;
d(2) = per_y*r*d(2)*symy;
d(3) = per_z*r*d(3)*symz;
[m, n, o] = size(phi);
dx = max([d(1), d(2), d(3)])/(max([(symx*per_x*(i-1)), ...
    (symy*per_y*(j-1)), (symz*per_z*(k-1))]));
[g{1}, g{2}, g{3}] = ndgrid(0:dx:d(1), 0:dx:d(2), 0:dx:d(3));
phi_new = zeros(per_x*(m-1)+1, per_y*(n-1)+1, per_z*(o-1)+1);
for i = 0:per_x-1
    for j = 0:per_y-1
        for k = 0:per_z-1
            for ii = 1:m
                for jj = 1:n
                    for kk = 1:o
                        phi_new((i*(m-1)+ii), ((j*(n-1)+jj)), ...
                            ((k*(o-1)+kk)))= phi(ii, jj, kk);
                    end
                end
            end
        end
    end
end
phi = phi_new;

% Reshaping Phi if res>=1
if res>=1
    dx = dx/res;
    [g, phi] = resolution(dx, g, phi, d(1), d(2), d(3));
    phi = r*phi;
end

if vis==1
    if ndims(phi)==2
        figure;
        surf(g{1}, g{2}, phi); title('Modified Phi'); axis equal;
        figure;
        title('Modified Design')
        contourf(g{1}, g{2}, phi, [0 0]); axis equal; 
        saveas(gcf,'Modified Design.jpeg');
    end
end

% Convert 2D to 3D (2.5d)
if ndims(phi)==2
    [phi, d(3), g] = convert_2d_to_3d(phi, d(1), d(2), t, di);
end

% Adding boundaries
phi_boundary = boundary(g, d(1), d(2), d(3));
phi = min(phi, phi_boundary);

disp('Preparing stl file...');
% Rotate the design for stl output only
gf{1} = g{x_new};
gf{2} = g{y_new};
gf{3} = g{z_new};
df(1)=d(x_new);
df(2)=d(y_new);
df(3)=d(z_new);
% Limit. Here, positive phi defines material, so +phi and +0.00001.
phi_stl = phi; % if positive/negative defines material (+/-phi).
limit = 0.00001; % the sign of the limit should match the sign of phi.
phi_stl(phi_stl<=limit) = -1;
phi_stl(phi_stl>limit) = 1;

% Smoothing the design
phi_stl = smooth3(phi_stl);

if vis==1
    figure;
    title('Design')
    show(gf{1}, gf{2}, gf{3}, phi_stl, df(1), df(2), df(3))
    saveas(gcf,'Design.jpeg');
end

% Getting the faces and vertices
[faces, vertices] = isosurface(gf{1}, gf{2}, gf{3}, phi_stl, 0);
triangles = triangulation(faces, vertices);

% % if vis==1
% %     figure
% %     title('Mesh')
% %     h=trimesh(triangles);
% %     set(h,'EdgeColor', [0 0 0.5]);
% %     axis equal;
% %     grid on
% %     axis([0 df(1)  0  df(2)  0 df(3)])
% %     xlabel('x');  ylabel('y');  zlabel('z');
% %     saveas(gcf,'Mesh.jpeg');
% % end

normals = faceNormal(triangles);
vertices=vertices'; faces=faces'; normals=normals';
% Transform vertices into single precision.
datav = single(vertices);
% Put all vertices of the faces in a single 2D matrix (3 x m).
datavxyz = datav(:,faces);
% Transform data to [V1, V2, V3] for each face in a (3 x 3 x m/3) matrix.
dataxyz = reshape(datavxyz,3,3,numel(datavxyz)/9);
% Reshape normals from (3 x m) to (3 x 1 x m) matrix.
normals = reshape(normals,3,1,numel(normals)/3);
% Include normal vectors: [n, V1, V2, V3] in a (3 x 4 x m/3) matrix.
datanxyz = [normals(:,:,:), dataxyz(:,:,:)];
% Transform datanxyz to 16-bit unsigned integer.
datanxyz = typecast(datanxyz(:), 'uint16');
% Reshape data. Each column will contain the information of one face.
data = reshape(datanxyz(:), 24, numel(datanxyz)/24);
% The 25th is the attribute byte count.
data(25,:) = 0;

% Generating stl file in Binary mode
fileID = fopen('.\stl.stl', 'w'); % name of stl file can be changed.
% Write a title up to 80 characters. (80 bytes)
fprintf(fileID, '%-80s',...
    'exported using phi2stl created by P. Vogiatzis (Advisor: S. Chen)');
% Write the number of faces in 32-bit unsigned integer.
fwrite(fileID, numel(data)/25, 'uint32');
% Write data.
fwrite(fileID, data, 'uint16');
fclose(fileID);
disp('...stl file ready');

% Saving new phi data
assignin('base', 'phi', phi);
clearvars -except phi
save('phi.mat')

function [r, res, vis, t] = inputs(varargin)
r = 1; % Default size 1mm, unless otherwise specified
res = 1; % Default resolution is the initial resolution
vis = 0; % By default, visualization is off
t = 20; % Default thickness for 2D designs is set to 20%
if nargin==1 % phi2stl(phi, r)
    r = varargin{1};
elseif nargin==2 % phi2stl(phi, r, res)
    r = varargin{1};
    res = varargin{2};
elseif nargin==3 % phi2stl(phi, r, res, vis)
    r = varargin{1};
    res = varargin{2};
    vis = varargin{3};
elseif nargin==4 % phi2stl(phi, r, res, vis)
    r = varargin{1};
    res = varargin{2};
    vis = varargin{3};
    t = varargin{4};
elseif nargin>0
    disp('too many inputs!')
end

function [phi, ti, g] = convert_2d_to_3d(phi, dimx, dimy, ti, di)
[m, n] = size(phi);
ti = di*ti/100;
dx = max([dimx, dimy])/((max([m, n])-1));
[gb{1}, gb{2}] = ndgrid(0:dx:dimx, 0:dx:dimy);
[g{1}, g{2}, g{3}] = ndgrid(0:dx:dimx, 0:dx:dimy, 0:dx:ti);

phi = interpn(gb{1}, gb{2}, phi, g{1}(:), g{2}(:));
phi = reshape(phi, size(g{1}));
phi1 = g{1};
phi2 =  max(g{1}(:))-g{1};
phi3 = g{2};
phi4 =  max(g{2}(:))-g{2};
phi5 = g{3};
phi6 = max(g{3}(:))-g{3};

phi12 = min(phi1, phi2);
phi123 = min(phi12, phi3);
phi1234 = min(phi123, phi4);
phi12345 = min(phi1234, phi5);
phi_boundary = min(phi12345, phi6);
phi = min(phi, phi_boundary);

function [g, phi] = resolution(dx_new, g2, phi, dimx, dimy, dimz)
[g{1}, g{2}, g{3}] = ndgrid(0:dx_new:dimx, 0:dx_new:dimy, 0:dx_new:dimz);
if ndims(phi)==2
    phic = interpn(g2{1}, g2{2}, phi, g{1}(:), g{2}(:));
else
    phic = interpn(g2{1}, g2{2}, g2{3}, phi, g{1}(:), g{2}(:), g{3}(:));
end
phi = reshape(phic, size(g{1}));

function [phi] = boundary(g, dimx, dimy, dimz)
phi1 = g{1};
phi2 = dimx-g{1};
phi3 = min(phi1, phi2);
phi4 = g{2};
phi5 = dimy-g{2};
phi6 = min(phi4, phi5);
phi7 = g{3};
phi8 = dimz-g{3};
phi9 = min(phi7, phi8);
phi10 = min(phi3, phi6);
phi = min(phi10, phi9);

function show(x, y, z, Phi, dimx, dimy, dimz)
h = patch(isosurface(x, y, z, Phi, 0));
set(h, 'FaceColor', [0 1 0], 'EdgeColor', 'none');
camlight left;
grid on
axis equal
axis([0 dimx  0  dimy  0 dimz])
xlabel('x');  ylabel('y');  zlabel('z');
alpha (h, 1);
axis off
% view([30,30])