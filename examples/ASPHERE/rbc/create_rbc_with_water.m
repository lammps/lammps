% Create a initial configuration for lammps simulations of a red blood cell immersed in water
% Written by Hongyan Yuan & Ralph Kfoury, 2015

function create_rbc_with_water


d0_water=2.7;                  % distance between water particles in the initial configuration 

d0_vesicle  = 0.97;            % particle distance

D_vesicle=50;                  % distance between the bilayer membrane particles

R_vesicle = D_vesicle/2;       % Diameter of the vesicle in initial configuration

gap = 2;                       % gap between bilayer and spectrin network in initial configuration

D_network = D_vesicle - gap*2; % Diameter of Spectrin Network in initial configuration 

b = 10;                        % extra space between box boundary and vesicle 

n_bead = 11;                   % Number of beads between nodes of spectrin network, *needs to be an odd number

d0_network  = 9;               % Distance between the nodes in the initial configuration or A-A distance 

shape = [1 1 1];               % Shape of bilayer particles, * No need to change this



%% Creates the Initial Configuration of RBC membrane with spectrin network
[natoms_vesicle,xyz_vesicle]=points_on_sphere(D_vesicle,d0_vesicle);

[natoms_anchor,xyz_network]=points_on_sphere(D_network,d0_network);

tri = convhulln(xyz_network);
[bond,~]=find_bond(tri,xyz_network);

n_tri = size(tri,1);
n_bond = n_tri*3/2;

new_n_bond = n_bond*(n_bead+1)+natoms_anchor+n_bond;
new_bond=zeros(new_n_bond,2);

for i=1:n_bond
    atom_i = bond(i,1);
    atom_j = bond(i,2);
    r_i = xyz_network(atom_i,:);
    r_j = xyz_network(atom_j,:);
    
    j=1;
    new_atom_id = natoms_anchor+j+(i-1)*n_bead;
    new_bond_id = j+(i-1)*(n_bead+1);
    new_bond(new_bond_id,1:2) = [atom_i,new_atom_id];
    xyz_network(new_atom_id,1:3) = r_i+j*(r_j-r_i)/(n_bead+1);
    
    for j=2:n_bead
        new_atom_id = natoms_anchor+j+(i-1)*n_bead;
        new_bond_id = j+(i-1)*(n_bead+1);
        new_bond(new_bond_id,1:2) = [new_atom_id-1,new_atom_id];
        xyz_network(new_atom_id,1:3) = r_i+j*(r_j-r_i)/(n_bead+1);
    end
    
    
        new_bond_id = new_bond_id+1;
        new_bond(new_bond_id,1:2) = [new_atom_id,atom_j];                        
end

vesicle_atom_type=ones(1,natoms_vesicle);

for i=1:natoms_anchor
   r=xyz_vesicle - ones(natoms_vesicle,1)*xyz_network(i,:);
   r=sum(r.*r,2);
   [~,ind]=sort(r);

   new_bond(n_bond*(n_bead+1)+i,:) = [ind(1),i+natoms_vesicle];

   vesicle_atom_type(ind(1)) = 2;
end


for i=natoms_anchor+1:natoms_anchor+n_bond
    k=(i-natoms_anchor-1)*n_bead+(n_bead+1)/2;
   r=xyz_vesicle - ones(natoms_vesicle,1)*xyz_network(natoms_anchor+k,:);
   r=sum(r.*r,2);
   [~,ind]=sort(r);

   new_bond(n_bond*(n_bead+1)+i,:) = [ind(1),natoms_anchor+k+natoms_vesicle];
    
   vesicle_atom_type(ind(1)) = 2;
end

natoms_network = natoms_anchor+n_bond*n_bead;

n_point = natoms_vesicle+natoms_network;

density = 1/(4/3*pi*shape(1)*shape(2)*shape(3));


%% creates a water box surrounding the RBC membrane and inside the RBC


x = -R_vesicle-b+d0_water : d0_water : R_vesicle+b-d0_water; 
nx=size(x,2);

inner_r2 = (R_vesicle-d0_water-gap*2)*(R_vesicle-d0_water-gap*2);
outer_r2 = (R_vesicle+d0_water)*(R_vesicle+d0_water);
water_in=0;
water_out=0;
for i=1:nx
   for j=1:nx
       for k=1:nx
           r=[x(i),x(j),x(k)];
           r2=r*r';
           if r2<inner_r2 
               water_in=water_in+1;
              xyz_water_in(water_in,1:3) = r;   
              
           elseif r2>outer_r2
               water_out=water_out+1;
              xyz_water_out(water_out,1:3) = r;               
           end
       end
   end
end
natoms_water=water_in+water_out;
n_point_with_water = n_point+natoms_water;




%% Creates the Read Data File neccessary for a LAMMPS Input File
fid=fopen(['read_data.rbc_D' num2str(D_vesicle) 'AA'  num2str(d0_network) 'm' num2str(n_bead+1) '_N' num2str(n_point_with_water) '_W_water_d0_water'], 'w');

fprintf(fid, 'vesicle with diameter = %12f', D_vesicle);
fprintf(fid, '\n');
fprintf(fid, '%8i atoms\n',n_point_with_water);
fprintf(fid, '%6i bonds\n',new_n_bond);
fprintf(fid, '%8i ellipsoids\n',natoms_vesicle);
fprintf(fid, '\n');
fprintf(fid, '%8i atom types\n',7);
fprintf(fid, '%8i bond types\n', 3);
fprintf(fid, '\n');
fprintf(fid, '%12f %12f xlo xhi\n',-R_vesicle-b,R_vesicle+b);
fprintf(fid, '%12f %12f ylo yhi\n',-R_vesicle-b,R_vesicle+b);
fprintf(fid, '%12f %12f zlo zhi\n',-R_vesicle-b,R_vesicle+b);
fprintf(fid, '\n');


fprintf(fid, '\n');
fprintf(fid, 'Atoms\n\n');


% First writes X,Y,Z coordinates of System to Read Data 
for j=1:natoms_vesicle

       type = vesicle_atom_type(j);
    ellipsoidflag=1;
    molecule_ID=0;
    fprintf(fid, '%8i %6i  %12f %12f %12f  %6i %12f %6i \n',j,type,xyz_vesicle(j,:),ellipsoidflag, density,molecule_ID);
end

for j=1:natoms_anchor
    type = 2+1;
    ellipsoidflag=0;
    molecule_ID=1;
    fprintf(fid, '%8i %6i  %12f %12f %12f  %6i %12f %6i \n',j+natoms_vesicle,type,xyz_network(j,:),ellipsoidflag, density,molecule_ID);
end

for j=natoms_anchor+1:natoms_network
    
    type = 3+1;
    if mod(j-natoms_anchor,n_bead)==(n_bead+1)/2
        type =4+1;
    end
    
    ellipsoidflag=0;
    molecule_ID=2;
    fprintf(fid, '%8i %6i  %12f %12f %12f  %6i %12f %6i \n',j+natoms_vesicle,type,xyz_network(j,:),ellipsoidflag, density,molecule_ID);
end



for j=1:water_in
    type = 5+1;
    ellipsoidflag=0;
    molecule_ID=3;
    fprintf(fid, '%8i %6i  %12f %12f %12f  %6i %12f %6i \n',j+n_point,type,xyz_water_in(j,:),ellipsoidflag, density,molecule_ID);
end

for j=1:water_out
    type = 6+1;
    ellipsoidflag=0;
    molecule_ID=4;
    fprintf(fid, '%8i %6i  %12f %12f %12f  %6i %12f %6i \n',j+n_point+water_in,type,xyz_water_out(j,:),ellipsoidflag, density,molecule_ID);
end


fprintf(fid, '\n');

fprintf(fid, 'Ellipsoids\n\n');

nx=[1 0 0];
for j=1:natoms_vesicle
    n=xyz_vesicle(j,:)/norm(xyz_vesicle(j,:));
    nv = cross(nx,n); 
    nv = nv/norm(nv) ;
    theta2=acos(sum(nx.*n))/2;
    
    fprintf(fid, '%8i %12f %12f %12f  %12f %12f %12f %12f\n',j,shape,cos(theta2),nv(1)*sin(theta2),nv(2)*sin(theta2),nv(3)*sin(theta2));
end
fprintf(fid, '\n\n');



% Outputs the Bond Information between particles in the system to the Read Data File Above
fprintf(fid, 'Bonds\n\n');
for i=1:n_bond*(n_bead+1)
    fprintf(fid, '%8i %8i %8i %8i\n',i,1,new_bond(i,1)+natoms_vesicle,new_bond(i,2)+natoms_vesicle); % Note "+natoms_vesicle"
end


for i=n_bond*(n_bead+1)+1:n_bond*(n_bead+1)+natoms_anchor
    fprintf(fid, '%8i %8i %8i %8i\n',i,2,new_bond(i,1),new_bond(i,2)); 
end

for i=n_bond*(n_bead+1)+natoms_anchor+1:n_bond*(n_bead+1)+natoms_anchor+n_bond
    fprintf(fid, '%8i %8i %8i %8i\n',i,3,new_bond(i,1),new_bond(i,2)); 
end

fclose(fid);



%% Generates a (PSF) Protein Structure File for the System above
fid=fopen(['rbc_D' num2str(D_vesicle) 'AA'  num2str(d0_network) 'm' num2str(n_bead+1) '_N' num2str(n_point_with_water)  '_W_water.psf'], 'w');

fprintf(fid, 'psf for red blood cell with diameter of %12f', D_vesicle);
fprintf(fid, '\n\n');

number_of_remarks = 1;
fprintf(fid, '%8i !NTITLE\n',number_of_remarks);
fprintf(fid, 'REMARKS put something here if necessary', D_vesicle);

% % fprintf(fid, '%6i bonds\n',new_n_bond);
% % fprintf(fid, '%8i ellipsoids\n',natoms_vesicle);
% % fprintf(fid, '\n');
% % fprintf(fid, '%8i atom types\n',5);
% % fprintf(fid, '%8i bond types\n', 3);
% % fprintf(fid, '\n');
% % fprintf(fid, '%12f %12f xlo xhi\n',-R_vesicle-b,R_vesicle+b);
% % fprintf(fid, '%12f %12f ylo yhi\n',-R_vesicle-b,R_vesicle+b);
% % fprintf(fid, '%12f %12f zlo zhi\n',-R_vesicle-b,R_vesicle+b);
fprintf(fid, '\n');


fprintf(fid, '\n');

fprintf(fid, '%8i !NATOM\n',n_point_with_water);

for j=1:natoms_vesicle     
    type = vesicle_atom_type(j);
    fprintf(fid,'%8i %1i %4i %6s %2s %4i %13f %15f %9i \n',j,1,1,'GLU','C',type,0,0,0);
end

for j=1:natoms_anchor
    type = 2+1;
    fprintf(fid,'%8i %1i %4i %6s %2s %4i %13f %15f %9i \n',j+natoms_vesicle,1,1,'ALA','C',type,0,0,0);
end

for j=natoms_anchor+1:natoms_network
    
    type = 3+1;
    if mod(j-natoms_anchor,n_bead)==(n_bead+1)/2
        type =4+1;
    end
        
    fprintf(fid, '%8i %1i %4i %6s %2s %4i %13f %15f %9i \n',j+natoms_vesicle,1,1,'MET','C',type,0,0,0);
end


for j=1:water_in
    type = 6;
    fprintf(fid, '%8i %1i %4i %6s %2s %4i %13f %15f %9i \n',j+n_point,1,1,'MET','C',type,0,0,0);
end

for j=1:water_out
    type = 7;
    fprintf(fid, '%8i %1i %4i %6s %2s %4i %13f %15f %9i \n',j+water_in+n_point,1,1,'MET','C',type,0,0,0);
end


fprintf(fid, '\n');


fprintf(fid, '%6i !NBOND: bonds\n',n_bond*(n_bead+1)+natoms_anchor+n_bond);

new_bond(1:n_bond*(n_bead+1),:) = new_bond(1:n_bond*(n_bead+1),:) + natoms_vesicle;

bond_one_column = reshape(new_bond',[],1);

for i=1 : (n_bond*(n_bead+1)+natoms_anchor+n_bond)*2

    fprintf(fid, '%8i %8i %8i %8i %8i %8i %8i %8i\n',bond_one_column);
end


fclose(fid);
