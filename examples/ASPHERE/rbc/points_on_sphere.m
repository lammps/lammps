% A function to be used to create a initial configuration for lammps simulations of a red blood cell 
% Written by Hongyan Yuan, 2011

function  [natoms,xyz]=points_on_sphere(D,d0)


rmin_t=d0;
R1=D/2;

n = 1;
x(n)=0.d0;
y(n)=0.d0;
z(n)=R1;
    
n = 2;

x(n)=0.d0;
y(n)=0.d0;
z(n)=-R1;

for i=2:round(pi/2.0/asin(rmin_t/2.0/R1))
  ti=(i-1)*pi/round(pi/2.0/asin(rmin_t/2.0/R1));
  for j=1:round(pi*sin(ti)/asin(rmin_t/2.0/R1))

    n=n+1;
    fj=(j-1)*2.0*pi/round(pi*sin(ti)/asin(rmin_t/2.0/R1));
    x(n)=R1*sin(ti)*cos(fj);
    y(n)=R1*sin(ti)*sin(fj);
    z(n)=R1*cos(ti);
  
  end
end

natoms = size(x,2);

xyz=[x;y;z]';

xyz=xyz+rand(size(xyz,1),3)/100;

% plot3(x,y,z,'.')
% axis tight equal off