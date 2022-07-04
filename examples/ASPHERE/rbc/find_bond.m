% A function to be used to create a initial configuration for lammps simulations of a red blood cell 
% Written by Hongyan Yuan, 2011

function [bond,bond_length]=find_bond(tri,xyz)

n_tri = size(tri,1);
n_bond = n_tri*3/2;

t=1;
    i=tri(t,1);
    j=tri(t,2);
    k=tri(t,3);
    
    bn=1;
    bond(bn,:)=[i j];
    
    bn=bn+1;
    bond(bn,:)=[j k];
        
    bn=bn+1;
    bond(bn,:)=[k i];

    
for t=2:n_tri
    i=tri(t,1);
    j=tri(t,2);
    k=tri(t,3);

    overlap=0;
    for a=1:bn
       if (i==bond(a,1) && j==bond(a,2)) || (j==bond(a,1) && i==bond(a,2))
           overlap=1;
       end
    end
    
    if overlap==0
      bn=bn+1;
      bond(bn,:)=[i j];        
    end
    
    overlap=0;
    for a=1:bn
       if (j==bond(a,1) && k==bond(a,2)) || (k==bond(a,1) && j==bond(a,2))
           overlap=1;
       end
    end
    
    if overlap==0
      bn=bn+1;
      bond(bn,:)=[j k];        
    end
    
    overlap=0;
    for a=1:bn
       if (i==bond(a,1) && k==bond(a,2)) || (k==bond(a,1) && i==bond(a,2))
           overlap=1;
       end
    end
    
    if overlap==0
      bn=bn+1;
      bond(bn,:)=[i k];        
    end    
end

if n_bond==bn
    
else
    
   output= 'something wrong!' 
end

bond_length=zeros(n_bond,1);
for i=1:n_bond
   r1=xyz(bond(i,1),:) ;
   r2=xyz(bond(i,2),:) ;   
   
   bond_length(i)=norm(r1-r2);   
end
    

