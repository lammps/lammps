DATE: 2017-11-20 UNITS: real CONTRIBUTOR: Michal Kanski michal.kanski@uj.edu.pl CITATION: Michal Kanski, Dawid Maciazek, Zbigniew Postawa, Adri C.T. van Duin, Chowdhury Ashraf, and Barbara J. Garrison Development of a Charge-Implicit ReaxFF Potential for Hydrocarbon Systems, J. Phys. Chem. Lett., 2018, 9, pp 359–363,  DOI: 10.1021/acs.jpclett.7b03155
 39       ! Number of general parameters                                        
   50.0000 !p(boc1)                                                             
    9.5469 !p(boc2)                                                             
   26.5405 !p(coa2)                                                             
    1.5105 !p(trip4)                                                            
    6.6630 !p(trip3)                                                            
   70.0000 !kc2                                                                 
    1.0588 !p(ovun6)                                                            
    4.6000 !p(trip2)                                                            
   12.1176 !p(ovun7)                                                            
   13.3056 !p(ovun8)                                                            
  -70.1292 !p(trip1)                                                            
    0.0000 !Lower Taper-radius (swa)                                            
    7.0000 !Upper Taper-radius (swb)                                            
    0.0000 !not used                                                            
   33.8667 !p(val7)                                                             
    6.0891 !p(lp1)                                                              
    1.0563 !p(val9)                                                             
    2.0384 !p(val10)                                                            
    6.1431 !not used                                                            
    6.9290 !p(pen2)                                                             
    0.3989 !p(pen3)                                                             
    3.9954 !p(pen4)                                                             
    0.0000 !not used                                                            
    5.7796 !p(tor2)                                                             
   10.0000 !p(tor3)                                                             
    1.9487 !p(tor4)                                                             
    0.0000 !not used                                                            
    2.1645 !p(cot2)                                                             
    1.5591 !p(vdW1)                                                             
    0.1000 !Cutoff for bond order*100 (cutoff)                                  
    2.1365 !p(coa4)                                                             
    0.6991 !p(ovun4)                                                            
   50.0000 !p(ovun3)                                                            
    1.8512 !p(val8)                                                             
    0.0000 !not used                                                            
    0.0000 !not used                                                            
    0.0000 !not used                                                            
    0.0000 !not used                                                            
    2.6962 !p(coa3)                                                             
  2    ! Nr of atoms; atomID;ro(sigma); Val;atom mass;Rvdw;Dij;gamma            
            alfa;gamma(w);Val(angle);p(ovun5);n.u.;chiEEM;etaEEM;n.u.           
            ro(pipi);p(lp2);Heat increment;p(boc4);p(boc3);p(boc5),n.u.;n.u.    
            p(ovun2);p(val3);n.u.;Val(boc);p(val5);n.u.;n.u.;n.u.               
 C    1.5807   4.0000  12.0000   2.1339   0.2705   0.9000   1.3561   4.0000     
      7.0578   2.2030   4.0000  33.2433  79.5548   1.0000   7.0000   0.0000     
      1.2857   0.0000 199.0303  25.5946  34.7987  33.8858   0.8563   0.0000     
     -3.2763   4.0000   1.0564   4.0000   1.1829   1.8737   0.9340  10.7610     
 H    0.8431   1.0000   1.0080   1.7993   0.0531   1.0206  -0.1000   1.0000     
      5.1365   2.3597   1.0000   0.0000 121.1250   1.0000   7.0000   1.0000     
     -0.1000   0.0000  62.4879   1.7831  10.0956   1.5343   1.0698   0.0000     
    -39.1433   2.0000   1.0338   1.0000   1.4937   1.4730   0.2395   9.1813     
  3      ! Nr of bonds; at1;at2;De(sigma);De(pi);De(pipi);p(be1);p(b            
                      p(be2);p(bo3);p(bo4);n.u.;p(bo1);p(bo2)                   
  1  1 143.8966  78.5316  96.1991  -0.6337  -1.0884   1.0000  19.6513   0.4644  
         1.0466  -0.6287   8.3354   1.0000  -0.1410   7.4059   1.0000   0.0000  
  1  2 142.2786   0.0000   0.0000  -0.7910   0.0000   1.0000   7.6804   0.6354  
         8.9091   1.0000   0.0000   1.0000  -0.1759   9.1526   0.0000   0.0000  
  2  2 166.9928   0.0000   0.0000  -0.2894   0.0000   1.0000   9.0000   0.7986  
        22.4216   1.0000   0.0000   1.0000  -0.1025   5.0250   0.0000   0.0000  
  1    ! Nr of off-diagonal terms. at1;at2;Dij;RvdW;alfa;ro(sigma);r            
  1  2   0.0384   1.4571   8.0036   1.4872  -1.0000  -1.0000                    
  6   ! Nr of angles. at1;at2;at3;Thetao,o;p(val1);p(val2);p(coa1);             
  1  1  1  66.9930  26.8666   1.3549  -0.9012   4.6851  18.8775   1.8087        
  1  1  2  65.9137   8.8776   9.1600   0.0000   0.3144   0.0000   4.3449        
  2  1  2  70.6653  32.2568   1.7885   0.0000   2.1528   0.0000   9.9760        
  1  2  2   0.0000  40.0000   1.0000   0.0000   3.2849   0.0000   9.0719        
  1  2  1   0.0000  33.4427   5.0679   0.0000   1.1188   0.0000  10.0000        
  2  2  2   0.0000  27.9213   5.8635   0.0000   0.0000   0.0000   1.0400        
  6    ! Nr of torsions. at1;at2;at3;at4;;V1;V2;V3;p(tor1);p(cot1);n            
  1  1  1  1   0.0500  17.1779   0.1145  -3.7630  -1.4900   0.0000   0.0000     
  1  1  1  2  -0.4619  27.8830   0.2611  -3.8205  -2.0485   0.0000   0.0000     
  2  1  1  2   0.1629  22.7282   0.3893  -2.9569  -7.8729   0.0000   0.0000     
  0  1  2  0   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000     
  0  2  2  0   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000   0.0000     
  0  1  1  0   0.0000  50.0000   0.3000  -4.0000  -2.0000   0.0000   0.0000     
  0    ! Nr of hydrogen bonds. at1;at2;at3;r(hb);p(hb1);p(hb2);p(hb3            
