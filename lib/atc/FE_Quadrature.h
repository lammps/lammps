#ifndef FE_QUADRATURE_H
#define FE_QUADRATURE_H

// ATC_Transfer headers
#include "MatrixLibrary.h"
#include "ATC_TypeDefs.h"

using namespace std;

namespace ATC {

  /**
   *  @class  FE_Quadrature
   *  @brief  Stores quadrature schemes
   */
  struct FE_Quadrature {
    int numIPs;
    int numFaceIPs;
    DENS_MAT ipCoords;
    vector<DENS_MAT> ipFaceCoords;
    DENS_MAT ipFace2DCoords;
    DENS_VEC ipWeights;
    DENS_VEC ipFaceWeights;

    FE_Quadrature(FeEltGeometry geo, FeIntQuadrature quad)
    {
      if (geo == HEXA) {
        switch (quad) {
          // Degenerate 1 point quadrature
          case GAUSS1: {
            // Set number of IPs and face IPs
            numIPs = 1;
            numFaceIPs = 1;

            // Size matrices and vectors thereof accordingly
            ipCoords.resize(3,numIPs);
            ipFaceCoords.assign(6,DENS_MAT(3,numFaceIPs));
            ipFace2DCoords.resize(2,numFaceIPs);
            ipWeights.reset(numIPs);
            ipFaceWeights.reset(numFaceIPs);

            // "Matrix" of integration point location
            ipCoords(0,0) = 0.0; 
            ipCoords(1,0) = 0.0; 
            ipCoords(2,0) = 0.0;

            // Integration point for each face
            ipFaceCoords[0](0,0) = -1.0; ipFaceCoords[3](0,0) =  0.0; 
            ipFaceCoords[0](1,0) =  0.0; ipFaceCoords[3](1,0) =  1.0; 
            ipFaceCoords[0](2,0) =  0.0; ipFaceCoords[3](2,0) =  0.0;
            //
            ipFaceCoords[1](0,0) =  1.0; ipFaceCoords[4](0,0) =  0.0; 
            ipFaceCoords[1](1,0) =  0.0; ipFaceCoords[4](1,0) =  0.0; 
            ipFaceCoords[1](2,0) =  0.0; ipFaceCoords[4](2,0) = -1.0;
            //
            ipFaceCoords[2](0,0) =  0.0; ipFaceCoords[5](0,0) =  0.0;
            ipFaceCoords[2](1,0) = -1.0; ipFaceCoords[5](1,0) =  0.0; 
            ipFaceCoords[2](2,0) =  0.0; ipFaceCoords[5](2,0) =  1.0;
            
            // 2D integration scheme for the faces
            ipFace2DCoords(0,0) = 0.0; 
            ipFace2DCoords(1,0) = 0.0;
            
            // Integration point weights
            ipWeights = 8.0;

            // Face integration point weights
            ipFaceWeights = 4.0;

            break;
          }
          // 8 point quadratures
          case NODAL:
          case GAUSS2: {
            // Set number of IPs and face IPs
            numIPs = 8;
            numFaceIPs = 4;

            // Size matrices and vectors thereof accordingly
            ipCoords.resize(3,numIPs);
            ipFaceCoords.assign(6,DENS_MAT(3,numFaceIPs));
            ipFace2DCoords.resize(2,numFaceIPs);
            ipWeights.reset(numIPs);
            ipFaceWeights.reset(numFaceIPs);

            // Dictates difference in node locations for nodal/GAUSS2
            double a = 1.0/sqrt(3.0);
            if (quad == NODAL) a = 1.0;
          
            // Matrix of integration point locations & follows local 
            // conn
            ipCoords(0,0) = -a; ipCoords(0,4) = -a; 
            ipCoords(1,0) = -a; ipCoords(1,4) = -a; 
            ipCoords(2,0) = -a; ipCoords(2,4) =  a;
            //
            ipCoords(0,1) =  a; ipCoords(0,5) =  a; 
            ipCoords(1,1) = -a; ipCoords(1,5) = -a; 
            ipCoords(2,1) = -a; ipCoords(2,5) =  a;
            //
            ipCoords(0,2) =  a; ipCoords(0,6) =  a; 
            ipCoords(1,2) =  a; ipCoords(1,6) =  a; 
            ipCoords(2,2) = -a; ipCoords(2,6) =  a;
            //
            ipCoords(0,3) = -a; ipCoords(0,7) = -a; 
            ipCoords(1,3) =  a; ipCoords(1,7) =  a; 
            ipCoords(2,3) = -a; ipCoords(2,7) =  a;

            // Integration points by face
            ipFaceCoords[0](0,0) = -1; ipFaceCoords[3](0,0) = -a; 
            ipFaceCoords[0](1,0) = -a; ipFaceCoords[3](1,0) =  1; 
            ipFaceCoords[0](2,0) = -a; ipFaceCoords[3](2,0) = -a;
            //                                                        
            ipFaceCoords[0](0,1) = -1; ipFaceCoords[3](0,1) =  a; 
            ipFaceCoords[0](1,1) =  a; ipFaceCoords[3](1,1) =  1; 
            ipFaceCoords[0](2,1) = -a; ipFaceCoords[3](2,1) = -a; 
            //                                                      
            ipFaceCoords[0](0,2) = -1; ipFaceCoords[3](0,2) =  a; 
            ipFaceCoords[0](1,2) =  a; ipFaceCoords[3](1,2) =  1; 
            ipFaceCoords[0](2,2) =  a; ipFaceCoords[3](2,2) =  a; 
            //                                                        
            ipFaceCoords[0](0,3) = -1; ipFaceCoords[3](0,3) = -a; 
            ipFaceCoords[0](1,3) = -a; ipFaceCoords[3](1,3) =  1; 
            ipFaceCoords[0](2,3) =  a; ipFaceCoords[3](2,3) =  a; 

            ipFaceCoords[1](0,0) =  1; ipFaceCoords[4](0,0) = -a;
            ipFaceCoords[1](1,0) = -a; ipFaceCoords[4](1,0) = -a;
            ipFaceCoords[1](2,0) = -a; ipFaceCoords[4](2,0) = -1;
            //
            ipFaceCoords[1](0,1) =  1; ipFaceCoords[4](0,1) =  a;
            ipFaceCoords[1](1,1) =  a; ipFaceCoords[4](1,1) = -a;
            ipFaceCoords[1](2,1) = -a; ipFaceCoords[4](2,1) = -1;
            //
            ipFaceCoords[1](0,2) =  1; ipFaceCoords[4](0,2) =  a;
            ipFaceCoords[1](1,2) =  a; ipFaceCoords[4](1,2) =  a;
            ipFaceCoords[1](2,2) =  a; ipFaceCoords[4](2,2) = -1;
            //
            ipFaceCoords[1](0,3) =  1; ipFaceCoords[4](0,3) = -a;
            ipFaceCoords[1](1,3) = -a; ipFaceCoords[4](1,3) =  a;
            ipFaceCoords[1](2,3) =  a; ipFaceCoords[4](2,3) = -1;
            
            ipFaceCoords[2](0,0) = -a; ipFaceCoords[5](0,0) = -a; 
            ipFaceCoords[2](1,0) = -1; ipFaceCoords[5](1,0) = -a; 
            ipFaceCoords[2](2,0) = -a; ipFaceCoords[5](2,0) =  1;
            //                                                       
            ipFaceCoords[2](0,1) =  a; ipFaceCoords[5](0,1) =  a; 
            ipFaceCoords[2](1,1) = -1; ipFaceCoords[5](1,1) = -a; 
            ipFaceCoords[2](2,1) = -a; ipFaceCoords[5](2,1) =  1; 
            //                         
            ipFaceCoords[2](0,2) =  a; ipFaceCoords[5](0,2) =  a; 
            ipFaceCoords[2](1,2) = -1; ipFaceCoords[5](1,2) =  a; 
            ipFaceCoords[2](2,2) =  a; ipFaceCoords[5](2,2) =  1; 
            //                                                       
            ipFaceCoords[2](0,3) = -a; ipFaceCoords[5](0,3) = -a; 
            ipFaceCoords[2](1,3) = -1; ipFaceCoords[5](1,3) =  a; 
            ipFaceCoords[2](2,3) =  a; ipFaceCoords[5](2,3) =  1; 
            
            // Integration points for all faces ignoring the 
            // redundant dim
            ipFace2DCoords(0,0) = -a; ipFace2DCoords(0,2) =  a; 
            ipFace2DCoords(1,0) = -a; ipFace2DCoords(1,2) =  a; 
            // 
            ipFace2DCoords(0,1) =  a; ipFace2DCoords(0,3) = -a; 
            ipFace2DCoords(1,1) = -a; ipFace2DCoords(1,3) =  a; 

            // Integration point weights
            ipWeights = 1.0;

            // Face integration point weights
            ipFaceWeights = 1.0;

            break;
          }
          // 6 point "face" quadrature
          case FACE: {
            printf("using face quad!\n");
            // Set number of IPs and face IPs
            numIPs = 6;
            numFaceIPs = 4;

            // Size matrices and vectors thereof accordingly
            ipCoords.resize(3,numIPs);
            ipFaceCoords.assign(6,DENS_MAT(3,numFaceIPs));
            ipFace2DCoords.resize(2,numFaceIPs);
            ipWeights.reset(numIPs);
            ipFaceWeights.reset(numFaceIPs);

            // Use GAUSS2 for faces for now...
            double a = 1.0/sqrt(3.0);
          
            // Matrix of integration point locations 
            ipCoords(0,0) =  1.0; ipCoords(0,3) =  0.0; 
            ipCoords(1,0) =  0.0; ipCoords(1,3) = -1.0; 
            ipCoords(2,0) =  0.0; ipCoords(2,3) =  0.0;
            //                                        
            ipCoords(0,1) = -1.0; ipCoords(0,4) =  0.0; 
            ipCoords(1,1) =  0.0; ipCoords(1,4) =  0.0; 
            ipCoords(2,1) =  0.0; ipCoords(2,4) =  1.0;
            //                                        
            ipCoords(0,2) =  0.0; ipCoords(0,5) =  0.0; 
            ipCoords(1,2) =  1.0; ipCoords(1,5) =  0.0; 
            ipCoords(2,2) =  0.0; ipCoords(2,5) = -1.0;

            // Integration points by face
            ipFaceCoords[0](0,0) = -1; ipFaceCoords[3](0,0) = -a; 
            ipFaceCoords[0](1,0) = -a; ipFaceCoords[3](1,0) =  1; 
            ipFaceCoords[0](2,0) = -a; ipFaceCoords[3](2,0) = -a;
            //                                                        
            ipFaceCoords[0](0,1) = -1; ipFaceCoords[3](0,1) = -a; 
            ipFaceCoords[0](1,1) =  a; ipFaceCoords[3](1,1) =  1; 
            ipFaceCoords[0](2,1) = -a; ipFaceCoords[3](2,1) =  a; 
            //                                                      
            ipFaceCoords[0](0,2) = -1; ipFaceCoords[3](0,2) =  a; 
            ipFaceCoords[0](1,2) =  a; ipFaceCoords[3](1,2) =  1; 
            ipFaceCoords[0](2,2) =  a; ipFaceCoords[3](2,2) =  a; 
            //                                                        
            ipFaceCoords[0](0,3) = -1; ipFaceCoords[3](0,3) =  a; 
            ipFaceCoords[0](1,3) = -a; ipFaceCoords[3](1,3) =  1; 
            ipFaceCoords[0](2,3) =  a; ipFaceCoords[3](2,3) = -a; 

            ipFaceCoords[1](0,0) =  1; ipFaceCoords[4](0,0) = -a;
            ipFaceCoords[1](1,0) = -a; ipFaceCoords[4](1,0) = -a;
            ipFaceCoords[1](2,0) = -a; ipFaceCoords[4](2,0) = -1;
            //
            ipFaceCoords[1](0,1) =  1; ipFaceCoords[4](0,1) =  a;
            ipFaceCoords[1](1,1) =  a; ipFaceCoords[4](1,1) = -a;
            ipFaceCoords[1](2,1) = -a; ipFaceCoords[4](2,1) = -1;
            //
            ipFaceCoords[1](0,2) =  1; ipFaceCoords[4](0,2) =  a;
            ipFaceCoords[1](1,2) =  a; ipFaceCoords[4](1,2) =  a;
            ipFaceCoords[1](2,2) =  a; ipFaceCoords[4](2,2) = -1;
            //
            ipFaceCoords[1](0,3) =  1; ipFaceCoords[4](0,3) = -a;
            ipFaceCoords[1](1,3) = -a; ipFaceCoords[4](1,3) =  a;
            ipFaceCoords[1](2,3) =  a; ipFaceCoords[4](2,3) = -1;
            
            ipFaceCoords[2](0,0) = -a; ipFaceCoords[5](0,0) = -a; 
            ipFaceCoords[2](1,0) = -1; ipFaceCoords[5](1,0) = -a; 
            ipFaceCoords[2](2,0) = -a; ipFaceCoords[5](2,0) =  1;
            //                                                       
            ipFaceCoords[2](0,1) = -a; ipFaceCoords[5](0,1) =  a; 
            ipFaceCoords[2](1,1) = -1; ipFaceCoords[5](1,1) = -a; 
            ipFaceCoords[2](2,1) =  a; ipFaceCoords[5](2,1) =  1; 
            //                         
            ipFaceCoords[2](0,2) =  a; ipFaceCoords[5](0,2) =  a; 
            ipFaceCoords[2](1,2) = -1; ipFaceCoords[5](1,2) =  a; 
            ipFaceCoords[2](2,2) =  a; ipFaceCoords[5](2,2) =  1; 
            //                                                       
            ipFaceCoords[2](0,3) =  a; ipFaceCoords[5](0,3) = -a; 
            ipFaceCoords[2](1,3) = -1; ipFaceCoords[5](1,3) =  a; 
            ipFaceCoords[2](2,3) = -a; ipFaceCoords[5](2,3) =  1; 
            
            // Integration points for all faces ignoring the 
            // redundant dim
            ipFace2DCoords(0,0) = -a; ipFace2DCoords(0,2) =  a; 
            ipFace2DCoords(1,0) = -a; ipFace2DCoords(1,2) =  a; 
            // 
            ipFace2DCoords(0,1) =  a; ipFace2DCoords(0,3) = -a; 
            ipFace2DCoords(1,1) = -a; ipFace2DCoords(1,3) =  a; 

            // Integration point weights
            ipWeights = 4.0/3.0;

            // Face integration point weights
            ipFaceWeights = 1.0;

            break;
          }
          // 27 point quadrature
          case GAUSS3: {
            // Set number of IPs and face IPs
            numIPs = 27;
            numFaceIPs = 9;

            // Size matrices and vectors thereof accordingly
            ipCoords.resize(3,numIPs);
            ipFaceCoords.assign(6,DENS_MAT(3,numFaceIPs));
            ipFace2DCoords.resize(2,numFaceIPs);
            ipWeights.reset(numIPs);
            ipFaceWeights.reset(numFaceIPs);

            double a = sqrt(3.0/5.0);
          
            // Matrix of integration point locations & follows local 
            // conn
            ipCoords(0,0)  = -a; ipCoords(0,16) = -a; 
            ipCoords(1,0)  = -a; ipCoords(1,16) = -a; 
            ipCoords(2,0)  = -a; ipCoords(2,16) =  0;
            //                                      
            ipCoords(0,1)  =  a; ipCoords(0,17) = -a; 
            ipCoords(1,1)  = -a; ipCoords(1,17) =  a; 
            ipCoords(2,1)  = -a; ipCoords(2,17) =  0;
            //                                     
            ipCoords(0,2)  =  a; ipCoords(0,18) =  a; 
            ipCoords(1,2)  =  a; ipCoords(1,18) = -a; 
            ipCoords(2,2)  = -a; ipCoords(2,18) =  0;
            //                                      
            ipCoords(0,3)  = -a; ipCoords(0,19) =  a; 
            ipCoords(1,3)  =  a; ipCoords(1,19) =  a; 
            ipCoords(2,3)  = -a; ipCoords(2,19) =  0;

            ipCoords(0,4)  = -a; ipCoords(0,20) =  0; 
            ipCoords(1,4)  = -a; ipCoords(1,20) =  0; 
            ipCoords(2,4)  =  a; ipCoords(2,20) = -a;
                                                   
            ipCoords(0,5)  =  a; ipCoords(0,21) =  0; 
            ipCoords(1,5)  = -a; ipCoords(1,21) =  0; 
            ipCoords(2,5)  =  a; ipCoords(2,21) =  a;
                                                  
            ipCoords(0,6)  =  a; ipCoords(0,22) =  0; 
            ipCoords(1,6)  =  a; ipCoords(1,22) = -a; 
            ipCoords(2,6)  =  a; ipCoords(2,22) =  0;
                                                   
            ipCoords(0,7)  = -a; ipCoords(0,23) =  0; 
            ipCoords(1,7)  =  a; ipCoords(1,23) =  a; 
            ipCoords(2,7)  =  a; ipCoords(2,23) =  0;

            ipCoords(0,8)  =  0; ipCoords(0,24) = -a;  
            ipCoords(1,8)  = -a; ipCoords(1,24) =  0; 
            ipCoords(2,8)  = -a; ipCoords(2,24) =  0;
                                                   
            ipCoords(0,9)  =  0; ipCoords(0,25) =  a; 
            ipCoords(1,9)  = -a; ipCoords(1,25) =  0; 
            ipCoords(2,9)  =  a; ipCoords(2,25) =  0;
                                                  
            ipCoords(0,10) =  0; ipCoords(0,26) =  0; 
            ipCoords(1,10) =  a; ipCoords(1,26) =  0; 
            ipCoords(2,10) = -a; ipCoords(2,26) =  0;
                                                   
            ipCoords(0,11) =  0;  
            ipCoords(1,11) =  a;  
            ipCoords(2,11) =  a; 

            ipCoords(0,12) = -a;
            ipCoords(1,12) =  0;
            ipCoords(2,12) = -a;
                               
            ipCoords(0,13) = -a;
            ipCoords(1,13) =  0;
            ipCoords(2,13) =  a;
                               
            ipCoords(0,14) =  a;
            ipCoords(1,14) =  0;
            ipCoords(2,14) = -a;
                               
            ipCoords(0,15) =  a;
            ipCoords(1,15) =  0;
            ipCoords(2,15) =  a;

            // Integration points by face
            ipFaceCoords[0](0,0) = -1; ipFaceCoords[0](0,5) = -1; 
            ipFaceCoords[0](1,0) = -a; ipFaceCoords[0](1,5) =  0; 
            ipFaceCoords[0](2,0) = -a; ipFaceCoords[0](2,5) =  a;
            //                                                        
            ipFaceCoords[0](0,1) = -1; ipFaceCoords[0](0,6) = -1; 
            ipFaceCoords[0](1,1) =  a; ipFaceCoords[0](1,6) = -a; 
            ipFaceCoords[0](2,1) = -a; ipFaceCoords[0](2,6) =  0; 
            //                                                      
            ipFaceCoords[0](0,2) = -1; ipFaceCoords[0](0,7) = -1; 
            ipFaceCoords[0](1,2) =  a; ipFaceCoords[0](1,7) =  a; 
            ipFaceCoords[0](2,2) =  a; ipFaceCoords[0](2,7) =  0; 
            //                                                        
            ipFaceCoords[0](0,3) = -1; ipFaceCoords[0](0,8) = -1; 
            ipFaceCoords[0](1,3) = -a; ipFaceCoords[0](1,8) =  0; 
            ipFaceCoords[0](2,3) =  a; ipFaceCoords[0](2,8) =  0; 
            //
            ipFaceCoords[0](0,4) = -1; 
            ipFaceCoords[0](1,4) =  0; 
            ipFaceCoords[0](2,4) = -a; 
            
            ipFaceCoords[1](0,0) =  1; ipFaceCoords[1](0,5) =  1; 
            ipFaceCoords[1](1,0) = -a; ipFaceCoords[1](1,5) =  0; 
            ipFaceCoords[1](2,0) = -a; ipFaceCoords[1](2,5) =  a;
            //                                                        
            ipFaceCoords[1](0,1) =  1; ipFaceCoords[1](0,6) =  1; 
            ipFaceCoords[1](1,1) =  a; ipFaceCoords[1](1,6) = -a; 
            ipFaceCoords[1](2,1) = -a; ipFaceCoords[1](2,6) =  0; 
            //                                                      
            ipFaceCoords[1](0,2) =  1; ipFaceCoords[1](0,7) =  1; 
            ipFaceCoords[1](1,2) =  a; ipFaceCoords[1](1,7) =  a; 
            ipFaceCoords[1](2,2) =  a; ipFaceCoords[1](2,7) =  0; 
            //                                                        
            ipFaceCoords[1](0,3) =  1; ipFaceCoords[1](0,8) =  1; 
            ipFaceCoords[1](1,3) = -a; ipFaceCoords[1](1,8) =  0; 
            ipFaceCoords[1](2,3) =  a; ipFaceCoords[1](2,8) =  0; 
            //
            ipFaceCoords[1](0,4) =  1; 
            ipFaceCoords[1](1,4) =  0; 
            ipFaceCoords[1](2,4) = -a; 
           
            ipFaceCoords[2](0,0) = -a; ipFaceCoords[2](0,5) =  0; 
            ipFaceCoords[2](1,0) = -1; ipFaceCoords[2](1,5) = -1; 
            ipFaceCoords[2](2,0) = -a; ipFaceCoords[2](2,5) =  a;
            //                                                    
            ipFaceCoords[2](0,1) = -a; ipFaceCoords[2](0,6) = -a; 
            ipFaceCoords[2](1,1) = -1; ipFaceCoords[2](1,6) = -1; 
            ipFaceCoords[2](2,1) =  a; ipFaceCoords[2](2,6) =  0; 
            //                                                    
            ipFaceCoords[2](0,2) =  a; ipFaceCoords[2](0,7) =  a; 
            ipFaceCoords[2](1,2) = -1; ipFaceCoords[2](1,7) = -1; 
            ipFaceCoords[2](2,2) =  a; ipFaceCoords[2](2,7) =  0; 
            //                                                    
            ipFaceCoords[2](0,3) =  a; ipFaceCoords[2](0,8) =  0; 
            ipFaceCoords[2](1,3) = -1; ipFaceCoords[2](1,8) = -1; 
            ipFaceCoords[2](2,3) = -a; ipFaceCoords[2](2,8) =  0; 
            //                         
            ipFaceCoords[2](0,4) =  0; 
            ipFaceCoords[2](1,4) = -1; 
            ipFaceCoords[2](2,4) = -a; 
            
            ipFaceCoords[3](0,0) = -a; ipFaceCoords[3](0,5) =  0; 
            ipFaceCoords[3](1,0) =  1; ipFaceCoords[3](1,5) =  1; 
            ipFaceCoords[3](2,0) = -a; ipFaceCoords[3](2,5) =  a;
            //                                                    
            ipFaceCoords[3](0,1) = -a; ipFaceCoords[3](0,6) = -a; 
            ipFaceCoords[3](1,1) =  1; ipFaceCoords[3](1,6) =  1; 
            ipFaceCoords[3](2,1) =  a; ipFaceCoords[3](2,6) =  0; 
            //                                                    
            ipFaceCoords[3](0,2) =  a; ipFaceCoords[3](0,7) =  a; 
            ipFaceCoords[3](1,2) =  1; ipFaceCoords[3](1,7) =  1; 
            ipFaceCoords[3](2,2) =  a; ipFaceCoords[3](2,7) =  0; 
            //                                                    
            ipFaceCoords[3](0,3) =  a; ipFaceCoords[3](0,8) =  0; 
            ipFaceCoords[3](1,3) =  1; ipFaceCoords[3](1,8) =  1; 
            ipFaceCoords[3](2,3) = -a; ipFaceCoords[3](2,8) =  0; 
            //                         
            ipFaceCoords[3](0,4) =  0; 
            ipFaceCoords[3](1,4) =  1; 
            ipFaceCoords[3](2,4) = -a; 
            
            ipFaceCoords[4](0,0) = -a; ipFaceCoords[4](0,5) =  0;
            ipFaceCoords[4](1,0) = -a; ipFaceCoords[4](1,5) =  a;
            ipFaceCoords[4](2,0) = -1; ipFaceCoords[4](2,5) = -1;
            //                          
            ipFaceCoords[4](0,1) =  a; ipFaceCoords[4](0,6) = -a;
            ipFaceCoords[4](1,1) = -a; ipFaceCoords[4](1,6) =  0;
            ipFaceCoords[4](2,1) = -1; ipFaceCoords[4](2,6) = -1;
            //                         
            ipFaceCoords[4](0,2) =  a; ipFaceCoords[4](0,7) =  a;
            ipFaceCoords[4](1,2) =  a; ipFaceCoords[4](1,7) =  0;
            ipFaceCoords[4](2,2) = -1; ipFaceCoords[4](2,7) = -1;
            //                                                    
            ipFaceCoords[4](0,3) = -a; ipFaceCoords[4](0,8) =  0;
            ipFaceCoords[4](1,3) =  a; ipFaceCoords[4](1,8) =  0;
            ipFaceCoords[4](2,3) = -1; ipFaceCoords[4](2,8) = -1;
            //
            ipFaceCoords[4](0,4) =  0;
            ipFaceCoords[4](1,4) = -a;
            ipFaceCoords[4](2,4) = -1;
            
            ipFaceCoords[5](0,0) = -a; ipFaceCoords[5](0,5) =  0;
            ipFaceCoords[5](1,0) = -a; ipFaceCoords[5](1,5) =  a;
            ipFaceCoords[5](2,0) =  1; ipFaceCoords[5](2,5) =  1;
            //                          
            ipFaceCoords[5](0,1) =  a; ipFaceCoords[5](0,6) = -a;
            ipFaceCoords[5](1,1) = -a; ipFaceCoords[5](1,6) =  0;
            ipFaceCoords[5](2,1) =  1; ipFaceCoords[5](2,6) =  1;
            //                         
            ipFaceCoords[5](0,2) =  a; ipFaceCoords[5](0,7) =  a;
            ipFaceCoords[5](1,2) =  a; ipFaceCoords[5](1,7) =  0;
            ipFaceCoords[5](2,2) =  1; ipFaceCoords[5](2,7) =  1;
            //                                                    
            ipFaceCoords[5](0,3) = -a; ipFaceCoords[5](0,8) =  0;
            ipFaceCoords[5](1,3) =  a; ipFaceCoords[5](1,8) =  0;
            ipFaceCoords[5](2,3) =  1; ipFaceCoords[5](2,8) =  1;
            //
            ipFaceCoords[5](0,4) =  0;
            ipFaceCoords[5](1,4) = -a;
            ipFaceCoords[5](2,4) =  1;
            
            // Integration points for all faces ignoring the 
            // redundant dim
            ipFace2DCoords(0,0) = -a; ipFace2DCoords(0,5) =  0; 
            ipFace2DCoords(1,0) = -a; ipFace2DCoords(1,5) =  a; 
            // 
            ipFace2DCoords(0,1) =  a; ipFace2DCoords(0,6) = -a; 
            ipFace2DCoords(1,1) = -a; ipFace2DCoords(1,6) =  0; 
            //
            ipFace2DCoords(0,2) =  a; ipFace2DCoords(0,7) =  a; 
            ipFace2DCoords(1,2) =  a; ipFace2DCoords(1,7) =  0; 
            // 
            ipFace2DCoords(0,3) = -a; ipFace2DCoords(0,8) =  0; 
            ipFace2DCoords(1,3) =  a; ipFace2DCoords(1,8) =  0; 
            // 
            ipFace2DCoords(0,4) =  0; 
            ipFace2DCoords(1,4) = -a; 

            // Integration point weights
            for (int i=0; i<numIPs; ++i) {
              if (i < 8)       ipWeights[i] = 125.0/729.0;
              else if (i < 20) ipWeights[i] = 200.0/729.0;
              else if (i < 26) ipWeights[i] = 320.0/729.0;
              else             ipWeights[i] = 512.0/729.0;
            }

            // Face integration point weights
            for (int i=0; i<numFaceIPs; ++i) {
              if (i < 4)      ipFaceWeights[i] = 25.0/81.0;
              else if (i < 8) ipFaceWeights[i] = 40.0/81.0;
              else            ipFaceWeights[i] = 64.0/81.0;
            }

            break;
          }
          // Error
          default: {
            throw ATC_Error("Unrecognized quadrature type "
                            "for element type HEXA.");
          }
        }
      } else if (geo == TETRA) {
        switch (quad) {
          // 4 point quadratures
          case NODAL:
          case GAUSS2: {
            // Set number of IPs and face IPs
            numIPs = 4;
            numFaceIPs = 3;

            // Size matrices and vectors thereof accordingly
            ipCoords.resize(4,numIPs);
            ipFaceCoords.assign(4,DENS_MAT(3,numFaceIPs));
            ipFace2DCoords.resize(3,numFaceIPs);
            ipWeights.reset(numIPs);
            ipFaceWeights.reset(numFaceIPs);

            double v1, v2, a1, a2;
            v1 = 0.585410196624969;
            v2 = 0.138196601125011;
            a1 = 2.0/3.0;
            a2 = 1.0/6.0;
            if (quad == NODAL) {
              // Nodal quadrature
              v1 = 1;
              v2 = 0;
              a1 = 1;
              a2 = 0;
            }
            
            // Integration point coordinates
            ipCoords.resize(4, numIPs);
            ipCoords(0,0) = v2; ipCoords(0,2) = v2; 
            ipCoords(1,0) = v2; ipCoords(1,2) = v1; 
            ipCoords(2,0) = v2; ipCoords(2,2) = v2; 
            ipCoords(3,0) = v1; ipCoords(3,2) = v2;
            //
            ipCoords(0,1) = v1; ipCoords(0,3) = v2; 
            ipCoords(1,1) = v2; ipCoords(1,3) = v2; 
            ipCoords(2,1) = v2; ipCoords(2,3) = v1; 
            ipCoords(3,1) = v2; ipCoords(3,3) = v2;
            
            // Integration points by face
            
            // ...face 0                  ...face 2
            ipFaceCoords[0](0,0) = a1; ipFaceCoords[2](0,0) = a2;
            ipFaceCoords[0](1,0) = a2; ipFaceCoords[2](1,0) = 0; 
            ipFaceCoords[0](2,0) = a2; ipFaceCoords[2](2,0) = a2;
            //                                                   
            ipFaceCoords[0](0,1) = a2; ipFaceCoords[2](0,1) = a1;
            ipFaceCoords[0](1,1) = a1; ipFaceCoords[2](1,1) = 0; 
            ipFaceCoords[0](2,1) = a2; ipFaceCoords[2](2,1) = a2;
            //                                                   
            ipFaceCoords[0](0,2) = a2; ipFaceCoords[2](0,2) = a2;
            ipFaceCoords[0](1,2) = a2; ipFaceCoords[2](1,2) = 0; 
            ipFaceCoords[0](2,2) = a1; ipFaceCoords[2](2,2) = a1;
                                                                 
            // ...face 1                  ...face 3
            ipFaceCoords[1](0,0) = 0;  ipFaceCoords[3](0,0) = a2;
            ipFaceCoords[1](1,0) = a1; ipFaceCoords[3](1,0) = a2;
            ipFaceCoords[1](2,0) = a2; ipFaceCoords[3](2,0) = 0;
            //                                                   
            ipFaceCoords[1](0,1) = 0;  ipFaceCoords[3](0,1) = a2;
            ipFaceCoords[1](1,1) = a2; ipFaceCoords[3](1,1) = a1;
            ipFaceCoords[1](2,1) = a2; ipFaceCoords[3](2,1) = 0;
            //                                                   
            ipFaceCoords[1](0,2) = 0;  ipFaceCoords[3](0,2) = a1;
            ipFaceCoords[1](1,2) = a2; ipFaceCoords[3](1,2) = a2;
            ipFaceCoords[1](2,2) = a1; ipFaceCoords[3](2,2) = 0;
                                                                                       
            // 2D integration points for faces
            ipFace2DCoords(0,0) = a2; ipFace2DCoords(0,2) = a2; 
            ipFace2DCoords(1,0) = a2; ipFace2DCoords(1,2) = a1; 
            ipFace2DCoords(2,0) = a1; ipFace2DCoords(2,2) = a2;  
            //
            ipFace2DCoords(0,1) = a1; 
            ipFace2DCoords(1,1) = a2; 
            ipFace2DCoords(2,1) = a2;

            // Integration point weights
            ipWeights = (1.0/6.0)/numIPs;
            
            // Integration point face weights
            ipFaceWeights = (1.0/2.0)/numFaceIPs;

            break;
          }
          case GAUSS3: {
            // Set number of IPs and face IPs
            numIPs = 5;
            numFaceIPs = 4;

            // Size matrices and vectors thereof accordingly
            ipCoords.resize(4,numIPs);
            ipFaceCoords.assign(4,DENS_MAT(3,numFaceIPs));
            ipFace2DCoords.resize(3,numFaceIPs);
            ipWeights.reset(numIPs);
            ipFaceWeights.reset(numFaceIPs);

            double v1, v2, v3, a1, a2, a3;
            /* These weights for calculating Gaussian Quadrature 
             * points are taken from the paper "Integration Points 
             * for Triangles and Tetrahedra Obtained from the 
             * Gaussian Quadrature Points for a Line" by K. S. Sunder
             * and R. A. Cookson, Computers and Structures, Vol 21, 
             * No. 5, 1985. */
            v1 = 1.0/4.0;
            v2 = 1.0/2.0;
            v3 = 1.0/6.0;
            a1 = 1.0/3.0;
            a2 = 3.0/5.0;
            a3 = 1.0/5.0;
            
            // Integration point coordinates
            ipCoords.resize(4, numIPs);
            ipCoords(0,0) = v1;
            ipCoords(1,0) = v1;
            ipCoords(2,0) = v1;
            ipCoords(3,0) = v1;
            //
            ipCoords(0,1) = v3; ipCoords(0,3) = v3; 
            ipCoords(1,1) = v3; ipCoords(1,3) = v2; 
            ipCoords(2,1) = v3; ipCoords(2,3) = v3; 
            ipCoords(3,1) = v2; ipCoords(3,3) = v3;
            //
            ipCoords(0,2) = v2; ipCoords(0,4) = v3; 
            ipCoords(1,2) = v3; ipCoords(1,4) = v3; 
            ipCoords(2,2) = v3; ipCoords(2,4) = v2; 
            ipCoords(3,2) = v3; ipCoords(3,4) = v3;
            
            // Integration points by face
            
            // ...face 0                  ...face 2
            ipFaceCoords[0](0,0) = a1; ipFaceCoords[2](0,0) = a1;
            ipFaceCoords[0](1,0) = a1; ipFaceCoords[2](1,0) = 0; 
            ipFaceCoords[0](2,0) = a1; ipFaceCoords[2](2,0) = a1;
            //                                                   
            ipFaceCoords[0](0,1) = a2; ipFaceCoords[2](0,1) = a3;
            ipFaceCoords[0](1,1) = a3; ipFaceCoords[2](1,1) = 0; 
            ipFaceCoords[0](2,1) = a3; ipFaceCoords[2](2,1) = a3;
            //                                                   
            ipFaceCoords[0](0,2) = a3; ipFaceCoords[2](0,2) = a2;
            ipFaceCoords[0](1,2) = a2; ipFaceCoords[2](1,2) = 0; 
            ipFaceCoords[0](2,2) = a3; ipFaceCoords[2](2,2) = a3;
            //                                                   
            ipFaceCoords[0](0,3) = a3; ipFaceCoords[2](0,3) = a3;
            ipFaceCoords[0](1,3) = a3; ipFaceCoords[2](1,3) = 0; 
            ipFaceCoords[0](2,3) = a2; ipFaceCoords[2](2,3) = a2;
                                                                 
            // ...face 1                  ...face 3
            ipFaceCoords[1](0,0) = 0;  ipFaceCoords[3](0,0) = a1;
            ipFaceCoords[1](1,0) = a1; ipFaceCoords[3](1,0) = a1;
            ipFaceCoords[1](2,0) = a1; ipFaceCoords[3](2,0) = 0;
            //                                                   
            ipFaceCoords[1](0,1) = 0;  ipFaceCoords[3](0,1) = a3;
            ipFaceCoords[1](1,1) = a2; ipFaceCoords[3](1,1) = a3;
            ipFaceCoords[1](2,1) = a3; ipFaceCoords[3](2,1) = 0;
            //                                                   
            ipFaceCoords[1](0,2) = 0;  ipFaceCoords[3](0,2) = a3;
            ipFaceCoords[1](1,2) = a3; ipFaceCoords[3](1,2) = a2;
            ipFaceCoords[1](2,2) = a3; ipFaceCoords[3](2,2) = 0;
            //                                                   
            ipFaceCoords[1](0,3) = 0;  ipFaceCoords[3](0,3) = a2;
            ipFaceCoords[1](1,3) = a3; ipFaceCoords[3](1,3) = a3;
            ipFaceCoords[1](2,3) = a2; ipFaceCoords[3](2,3) = 0;
                                                                                       
            // 2D integration points for faces
            //
            ipFace2DCoords(0,0) = a1; ipFace2DCoords(0,2) = a2; 
            ipFace2DCoords(1,0) = a1; ipFace2DCoords(1,2) = a3; 
            ipFace2DCoords(2,0) = a1; ipFace2DCoords(2,2) = a3;
            //
            ipFace2DCoords(0,1) = a3; ipFace2DCoords(0,3) = a3; 
            ipFace2DCoords(1,1) = a3; ipFace2DCoords(1,3) = a2; 
            ipFace2DCoords(2,1) = a2; ipFace2DCoords(2,3) = a3;  

            for (int i=0; i<numIPs; ++i) {
              if (i < 1)       ipWeights[i] = -4.0/5.0;
              else             ipWeights[i] = 9.0/20.0;
            }

            // Face integration point weights
            for (int i=0; i<numFaceIPs; ++i) {
              if (i < 1)      ipFaceWeights[i] = -9.0/16.0;
              else            ipFaceWeights[i] = 25.0/48.0;
            }

            break;
          }
          // Error
          default: {
            throw ATC_Error("Unrecognized quadrature type "
                            "for element type TETRA.");
          }
        }
      }
    }
  };

}; // namespace ATC

#endif // FE_QUADRATURE_H
