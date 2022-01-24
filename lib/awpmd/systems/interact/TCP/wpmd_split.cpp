# include "wpmd_split.h"
//# include "erf.h"


void AWPMD_split::resize(int flag){
  for(int s=0;s<2;s++){
    wf_norm[s].resize(ne[s]);
    if(flag&(0x8|0x4)){ //electron forces needed
      wf_norm_der[s].resize(nvar[s]);
      ovl_der[s].resize(nvar[s]);
      E_der[s].resize(nvar[s]);

      if(flag&(0x8|0x4) || norm_needed){ //electron forces or norm matrix needed
        if(approx==HARTREE){ // L and Norm are needed in block form
          Lh[s].resize(nvar[s]);
          if(norm_needed){
            Normh[s].resize(ne[s]);
            for(int i=0;i<ne[s];i++)
              Normh[s][i].init(10*nspl[s][i]);
          }
        }
        else if(norm_needed){
          Norm[s].init(nvar[s]);
        }
      }
    }
  }
  AWPMD::resize(flag);
}


int AWPMD_split::add_split(Vector_3 &x, Vector_3 &v, double &w, double &pw, Vector_2 &c, double mass, double q, int tag){
  if(!spl_add){
    nspl[s_add].push_back(1);
    ne[s_add]++;
  }
  else{
    nspl[s_add][ne[s_add]-1]++; // incrementing the WP number for the last electron
  }
  spl_add++;
  nwp[s_add]++;
  nvar[s_add]+=10;
  wp[s_add].push_back(create_wp(x,v,w,pw,mass));
  split_c[s_add].push_back(c);
  qe[s_add].push_back(q);

  if(tag==0)
    tag=spl_add;
  partition1[s_add].push_back(tag);

  return spl_add-1;
}

int AWPMD_split::set_electrons(int s, int nel, Vector_3P x, Vector_3P v, double* w, double* pw, Vector_2 *c, int *splits, double mass, double *q)
{
  if(s < 0 || s > 1)
    return LOGERR(-1,logfmt("AWPMD_split.set_electrons: invaid spin setting (%d)!",s),LINFO);

  // calculating the total n
  nvar[s]=0;
  int n=0;
  for(int i=0;i<nel;i++){
    n+=splits[i];
    nvar[s]+=10*splits[i]; // number of dynamic variables per wp: x[3],p[3],w,pw,c_re,c_im
  }
  nwp[s]=n;

  norm_matrix_state[s] = NORM_UNDEFINED;
  ne[s]=nel;
  wp[s].resize(n);

  split_c[s].resize(n);
  split_c[s].assign(c,c+n);


  nspl[s].resize(nel);
  nspl[s].assign(splits,splits+nel);

  partition1[s].clear();
  for(int i=0;i<n;i++){

    /*if(constraint==FIX){
      w[i]=w0;
      pw[i]=0.;
    }

    double rw;
    if(Lextra>0){ // width PBC, keeping the width are within [0,Lextra]
      w[i]=fmod(w[i],Lextra);
      if(w[i]<0) w[i]+=Lextra;
      rw=w[i]; // WP width for energy evaluation is within [0, L/2]
      if(rw > Lextra/2) rw = Lextra - rw;
    }
    else
      rw=w[i];

    wp[s][i].init(rw,x[i],v[i]*m_electron/h_plank,pw[i]/h_plank);*/
    wp[s][i]=create_wp(x[i],v[i],w[i],pw[i],mass);
    //printf("%15d %15g\n",i,rw);
    // assign default partition
    partition1[s].push_back(i+1);
  }

  // assign electronic charge
  if(q)
    qe[s].assign(q,q+nwp[s]);
  else
    qe[s].assign(nwp[s],-1);



  return 1;
}



void AWPMD_split::eterm_deriv(int ic1,int s1,int c1, int j1,int ic2,int s2, int c2, int k2,cdouble pref,
                              const OverlapDeriv &o,cdouble v,cdouble dv_aj_conj,
                              cdouble dv_ak,cVector_3 dv_bj_conj, cVector_3 dv_bk){
  cdouble cj(split_c[s1][ic1+j1][0],split_c[s1][ic1+j1][1]);
  cdouble ck(split_c[s2][ic2+k2][0],split_c[s2][ic2+k2][1]);
  int indw1=8*ic1, indw2=8*ic2;
  int indn1=(nvar[s1]/10)*8+2*ic1, indn2=(nvar[s2]/10)*8+2*ic2;
  cdouble part_jk=conj(cj)*ck;

  int M= 1; //(j1==k2 ? 1 : 2);

  // over a_k_re
  E_der[s2][indw2+8*k2]+=   M*real(pref*part_jk*(o.da2_re()*v+o.I0*dv_ak));
  // over a_k_im
  E_der[s2][indw2+8*k2+1]+=   M*real(pref*part_jk*(o.da2_im()*v+i_unit*o.I0*dv_ak));
  // over a_j_re
  E_der[s1][indw1+8*j1]+=   M*real(pref*part_jk*(o.da1_re()*v+o.I0*dv_aj_conj));
  // over a_j_im
  E_der[s1][indw1+8*j1+1]+= M*real(pref*part_jk*(o.da1_im()*v-i_unit*o.I0*dv_aj_conj));

  for(int i=0;i<3;i++){
    // over b_k_re
    E_der[s2][indw2+8*k2+2+2*i]+= M*real(pref*part_jk*(o.db2_re(i)*v+o.I0*dv_bk[i]));
    // over b_k_im
    E_der[s2][indw2+8*k2+2+2*i+1]+= M*real(pref*part_jk*(o.db2_im(i)*v+i_unit*o.I0*dv_bk[i]));
    // over b_j_re
    E_der[s1][indw1+8*j1+2+2*i]+= M*real(pref*part_jk*(o.db1_re(i)*v+o.I0*dv_bj_conj[i]));
    // over b_j_im
    E_der[s1][indw1+8*j1+2+2*i+1]+= M*real(pref*part_jk*(o.db1_im(i)*v-i_unit*o.I0*dv_bj_conj[i]));
  }


  // over ck_re
  E_der[s2][indn2+2*k2]+=M*real(pref*conj(cj)*o.I0*v);
  // over ck_im
  E_der[s2][indn2+2*k2+1]+=M*real(pref*i_unit*conj(cj)*o.I0*v);
  // over cj_re
  E_der[s1][indn1+2*j1]+=M*real(pref*ck*o.I0*v);
  // over cj_im
  E_der[s1][indn1+2*j1+1]+=M*real(-pref*i_unit*ck*o.I0*v);

  double t= -M*real(pref*part_jk*o.I0*v);
  // nonlocal terms: TODO: make a separate global loop for summation of nonlocal terms
  for(int j=0;j<nspl[s1][c1];j++){
    for(int i=0;i<8;i++)
      E_der[s1][indw1+8*j+i]+=t*wf_norm_der[s1][indw1+8*j+i];
    E_der[s1][indn1+2*j  ]+=t*wf_norm_der[s1][indn1+2*j  ];
    E_der[s1][indn1+2*j+1]+=t*wf_norm_der[s1][indn1+2*j+1];
  }
  for(int k=0;k<nspl[s2][c2];k++){
    for(int i=0;i<8;i++)
      E_der[s2][indw2+8*k+i]+=t*wf_norm_der[s2][indw2+8*k+i];
    E_der[s2][indn2+2*k]+=   t*wf_norm_der[s2][indn2+2*k  ];
    E_der[s2][indn2+2*k+1]+= t*wf_norm_der[s2][indn2+2*k+1];
  }

}


void AWPMD_split::calc_norms(int flag){
  if(flag&0x4){ // electron forces requested
    for(int s1=0;s1<2;s1++){ // clearing norm derivatives
      for(int i=0;i<nvar[s1];i++){
        wf_norm_der[s1][i]=0;
        E_der[s1][i]=0;
        ovl_der[s1][i]=0;
      }
    }
  }
  // calculating block norms and derivatives
  for(int s1=0;s1<2;s1++){
    int ic1=0; // starting index of the wp for current electron
    int indw1=0; // starting index of the electron wp coordinates
    int indn1=(nvar[s1]/10)*8; // starting index of the electron norm coordinates

    for(int c1=0;c1<ne[s1];c1++){

      // calculating the block norm
      wf_norm[s1][c1]=0.;
      for(int j1=0;j1<nspl[s1][c1];j1++){
        double cj_re=split_c[s1][ic1+j1][0];
        double cj_im=split_c[s1][ic1+j1][1];
        cdouble ccj=cdouble(cj_re,-cj_im);
        double part_jj=norm(ccj);
        wf_norm[s1][c1]+=part_jj;
        WavePacket& wj = wp[s1][ic1+j1];
        OverlapDeriv o;
        if(flag&(0x8|0x4)){ //electron forces needed
          wf_norm_der[s1][indn1+2*j1]+=2*cj_re;  // over cj_re
          wf_norm_der[s1][indn1+2*j1+1]+=2*cj_im; // over cj_im
          o.set1(wj);// conjugate: mu -> -mu, v -> -v !!!
        }


        for(int k1=j1+1;k1<nspl[s1][c1];k1++){
          double ck_re=split_c[s1][ic1+k1][0];
          double ck_im=split_c[s1][ic1+k1][1];
          cdouble ck=cdouble(ck_re,ck_im);

          WavePacket wk=wp[s1][ic1+k1];
          if(pbc)
            move_to_image(wj,wk);
          WavePacket wjk=conj(wj)*wk;
          cdouble I0=wjk.integral();

          cdouble part_jk=ccj*ck;
          wf_norm[s1][c1]+=2*real(part_jk*I0);


          if(flag&(0x8|0x4)){ //electron forces needed
            o.set2(wk,&I0);


            wf_norm_der[s1][indw1+8*k1]+=   2*real(part_jk*o.da2_re());           // over a_k_re
            wf_norm_der[s1][indw1+8*k1+1]+=   2*real(part_jk*o.da2_im());    // over a_k_im

            wf_norm_der[s1][indw1+8*j1]+=   2*real(part_jk*o.da1_re());             // over a_j_re
            wf_norm_der[s1][indw1+8*j1+1]+= 2*real(part_jk*o.da1_im());    // over a_j_im

            for(int i=0;i<3;i++){
              wf_norm_der[s1][indw1+8*k1+2+2*i]+= 2*real(part_jk*o.db2_re(i));           // over b_k_re
              wf_norm_der[s1][indw1+8*k1+2+2*i+1]+= 2*real(part_jk*o.db2_im(i));  // over b_k_im

              wf_norm_der[s1][indw1+8*j1+2+2*i]+= 2*real(part_jk*o.db1_re(i));           // over b_j_re
              wf_norm_der[s1][indw1+8*j1+2+2*i+1]+= 2*real(part_jk*o.db1_im(i));  // over b_j_im
            }

            wf_norm_der[s1][indn1+2*j1]+=2*real(ck*I0);  // over cj_re
            wf_norm_der[s1][indn1+2*j1+1]+=2*real(-i_unit*ck*I0);  // over cj_im
            wf_norm_der[s1][indn1+2*k1]+=2*real(ccj*I0); // over ck_re
            wf_norm_der[s1][indn1+2*k1+1]+=2*real(i_unit*ccj*I0); // over ck_im


            // overlap derivatives (for norm matrix)
            ovl_der[s1][indw1+8*k1]+=   ccj*o.da2_re();           // over a_k_re
            ovl_der[s1][indw1+8*k1+1]+=   ccj*o.da2_im();    // over a_k_im

            ovl_der[s1][indw1+8*j1]+=   ck*o.da1_re();             // over a_j_re
            ovl_der[s1][indw1+8*j1+1]+= ck*o.da1_im();    // over a_j_im

            for(int i=0;i<3;i++){
              ovl_der[s1][indw1+8*k1+2+2*i]+= ccj*o.db2_re(i);           // over b_k_re
              ovl_der[s1][indw1+8*k1+2+2*i+1]+= ccj*o.db2_im(i);  // over b_k_im

              ovl_der[s1][indw1+8*j1+2+2*i]+= ck*o.db1_re(i);           // over b_j_re
              ovl_der[s1][indw1+8*j1+2+2*i+1]+= ck*o.db1_im(i);  // over b_j_im
            }

            ovl_der[s1][indn1+2*j1]+=ck*I0;  // over cj_re
            ovl_der[s1][indn1+2*j1+1]+=-i_unit*ck*I0;  // over cj_im
            ovl_der[s1][indn1+2*k1]+=ccj*I0; // over ck_re
            ovl_der[s1][indn1+2*k1+1]+=i_unit*ccj*I0; // over ck_im
          }
        } // k1
      }// j1
      // normalizing the norm derivative
      for(int j1=0;j1<nspl[s1][c1];j1++){
        for(int i=0;i<8;i++) // wp parameters
          wf_norm_der[s1][indw1+8*j1+i]/=2*wf_norm[s1][c1];
        // c
        wf_norm_der[s1][indn1+2*j1]/=2*wf_norm[s1][c1];
        wf_norm_der[s1][indn1+2*j1+1]/=2*wf_norm[s1][c1];
      }
      //wf_norm[s1][c1]=1;
      //wf_norm[s1][c1]=sqrt(wf_norm[s1][c1]);
      ic1+=nspl[s1][c1];
      indw1+=8*nspl[s1][c1]; // 8 variables in each wavepacket
      indn1+=2*nspl[s1][c1]; // 2 variables in each wp norm
    }// c1
  }// s1
}

void AWPMD_split::clear_forces(int flag,Vector_3P fi, Vector_3P fe_x,
                               Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c){
  if(flag&0x1){
    for(int i=0;i<ni;i++)
      fi[i]=Vector_3(0.);
  }
  if(flag&0x4 && !(flag&0x10)){ // electron forces requested in physical representation
    int iv1=0;
    for(int s1=0;s1<2;s1++){ // clearing forces
      for(int c1=0;c1<ne[s1];c1++){
        for(int j1=0;j1<nspl[s1][c1];j1++){
          fe_x[iv1+j1]=Vector_3(0,0,0);
          fe_p[iv1+j1]=Vector_3(0,0,0);
          fe_w[iv1+j1]=0;
          fe_pw[iv1+j1]=0;
          fe_c[iv1+j1]=Vector_2(0,0);
        }
        iv1+=nspl[s1][c1];
      }
    }
  }
}

void AWPMD_split::get_el_forces(int flag, Vector_3P fe_x,
                               Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c){
  if(flag&0x4) //need to replace the forces
    clear_forces(0x4,nullptr,fe_x,fe_p,fe_w,fe_pw,fe_c);

  // recalculating derivatives
  if(flag&(0x8|0x4)){ //electron forces needed
    int iv1=0;
    for(int s1=0;s1<2;s1++){
      int ic1=0; // starting index of the wp for current electron
      for(int c1=0;c1<ne[s1];c1++){
        int indw1=8*ic1;
        int indn1=(nvar[s1]/10)*8+2*ic1;
        for(int k1=0;k1<nspl[s1][c1];k1++){
          WavePacket wk=wp[s1][ic1+k1];
          /*double w=wk.get_width();
          Vector_3 r=wk.get_r();
          double t=3/(2*w*w*w);
          fe_w[ic1+k1]+= t*E_der[s1][indw1+8*k1]+imag(wk.a)*E_der[s1][indw1+8*k1+1]/w;
          fe_pw[ic1+k1]+=E_der[s1][indw1+8*k1+1]/(2*w*h_plank);
          for(int i=0;i<3;i++){
            fe_x[ic1+k1][i]+= -2*real(wk.a)*E_der[s1][indw1+8*k1+2+2*i]-2*imag(wk.a)*E_der[s1][indw1+8*k1+2+2*i+1];
            fe_p[ic1+k1][i]+= (-E_der[s1][indw1+8*k1+2+2*i+1])*(m_electron/h_plank); // *(h_plank/m_electron);
            fe_pw[ic1+k1]+=(r[i]*E_der[s1][indw1+8*k1+2+2*i+1]/w)/h_plank;
            fe_w[ic1+k1]+=2*r[i]*(t*E_der[s1][indw1+8*k1+2+2*i]+imag(wk.a)*E_der[s1][indw1+8*k1+2+2*i+1]/w);
          }*/
          wk.int2phys_der< minus >(E_der[s1].begin()+indw1,(double *)&fe_x[iv1+k1],(double *)&fe_p[iv1+k1],&fe_w[iv1+k1],&fe_pw[iv1+k1], 1./one_h);
          fe_c[iv1+k1]+=-Vector_2(E_der[s1][indn1+2*k1],E_der[s1][indn1+2*k1+1]);
        }// k1
        ic1+=nspl[s1][c1]; // incrementing block1 wp address
        iv1+=nspl[s1][c1]; // incrementing global variable address
      }
    }
  }
}

//e same as interaction, but using Hartee factorization (no antisymmetrization)
int  AWPMD_split::interaction_hartree(int flag, Vector_3P fi, Vector_3P fe_x,
                                      Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c){

  // resize arrays if needed
  enum APPROX tmp=HARTREE;
  swap(tmp,approx); // do not need large matrices
  resize(flag);
  swap(tmp,approx);

  // clearing forces
  clear_forces(flag,fi,fe_x,fe_p,fe_w,fe_pw,fe_c);
  // calculate block norms and (optionally) derivatives
  calc_norms(flag);

  Eee = Ew = 0.;
  for(int s1=0;s1<2;s1++){
    Ee[s1]=0.;
    Eei[s1]=0.;
    int ic1=0; // starting index of the wp for current electron

    for(int c1=0;c1<ne[s1];c1++){
      // calculating single-electron quantities within block
      double pref=-h2_me/(2*wf_norm[s1][c1]); // ekin
      double pref_ei=coul_pref/wf_norm[s1][c1];

      for(int j1=0;j1<nspl[s1][c1];j1++){
        cdouble cj(split_c[s1][ic1+j1][0],split_c[s1][ic1+j1][1]);
        WavePacket wj=wp[s1][ic1+j1];

        OverlapDeriv o;
        if(flag&(0x8|0x4)) //electron forces needed
          o.set1(wj);

        for(int k1=j1;k1<nspl[s1][c1];k1++){
          int M1a=(j1==k1 ? 1: 2);
          double M1e, M1f;
          _mytie(M1e,M1f)=check_part1(s1,ic1+j1,ic1+k1)*M1a;

          cdouble ck(split_c[s1][ic1+k1][0],split_c[s1][ic1+k1][1]);

          // electrons kinetic energy
          WavePacket wk=wp[s1][ic1+k1];
          if(pbc)
            move_to_image(wj,wk);

          WavePacket wjk=conj(wj)*wk;
          cdouble I0 = wjk.integral();
          cdouble part_jk=conj(cj)*ck;

          if(M1e){
            cVector_3 v1=conj(wj.b)*wk.a-wk.b*conj(wj.a);
            cdouble v=(v1*v1)/wjk.a;
            v-=6.*wk.a*conj(wj.a);
            v/=wjk.a;

            //v=1.;

            // kinetic energy contribution
            Ee[s1] += M1e*real(part_jk*I0*v)*pref;  // Ejk+Ekj=Ejk+conj(Ejk), j!=k or Ejj



            if(flag&(0x8|0x4)){ //electron forces needed
              cVector_3 tv=wk.b*conj(wj.a)-conj(wj.b)*wk.a;
              cdouble ajk2=wjk.a*wjk.a;
              cdouble ajk3=ajk2*wjk.a;
              cdouble dv_aj_conj=     -2*wk.a*(3*wjk.a*wk.a-tv*wjk.b)/ajk3;
              cdouble dv_ak=-2*conj(wj.a)*(3*wjk.a*conj(wj.a)+tv*wjk.b)/ajk3;
              cVector_3 dv_bj_conj=-2*tv*wk.a/ajk2;
              cVector_3 dv_bk=2*tv*conj(wj.a)/ajk2;

              /*cdouble dv_aj_conj=0.;
              cdouble dv_ak=0.;
              cVector_3 dv_bj_conj;
              cVector_3 dv_bk;*/



              o.set2(wk,&I0);
              // calculate full derivative of the term pref*conj(cj)*ck*Ijk*vjk/sqrt(nrm(s1)*nrm(s2))
              // denominator must be included in pref
              eterm_deriv(ic1,s1,c1,j1,ic1,s1,c1,k1,M1f*pref,o,v,dv_aj_conj,dv_ak,dv_bj_conj,dv_bk);
            }
          }// M1e

          cVector_3 djk=wjk.b/(2.*wjk.a);
          // e-i energy
          //cdouble sum(0.,0.);
          for(int i=0;i<ni;i++){  // ions loop
            double M1ie, M1if;
            _mytie(M1ie,M1if)=check_part1ei(s1,ic1+j1,ic1+k1,i)*M1a;
            if(!M1ie)
              continue;

            cVector_3 gjki=djk - cVector_3(xi[i]);

            if(pbc) // correcting the real part (distance) according to PBC
              gjki=rcell1(gjki,cell,pbc);
              //-Igor- gkli=cVector_3(real(gkli).rcell1(cell,pbc),imag(gkli));

            cdouble ngjki=gjki.norm();
            cdouble sqajk=sqrt(wjk.a);
            //cdouble ttt = cerf_div(ngjki,c);
            cdouble v=cerf_div(ngjki,sqajk);
            double pref_eiq=pref_ei*qi[i]*(qe[s1][ic1+j1]+qe[s1][ic1+k1])*0.5;
            cdouble dE=pref_eiq*I0*v*part_jk;

            // ions energy contribution
            Eei[s1] += M1ie*real(dE);
            //sum+=dE;

            cdouble t1, t2;
            bool zero_force=(fabs(real(ngjki))+fabs(imag(ngjki))<1e-10);
            if(flag&(0x8|0x4|0x3) ){ //some derivatives needed
              cdouble arg=ngjki*sqajk;
              t1=two_over_sqr_pi*exp(-arg*arg);
              t2=t1/(2*sqajk);
              t1*=sqajk;
            }

            if(flag&0x3 && !zero_force){// calculate forces on ions
              cdouble dEw=pref_eiq*I0*t1*part_jk;
              dEw=(dE-dEw)/ngjki;
              Vector_3 dir=-real(gjki);
              dir.normalize();
              fi[i]+=M1if*real( dEw)*dir;
            }
            if(flag&(0x8|0x4)){ //electron forces needed
              cdouble dv_ak= (zero_force ?  0.: (djk*gjki)*(v-t1)/(wjk.a*ngjki*ngjki)) +t2;
              cVector_3 dv_bk= zero_force ? cVector_3() :  -gjki*(v-t1)/(2*wjk.a*ngjki*ngjki);
              eterm_deriv(ic1,s1,c1,j1,ic1,s1,c1,k1,M1if*pref_eiq,o,v,dv_ak,dv_ak,dv_bk,dv_bk);
            }
          }

          // extra constraint energy
          if(j1==k1 && constraint == HARM && M1e){
            cdouble v=conj(wj.a)*wj.a;
            Ew += harm_w0_4 * real(part_jk*v)/wf_norm[s1][c1];
            if(flag&(0x8|0x4)){ //electron forces needed
              cdouble dv_ak=     conj(wj.a);
              cdouble dv_aj_conj =wj.a;
              eterm_deriv(ic1,s1,c1,j1,ic1,s1,c1,k1,harm_w0_4/wf_norm[s1][c1],o,v,dv_aj_conj,dv_ak,cVector_3(),cVector_3());
            }
          }
# if 1
          // second block
          // e-e interaction
          for(int s2=s1;s2<2;s2++){
            int ic2=0; // starting index of the wp for current electron

            for(int c2=0 ;c2<ne[s2];ic2+=nspl[s2][c2],c2++){ // incrementing block2 wp address
              if(s1==s2 && c2<=c1)
                continue;

              double pref_ee=coul_pref/(wf_norm[s1][c1]*wf_norm[s2][c2]);
              for(int j2=0;j2<nspl[s2][c2];j2++){

                cdouble cj2(split_c[s2][ic2+j2][0],split_c[s2][ic2+j2][1]);
                WavePacket& wj2=wp[s2][ic2+j2];

                OverlapDeriv o2;
                if(flag&(0x8|0x4)) //electron forces needed
                  o2.set1(wj2);

                for(int k2=j2;k2<nspl[s2][c2];k2++){
                  double M2e, M2f;
                  _mytie(M2e,M2f)=check_part1(s1,ic1+j1,ic1+k1);

                  cdouble ck2(split_c[s2][ic2+k2][0],split_c[s2][ic2+k2][1]);

                  WavePacket wk2=wp[s2][ic2+k2];
                  if(pbc)
                    move_to_image(wj2,wk2);
                  WavePacket wjk2=conj(wj2)*wk2;
                  cdouble I02 = wjk2.integral();


                  cdouble part_jk2=conj(cj2)*ck2;

                  cVector_3 djk2=wjk2.b/(2*wjk2.a);
                  cVector_3 ddv=djk-djk2;
                  cdouble dd=ddv.norm();
                  cdouble aa=1./sqrt(1./wjk.a+1./wjk2.a);
                  //double ww1=wj.get_width();
                  //double ww2=wj2.get_width();
                  cdouble v=cerf_div(dd,aa);
                  double pref_eeq=pref_ee*(qe[s1][ic1+j1]+qe[s1][ic1+k1])*(qe[s2][ic2+j2]+qe[s2][ic2+k2])*0.25;
                  cdouble Vj1j2k1k2=pref_eeq*I0*I02*v*part_jk*part_jk2;
                  Eee+=M1e*M2e*real(Vj1j2k1k2);
                  //Eee+=(j1==k1 || j2==k2 ? 1: 2)*real(Vj1j2k1k2);


                  cdouble t1, t2;
                  bool zero_force=(fabs(real(dd))+fabs(imag(dd))<1e-10);
                  if(flag&(0x8|0x4) ){ //electron forces needed
                    cdouble arg=dd*aa;
                    t1=two_over_sqr_pi*exp(-arg*arg)*aa;
                    t2=t1*aa*aa/2;


                    cdouble dv_ak1= (zero_force ?  0.: (djk*ddv)*(v-t1)/(wjk.a*dd*dd)) +t2/(wjk.a*wjk.a);
                    cVector_3 dv_bk1= zero_force ? cVector_3() :  -ddv*(v-t1)/(2*wjk.a*dd*dd);
                    eterm_deriv(ic1,s1,c1,j1,ic1,s1,c1,k1,M1f*M2f*pref_eeq*I02*part_jk2,o,v,dv_ak1,dv_ak1,dv_bk1,dv_bk1);

                    o2.set2(wk2,&I02);
                    cdouble dv_ak2= (zero_force ?  0.: (-djk2*ddv)*(v-t1)/(wjk2.a*dd*dd)) +t2/(wjk2.a*wjk2.a);
                    cVector_3 dv_bk2= zero_force ? cVector_3() :  ddv*(v-t1)/(2*wjk2.a*dd*dd);
                    eterm_deriv(ic2,s2,c2,j2,ic2,s2,c2,k2,M1f*M2f*pref_eeq*I0*part_jk,o2,v,dv_ak2,dv_ak2,dv_bk2,dv_bk2);
                  }

                }// k2
              } // j2

            } //c2
          }// s2
# endif
        } // k1
      }// j1
      ic1+=nspl[s1][c1]; // incrementing block1 wp address
    }// c1
  } // s1

  // transforming the forces to physical coordinates
  if(flag&(0x8|0x4) && !(flag&0x10))
    get_el_forces((flag&(~0x4))|0x8,fe_x,fe_p,fe_w,fe_pw,fe_c); // flag change: electronic forces were cleared already

  if(calc_ii)
    interaction_ii(flag,fi);

  return 1;
}


///\en Calcualtes the overlap between two electrons taking all split WPs into account.
///    Norms must be pre-calculated.
cdouble  AWPMD_split::overlap(int ic1, int s1, int c1,int ic2, int s2, int c2){
  cdouble sum(0,0);
  for(int j1=0;j1<nspl[s1][c1];j1++){
    double cj_re=split_c[s1][ic1+j1][0];
    double cj_im=split_c[s1][ic1+j1][1];
    cdouble ccj=cdouble(cj_re,-cj_im);
    WavePacket wj=wp[s1][ic1+j1];
    for(int k2=0;k2<nspl[s2][c2];k2++){
      double ck_re=split_c[s2][ic2+k2][0];
      double ck_im=split_c[s2][ic2+k2][1];
      cdouble ck=cdouble(ck_re,ck_im);

      WavePacket wk=wp[s2][ic2+k2];
      if(pbc)
        move_to_image(wj,wk);
      WavePacket wjk=conj(wj)*wk;
      sum+=ccj*ck*wjk.integral();
    }
  }
  return sum/sqrt(wf_norm[s1][c1]*wf_norm[s2][c2]);
}


/// adds the derivatives of Y in the term v*Y[s](c2,c1)
void AWPMD_split::y_deriv(cdouble v,int s,int c2, int c1){
  int ic=0;
  for(int c=0;c<ne[s];c++){
    for(int j=0;j<nspl[s][c];j++){
      E_der[s][ic]+=real((-M[s](c2,ic)*Y[s](c,c1)-Y[s](c2,c)*conj(M[s](c1,ic)))*v);
      ic++;
    }
  }
}

///\en Calculates interaction in the system of ni ions + electrons
/// the electonic subsystem must be previously setup by set_electrons, ionic by set_ions
/// 0x1   -- give back ion forces \n
/// 0x2   -- add ion forces to the existing set \n
/// 0x4   -- calculate electronic forces  \n
/// 0x8   -- add electronic forces to the existing arrays \n
/// 0x10  -- calculate internal electronic derivatives only: \n
///          will not update electronic force arrays, which may be null pointers, \n
///          the forces may be obtained then using \ref get_el_forces() for all WPs \n
///          or separately for each WP using \ref get_wp_force()
/// if PBCs are used the coords must be within a range [0, cell)
int AWPMD_split::interaction(int flag, Vector_3P fi, Vector_3P fe_x,
                             Vector_3P fe_p, double *fe_w, double *fe_pw, Vector_2P fe_c){
  //if(approx==HARTREE)
    //return interaction_hartree(flag,fi,fe_x,fe_p,fe_w,fe_pw,fe_c);

  //0. resize arrays if needed
  resize(flag);
  //1. clearing forces
  clear_forces(flag,fi,fe_x,fe_p,fe_w,fe_pw,fe_c);
  // calculate block norms and (optionally) derivatives
  calc_norms(flag);

  //2. calculating overlap matrix
  for(int s=0;s<2;s++){
    int nes = ne[s];
    if(nes == 0) continue;

    int ik=0;
    for(int k=0;k<nes;k++){

      Y[s].set(k,k,1.);  // Diagonal elements (=1)
      Oflg[s](k,k) = 1;
      int il=0;
      for(int l=0/*k+1*/;l<nes;il+=nspl[s][l],l++){ // incrementing block2 wp address
        if(l<k+1)
          continue;
        cdouble Okl=overlap(ik,s,k,il,s,l);
        Y[s].set(k,l,Okl);
        Oflg[s](k,l) = norm(Okl) > ovl_tolerance;
      }
      ik+=nspl[s][k]; // incrementing block1 wp address
    }
    O[s] = Y[s];  // save overlap matrix

    //3. inverting the overlap matrix
    int info=0;
    if(nes && approx!=HARTREE){
      /*FILE *f1=fopen(logfmt("matrO_%d.d",s),"wt");
      fileout(f1,Y[s],"%15g");
      fclose(f1);8*/

      ZPPTRF("L",&nes,Y[s].arr,&info);
      // analyze return code here
      if(info<0)
        return LOGERR(info,logfmt("AWPMD.interacton: call to ZPTRF failed (exitcode %d)!",info),LINFO);
      ZPPTRI("L",&nes,Y[s].arr,&info);
      if(info<0)
        return LOGERR(info,logfmt("AWPMD.interacton: call to ZPTRI failed (exitcode %d)!",info),LINFO);


      /*f1=fopen(logfmt("matrY_%d.d",s),"wt");
      fileout(f1,Y[s],"%15g");
      fclose(f1);*/
    }

    // Clearing matrices for electronic forces
    /*
    if(flag&0x4){
      Te[s].Set(0.);
      Tei[s].Set(0.);
    }*/
  }
# if 1
  // calculating the M matrix
  if(norm_needed || ( flag&(0x8|0x4) && approx!=HARTREE ) ){
    for(int s1=0;s1<2;s1++){
      if(approx!=HARTREE){
        M[s1].Set(0.);
        L[s1].Set(0.);
        if(norm_needed)
          Norm[s1].Set(0.);
      }
      else
        Lh[s1].assign(nvar[s1],0.);

      int indn=(nvar[s1]/10)*8;
      int ic1=0;
      for(int c1=0;c1<ne[s1];c1++){
        // L[s](c1,c1*)=0
        /*for(int j=0;j<nspl[s1][c1];j++){
          for(int i=0;i<10;i++)
            L[s1](c1,10*(ic1+j)+i)=0.;
        }*/

        int ic2=0;
        for(int c2=0;c2<ne[s1];ic2+=nspl[s1][c2],c2++){
          if(approx==HARTREE && c1!=c2) // only diagonals for Hartree
            continue;
          if(c2<c1) // taken assymmetry into account
            continue;
          if( !Oflg[s1](c1,c2) )
            continue; // non-overlapping WPs


          if(approx==HARTREE && norm_needed) // initializing the matrices with zero
            Normh[s1][c1].Set(0);


          double sq_norm12=sqrt(wf_norm[s1][c1]*wf_norm[s1][c2]);
          double pref=1./(sq_norm12);


          // WP blocks:
          for(int j1=0;j1<nspl[s1][c1];j1++){
            cdouble cj1(split_c[s1][ic1+j1][0],split_c[s1][ic1+j1][1]);
            WavePacket wj1=wp[s1][ic1+j1];

            OverlapDeriv o12;
            o12.set1(wj1);

           for(int k2=(c1==c2 ? j1: 0); k2<nspl[s1][c2];k2++){
              double M12= (c1==c2 && j1==k2) ? 0.5 : 1;

              cdouble ck2(split_c[s1][ic2+k2][0],split_c[s1][ic2+k2][1]);

              // electrons kinetic energy
              WavePacket wk2=wp[s1][ic2+k2];
              if(pbc)
                move_to_image(wj1,wk2);

              WavePacket wjk12=conj(wj1)*wk2;
              cdouble I012 = wjk12.integral();
              cdouble part_jk12=conj(cj1)*ck2;

              o12.set2(wk2,&I012);


              cdouble der_k[10], der_j[10];
              // over a_k_re
              der_k[0]=   part_jk12*o12.da2_re();
              // over a_k_im
              der_k[1]=   part_jk12*o12.da2_im();
              // over a_j_re
              der_j[0]=   part_jk12*o12.da1_re();
              // over a_j_im
              der_j[1]=   part_jk12*o12.da1_im();

              for(int i=0;i<3;i++){
                // over b_k_re
                der_k[2+2*i]=   part_jk12*o12.db2_re(i);
                // over b_k_im
                der_k[2+2*i+1]= part_jk12*o12.db2_im(i);
                // over b_j_re
                der_j[2+2*i]=   part_jk12*o12.db1_re(i);
                // over b_j_im
                der_j[2+2*i+1]= part_jk12*o12.db1_im(i);
              }

              // over ck_re
              der_k[8]=conj(cj1)*I012;
              // over ck_im
              der_k[9]=i_unit*conj(cj1)*I012;
              // over cj_re
              der_j[8]=ck2*I012;
              // over cj_im
              der_j[9]=-i_unit*ck2*I012;

              if(j1==0){ // add all together instead of adding by parts
                cdouble t=-O[s1](c1,c2)/pref;   //-part_jk12*I012;
                for(int i=0;i<8;i++){
                  der_j[i]+=t*wf_norm_der[s1][8*(ic1+j1)+i];
                  der_k[i]+=t*wf_norm_der[s1][8*(ic2+k2)+i];
                }
                der_j[8]+=t*wf_norm_der[s1][indn+2*(c1+j1)];
                der_j[9]+=t*wf_norm_der[s1][indn+2*(c1+j1)+1];
                der_k[8]+=t*wf_norm_der[s1][indn+2*(c2+k2)];
                der_k[9]+=t*wf_norm_der[s1][indn+2*(c2+k2)+1];
              }

              if(approx!=HARTREE){
                for(int i=0;i<10;i++){
                  L[s1](c1,10*(ic2+k2)+i)+=M12*pref*der_k[i];
                  L[s1](c2,10*(ic1+j1)+i)+=M12*pref*conj(der_j[i]);
                }

                // M=Y*L
                for(int g=0;g<ne[s1];g++){
                  for(int i=0;i<10;i++){
                    M[s1](g,10*(ic2+k2)+i)+=M12*pref*Y[s1](g,c1)*der_k[i];
                    M[s1](g,10*(ic1+j1)+i)+=M12*pref*Y[s1](g,c2)*conj(der_j[i]);
                  }
                }//g
              }
              else{ // HARTREE
                for(int i=0;i<10;i++){
                  Lh[s1][10*(ic2+k2)+i]+=M12*pref*der_k[i];
                  Lh[s1][10*(ic1+j1)+i]+=M12*pref*der_k[i];
                }
              }
              if(norm_needed){ // filling part of norm matrix
                o12.calc_der_overlap(false,conj(cj1),ck2);
                if(approx!=HARTREE){
                  for(int i=0;i<10;i++)  // 10x10
                    for(int j=0;j<10;j++)
                      Norm[s1](10*(ic1+j1)+i,10*(ic2+k2)+j)=-2*imag(o12.IDD((int)i,(int)j))/one_h;
                }
                else{
                  for(int i=0;i<10;i++)  // 10x10
                    for(int j=0;j<10;j++)
                      Normh[s1][c1](10*j1+i,10*k2+j)=-2*imag(o12.IDD((int)i,(int)j))/one_h;
                }
              }
            } // k2
          } //j1
        }// c2
        ic1+=nspl[s1][c1]; // incrementing block1 wp address
      }// c1
      if(norm_needed){  // filling the rest of norm_matrix
        if(approx!=HARTREE){
          int ic1=0;
          for(int c1=0;c1<ne[s1];c1++){
            int ic2=0;
            for(int c2=0;c2<ne[s1];ic2+=nspl[s1][c2],c2++){
              if(c2<c1)
                continue;
              for(int j1=0;j1<10*nspl[s1][c1];j1++){
                for(int k2=(c1==c2 ? j1: 0); k2<10*nspl[s1][c2];k2++){
                  cdouble y=0.;
                  for(int g=0;g<ne[s1];g++)
                    y+=conj(L[s1](g,10*ic1+j1))*M[s1](g,10*ic2+k2);
                  cdouble elm=-conj(L[s1](c2,10*ic1+j1))*wf_norm_der[s1][10*ic2+k2]-L[s1](c1,10*ic2+k2)*wf_norm_der[s1][10*ic1+j1]-O[s1](c1,c2)*wf_norm_der[s1][10*ic2+k2]*wf_norm_der[s1][10*ic1+j1];
                  elm-=y;
                  elm*=Y[s1](c2,c1);
                  Norm[s1](10*ic1+j1,10*ic2+k2)+=-2*imag(elm)/one_h;
                }// k2
              } // j1
            }// c2
            ic1+=nspl[s1][c1]; // incrementing block1 wp address
          }
        }
        else{ // HARTREE
          int ic1=0;
          for(int c1=0;c1<ne[s1];c1++){
            for(int j1=0;j1<10*nspl[s1][c1];j1++){
              for(int k2=j1; k2<10*nspl[s1][c1];k2++){
                cdouble y=conj(Lh[s1][10*ic1+j1])*Lh[s1][10*ic1+k2];
                cdouble elm=-conj(Lh[s1][10*ic1+j1])*wf_norm_der[s1][10*ic1+k2]-Lh[s1][10*ic1+k2]*wf_norm_der[s1][10*ic1+j1]-wf_norm_der[s1][10*ic1+k2]*wf_norm_der[s1][10*ic1+j1];
                elm-=y;
                Normh[s1][c1](j1,k2)+=-2*imag(elm)/one_h;
              }
            }
            ic1+=nspl[s1][c1]; // incrementing block1 wp address
          }
        }
      } // norm
    } // s1
  } // flag
# endif
  Vector_3 ndr;
  //  calculating single particle contribution
  Edc = Edk = Eee = Ew = 0.;
  for(int s1=0;s1<2;s1++){
    Ee[s1]=Eei[s1]=0.;
    int ic1=0;
    for(int c1=0;c1<ne[s1];c1++){

      //double Ee1=0., Ew1=0., Eei1=0.;

      int ic2=0;
      for(int c2=0;c2<ne[s1];ic2+=nspl[s1][c2],c2++){
        if( !Oflg[s1](c1,c2) )
            continue; // non-overlapping WPs

        if(c2<c1) // taken as M factor into account
          continue;

        double sq_norm12=sqrt(wf_norm[s1][c1]*wf_norm[s1][c2]);
        double pref=-h2_me/(2*sq_norm12); // ekin
        double pref_ei=coul_pref/sq_norm12;

        cdouble yy;
        if(approx==HARTREE)
          yy=1.;
        else
          yy=Y[s1](c2,c1);
        cdouble Tc1c2=0.;

        if((approx!=HARTREE || c2==c1) && Oflg[s1](c1,c2)){ // only diagonal terms for Hartree, and overlap nonzero
          // WP blocks:
          for(int j1=0;j1<nspl[s1][c1];j1++){
            cdouble cj1(split_c[s1][ic1+j1][0],split_c[s1][ic1+j1][1]);
            WavePacket wj1=wp[s1][ic1+j1];

            OverlapDeriv o12;
            if(flag&(0x8|0x4)) //electron forces needed
              o12.set1(wj1);

            for(int k2=(c1==c2 ? j1: 0); k2<nspl[s1][c2];k2++){
              int M12=(c1==c2 && j1==k2 ? 1: 2);
              double M12pe, M12pf;
              _mytie(M12pe,M12pf)=check_part1(s1,ic1+j1,ic2+k2)*M12;



              cdouble ck2(split_c[s1][ic2+k2][0],split_c[s1][ic2+k2][1]);


              // electrons kinetic energy
              WavePacket wk2=wp[s1][ic2+k2];
              if(pbc)
                move_to_image(wj1,wk2);

              WavePacket wjk12=conj(wj1)*wk2;
              cdouble I012 = wjk12.integral();
              cdouble part_jk12=conj(cj1)*ck2;
              cVector_3 djk12=wjk12.b/(2.*wjk12.a);

              // kinetic energy contribution
              if(M12pf){
                cVector_3 v1=conj(wj1.b)*wk2.a-wk2.b*conj(wj1.a);
                cdouble v=(v1*v1)/wjk12.a;
                v-=6.*wk2.a*conj(wj1.a);
                v/=wjk12.a;
                //v=1.;

                double dE=M12pe*real(part_jk12*I012*v*yy)*pref;
                //double dE=real(yy);
                //Tc1c2+=1.;
                Ee[s1] += dE;
                Eep[s1][ic1+j1]+= 0.5*dE; //per particle energy
                Eep[s1][ic2+k2]+= 0.5*dE;

                // diagonal term
                if(M12==1)
                  Edk+=dE;

                if(flag&(0x8|0x4)){ //electron forces needed
                  cVector_3 tv=wk2.b*conj(wj1.a)-conj(wj1.b)*wk2.a;
                  cdouble ajk2=wjk12.a*wjk12.a;
                  cdouble ajk3=ajk2*wjk12.a;
                  cdouble dv_aj_conj=     -2*wk2.a*(3*wjk12.a*wk2.a-tv*wjk12.b)/ajk3;
                  cdouble dv_ak=-2*conj(wj1.a)*(3*wjk12.a*conj(wj1.a)+tv*wjk12.b)/ajk3;
                  cVector_3 dv_bj_conj=-2*tv*wk2.a/ajk2;
                  cVector_3 dv_bk=2*tv*conj(wj1.a)/ajk2;


                  o12.set2(wk2,&I012);
                  // calculate full derivative of the term pref*conj(cj)*ck*Ijk*vjk/sqrt(nrm(s1)*nrm(s2))
                  // denominator must be included in pref
                  eterm_deriv(ic1,s1,c1,j1,ic2,s1,c2,k2,M12pf*pref*yy,o12,v,dv_aj_conj,dv_ak,dv_bj_conj,dv_bk);

                  Tc1c2+=M12pf*part_jk12*I012*v*pref; // all without Y
                }
              } // M12pe
              // e-i energy
              cdouble sum=0.;
              for(int i=0;i<ni;i++){  // ions loop
                double M12pie, M12pif;
                _mytie(M12pie,M12pif)=check_part1ei(s1,ic1+j1,ic2+k2,i)*M12;
                if(!M12pif)
                  continue;
                cVector_3 gjki=djk12 - cVector_3(xi[i]);

                if(pbc) // correcting the real part (distance) according to PBC
                  gjki=rcell1(gjki,cell,pbc);
                  //-Igor- gkli=cVector_3(real(gkli).rcell1(cell,pbc),imag(gkli));

                cdouble ngjki=gjki.norm();
                cdouble sqajk=sqrt(wjk12.a);
                //cdouble ttt = cerf_div(ngjki,c);
                cdouble v=cerf_div(ngjki,sqajk);
                double pref_eiq=pref_ei*qi[i]*(qe[s1][ic1+j1]+qe[s1][ic2+k2])*0.5;
                cdouble dE=pref_eiq*I012*v*part_jk12;
                sum+=M12pie*dE; // sum without Y
                dE*=yy;

                Eeip[s1][ic1+j1]+= 0.25*real(dE); //per particle energy
                Eeip[s1][ic2+k2]+= 0.25*real(dE);
                Eiep[i]+=0.5*real(dE);

                cdouble t1, t2;
                bool zero_force=(fabs(real(ngjki))+fabs(imag(ngjki))<1e-10);
                if(flag&(0x8|0x4|0x3) ){ //some derivatives needed
                  cdouble arg=ngjki*sqajk;
                  t1=two_over_sqr_pi*exp(-arg*arg);
                  t2=t1/(2*sqajk);
                  t1*=sqajk;
                }

                if(flag&0x3 && !zero_force){// calculate forces on ions
                  cdouble dEw=(pref_eiq*yy)*I012*t1*part_jk12;
                  dEw=(dE-dEw)/ngjki;
                  Vector_3 dir=-real(gjki);
                  dir.normalize();
                  fi[i]+=M12pif*real( dEw)*dir;
                }
                if(flag&(0x8|0x4)){ //electron forces needed
                  cdouble dv_ak= (zero_force ?  0.: (djk12*gjki)*(v-t1)/(wjk12.a*ngjki*ngjki)) +t2;
                  cVector_3 dv_bk= zero_force ? cVector_3() :  -gjki*(v-t1)/(2*wjk12.a*ngjki*ngjki);
                  eterm_deriv(ic1,s1,c1,j1,ic2,s1,c2,k2,M12pif*pref_eiq*yy,o12,v,dv_ak,dv_ak,dv_bk,dv_bk);
                }
              }
              // ions energy contribution
              Eei[s1] += real(sum*yy);

              // diagonal term
              if(M12==1)
                Edc+= real(sum*yy);;

              if(flag&(0x8|0x4)) //electron forces needed
                Tc1c2+=sum;

              // extra constraint energy, will only arrive here when M12p=1
              if((c1==c2 && j1==k2) && constraint == HARM){
                cdouble v=conj(wj1.a)*wk2.a;
                double dE= M12pe * harm_w0_4 * real(part_jk12*v)/wf_norm[s1][c1];
                Ew += dE;

                Ewp[s1][ic1+j1]+= dE; //per particle energy


                if(flag&(0x8|0x4)){ //electron forces needed
                  cdouble dv_ak=     conj(wj1.a);
                  cdouble dv_aj_conj =wk2.a;
                  eterm_deriv(ic1,s1,c1,j1,ic1,s1,c1,k2,M12pf*harm_w0_4/wf_norm[s1][c1],o12,v,dv_aj_conj,dv_ak,cVector_3(),cVector_3());
                }
              }
            }// k2
          }// j1
          if(flag&(0x8|0x4) && approx!=HARTREE){ //electron forces needed
            // adding Y derivative term for all variables (te and tei terms)
            y_deriv(Tc1c2,s1,c2,c1);
          }
        } // !HARTREE or spins or overlap



# if 1 // pair by pair sum
        // second block
        // e-e interaction
        for(int s2=s1;s2<2;s2++){
          if(approx==HARTREE && c1!=c2) // only Vkmkm terms for Hartree
               continue;
          if(/*s1==s2 &&*/ c2<c1) // pair selection term
             continue;

          int ic3=0; // starting index of the wp for current electron
          for(int c3=0 ;c3<ne[s2];ic3+=nspl[s2][c3],c3++){ // incrementing block2 wp address
            if(approx==HARTREE && s1==s2 && c3==c1) // [Vkkmn contribution for same spin is 0], also  no Vkkkk terms for Hartree (=)
              continue;
            if(s1==s2 && c3<=c1) // and general Coulomb sum: i<j (<) //??
              continue;

            int ic4=0;
            for(int c4=0; c4<ne[s2];ic4+=nspl[s2][c4],c4++){
              if(approx==HARTREE && c4!=c3) // only Vkmkm terms for Hartree
                continue;
              //if(s1==s2 && c4==c2) // Vklmm contribution for same spin is 0 for antisymmetrized approximations, also Vmmmm is 0 for Hartree
              //  continue;
              if(s1==s2 && c4<=c2) // pair selection term: Vklmm
                continue;
              if(/*s1==s2 &&*/ c4<c3) // pair selection term
                continue;

              double sq_norm34=sqrt(wf_norm[s2][c3]*wf_norm[s2][c4]);
              double pref_ee=coul_pref/(sq_norm12*sq_norm34);


              cdouble yy;
              double K=1.;
              if(approx==HARTREE){
                yy=1.;
              }
              else{
                if(s1==s2){ // same spin antisymmetrized term
                  yy=Y[s1](c2,c1)*Y[s1](c4,c3)-Y[s1](c2,c3)*Y[s1](c4,c1);  // check the order of c
                }
                else{ // different spin atisymmetrized term
                  //if(approx==UHF)
                    //K=2.;
                  yy=K*Y[s2](c4,c3)*Y[s1](c2,c1);
                }
                //yy=1.;
                //yy=Y[s1](c2,c1)*Y[s1](c4,c3);
                //yy=Y[s1](c2,c3)*Y[s1](c4,c1);
                //yy=Y[s1](c2,c1)*Y[s1](c4,c3)-Y[s1](c2,c3)*Y[s1](c4,c1);
              }

              cdouble Vy=0.;
              // WP blocks: Vc1c3c2c4
              for(int j1=0;j1<nspl[s1][c1];j1++){
                cdouble cj1(split_c[s1][ic1+j1][0],split_c[s1][ic1+j1][1]);
                WavePacket wj1=wp[s1][ic1+j1];

                OverlapDeriv o12,o14;
                if(flag&(0x8|0x4)){ //electron forces needed
                  o12.set1(wj1);
                  o14.set1(wj1);
                }

                for(int k2=(approx==HARTREE ? j1: 0); k2<nspl[s1][c2];k2++){
                  int M12=(c1==c2 && j1==k2 ? 1: 2);


                  cdouble ck2(split_c[s1][ic2+k2][0],split_c[s1][ic2+k2][1]);


                  WavePacket wk2=wp[s1][ic2+k2];
                  if(pbc)
                    move_to_image(wj1,wk2);

                  WavePacket wjk12=conj(wj1)*wk2;
                  cdouble I012 = wjk12.integral();
                  cdouble part_jk12=conj(cj1)*ck2;
                  cVector_3 djk12=wjk12.b/(2.*wjk12.a);

                  if(flag&(0x8|0x4)) //electron forces needed
                    o12.set2(wk2,&I012);


                  for(int j3=0;j3<nspl[s2][c3];j3++){

                    cdouble cj3(split_c[s2][ic3+j3][0],split_c[s2][ic3+j3][1]);
                    WavePacket& wj3=wp[s2][ic3+j3];

                    OverlapDeriv o34, o32;
                    if(flag&(0x8|0x4)){ //electron forces needed
                      o34.set1(wj3);
                      o32.set1(wj3);
                      //o32.set2(wk2);
                    }

                    // 3-2
                    WavePacket wjk32;
                    cdouble I032=0., part_jk32;
                    cVector_3 djk32;
                    if(s1==s2 || approx==UHF){
                      WavePacket wk2=wp[s2][ic2+k2];
                      if(pbc)
                         move_to_image(wj3,wk2);
                      wjk32=conj(wj3)*wk2;
                      I032 = wjk32.integral();
                      part_jk32=conj(cj3)*ck2;
                      djk32=wjk32.b/(2*wjk32.a);
                      if(flag&(0x8|0x4)) //electron forces needed
                        o32.set2(wk2,&I032);
                    }


                    for(int k4=(approx==HARTREE ? j3: 0);k4<nspl[s2][c4];k4++){
                      int M34=(c3==c4 && j3==k4 ? 1: 2);



                      double M0;
                      if(approx==HARTREE)
                        M0=M12*M34;
                      else{
                        M0= (c1==c2 && c3==c4 ? 1: 2); // will have exchange term for different pairs instead of M12*M34 factor
                      }
                      double Me, Mf;
                      _mytie(Me,Mf)=check_part1(s1,ic1+j1,ic2+k2)*M0;
                      if(!Mf)
                        continue;

                      cdouble ck4(split_c[s2][ic4+k4][0],split_c[s2][ic4+k4][1]);

                      // 3-4
                      WavePacket wk4=wp[s2][ic4+k4];
# if 1
                      if(pbc)
                        move_to_image(wj3,wk4);
                      WavePacket wjk34=conj(wj3)*wk4;
                      cdouble I034 = wjk34.integral();
                      cdouble part_jk34=conj(cj3)*ck4;
                      cVector_3 djk34=wjk34.b/(2*wjk34.a);

                      cVector_3 ddv=djk12-djk34;
                      cdouble dd=ddv.norm();
                      cdouble aa=1./sqrt(1./wjk12.a+1./wjk34.a);
                      cdouble v=cerf_div(dd,aa);
                      double pref_eeq=pref_ee*(qe[s1][ic1+j1]+qe[s1][ic2+k2])*(qe[s2][ic3+j3]+qe[s2][ic4+k4])*0.25;
                      cdouble Vj1j3k2k4=pref_eeq*I012*I034*v*part_jk12*part_jk34;  // Vklmn
                      double dE=Me*real(Vj1j3k2k4*yy);
                      Eee+=dE;

                      Eeep[s1][ic1+j1]+=0.25*dE; // per-particle energy
                      Eeep[s1][ic2+k2]+=0.25*dE;
                      Eeep[s2][ic3+j3]+=0.25*dE;
                      Eeep[s2][ic4+k4]+=0.25*dE;

                      // diagonal term
                      if(M12==1 && M34 ==1)
                        Edc+=dE;

                      bool zero_force=(fabs(real(dd))+fabs(imag(dd))<1e-10);
                      if(flag&(0x8|0x4) ){ //electron forces needed
                        cdouble arg=dd*aa;
                        cdouble t1=two_over_sqr_pi*exp(-arg*arg)*aa;
                        cdouble t2=t1*aa*aa/2;


                        cdouble dv_ak1= (zero_force ?  0.: (djk12*ddv)*(v-t1)/(wjk12.a*dd*dd)) +t2/(wjk12.a*wjk12.a);
                        cVector_3 dv_bk1= zero_force ? cVector_3() :  -ddv*(v-t1)/(2*wjk12.a*dd*dd);
                        eterm_deriv(ic1,s1,c1,j1,ic2,s1,c2,k2,Mf*pref_eeq*I034*part_jk34*yy,o12,v,dv_ak1,dv_ak1,dv_bk1,dv_bk1);

                        o34.set2(wk4,&I034);
                        cdouble dv_ak2= (zero_force ?  0.: (-djk34*ddv)*(v-t1)/(wjk34.a*dd*dd)) +t2/(wjk34.a*wjk34.a);
                        cVector_3 dv_bk2= zero_force ? cVector_3() :  ddv*(v-t1)/(2*wjk34.a*dd*dd);
                        eterm_deriv(ic3,s2,c3,j3,ic4,s2,c4,k4,Mf*pref_eeq*I012*part_jk12*yy,o34,v,dv_ak2,dv_ak2,dv_bk2,dv_bk2);

                        Vy+=Me*Vj1j3k2k4; // all without yy
                      }// force
# endif
# if 1
                      // 1-4
                      if(approx!=HARTREE && (approx==UHF || s1==s2)){
                        wk4=wp[s2][ic4+k4];
                        if(pbc)
                          move_to_image(wj1,wk4);
                        WavePacket wjk14=conj(wj1)*wk4;
                        cdouble I014 = wjk14.integral();
                        cdouble part_jk14=conj(cj1)*ck4;
                        cVector_3 djk14=wjk14.b/(2*wjk14.a);

                        cVector_3 ddv=djk32-djk14;
                        cdouble dd=ddv.norm();
                        cdouble aa=1./sqrt(1./wjk32.a+1./wjk14.a);
                        cdouble v=cerf_div(dd,aa);


                        double pref_eeq=pref_ee*(qe[s1][ic1+j1]+qe[s2][ic4+k4])*(qe[s2][ic3+j3]+qe[s1][ic2+k2])*0.25;
                        cdouble Vj1j3k4k2=pref_eeq*I032*I014*v*part_jk32*part_jk14;  // Vklmn
                        double dE=-Me*real(Vj1j3k4k2*yy);
                        Eee+=dE;

                        Eeep[s1][ic1+j1]+=0.25*dE; // per-particle energy
                        Eeep[s1][ic2+k2]+=0.25*dE;
                        Eeep[s2][ic3+j3]+=0.25*dE;
                        Eeep[s2][ic4+k4]+=0.25*dE;


                        bool zero_force=(fabs(real(dd))+fabs(imag(dd))<1e-10);
                        if(flag&(0x8|0x4) ){ //electron forces needed
                          cdouble arg=dd*aa;
                          cdouble t1=two_over_sqr_pi*exp(-arg*arg)*aa;
                          cdouble t2=t1*aa*aa/2;

                          o14.set2(wk4,&I014);
                          cdouble dv_ak1= (zero_force ?  0.: (-djk14*ddv)*(v-t1)/(wjk14.a*dd*dd)) +t2/(wjk14.a*wjk14.a);
                          cVector_3 dv_bk1= zero_force ? cVector_3() :  ddv*(v-t1)/(2*wjk14.a*dd*dd);
                          eterm_deriv(ic1,s1,c1,j1,ic4,s2,c4,k4,-Mf*pref_eeq*I032*part_jk32*yy,o14,v,dv_ak1,dv_ak1,dv_bk1,dv_bk1);


                          cdouble dv_ak2= (zero_force ?  0.: (djk32*ddv)*(v-t1)/(wjk32.a*dd*dd)) +t2/(wjk32.a*wjk32.a);
                          cVector_3 dv_bk2= zero_force ? cVector_3() :  -ddv*(v-t1)/(2*wjk32.a*dd*dd);
                          eterm_deriv(ic3,s2,c3,j3,ic2,s1,c2,k2,-Mf*pref_eeq*I014*part_jk14*yy,o32,v,dv_ak2,dv_ak2,dv_bk2,dv_bk2);

                          Vy+=-Me*Vj1j3k4k2; // all without yy
                        }// force
                      } // s1==s2
# endif
                    }// k4
                  }// j3
                } // k2
              } // j1

              // derivative of yy term
              if(approx!=HARTREE && flag&(0x8|0x4)){ //electron forces needed
                //yy=Y[s1](c2,c1)*Y[s1](c4,c3)-Y[s1](c2,c3)*Y[s1](c4,c1);
# if 1
                if(s1==s2){
                  // yy
                  y_deriv(Vy*Y[s1](c4,c3),s1,c2,c1);
                  y_deriv(Vy*Y[s1](c2,c1),s1,c4,c3);
                  y_deriv(-Vy*Y[s1](c4,c1),s1,c2,c3);
                  y_deriv(-Vy*Y[s1](c2,c3),s1,c4,c1);
                }
                else{
                  y_deriv(K*Vy*Y[s2](c4,c3),s1,c2,c1);
                  y_deriv(K*Vy*Y[s1](c2,c1),s2,c4,c3);
                }
# endif
               // y_deriv(Vy*Y[s1](c4,c1),s1,c2,c3);
               // y_deriv(Vy*Y[s1](c2,c3),s1,c4,c1);
              }
            } // c4
          } // c3
        } // s2
# endif
      }// c2
      ic1+=nspl[s1][c1]; // incrementing block1 wp address
    }// c1
  } // s1

  // transforming the forces to physical coordinates
  if(flag&(0x8|0x4) && !(flag&0x10))
    get_el_forces((flag&(~0x4))|0x8,fe_x,fe_p,fe_w,fe_pw,fe_c); // flag change: electronic forces were cleared already

  if(calc_ii)
    interaction_ii(flag,fi);

  return 1;
}
