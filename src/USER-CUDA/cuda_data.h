/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef _CUDA_DATA_H_
#define _CUDA_DATA_H_


enum copy_mode {x, xx, xy, yx, xyz, xzy}; // yxz, yzx, zxy, zyx not yet implemented since they were not needed yet
//xx==x in atom_vec x is a member therefore copymode x produces compile errors
#include "cuda_shared.h"
#include "cuda_wrapper_cu.h"
#include "cuda_data_cu.h"
#include <ctime>

#include <cstdio>
#include <typeinfo>
template <typename host_type, typename dev_type, copy_mode mode>
class cCudaData
{
        protected:
        void** buffer;
        int* buf_size;
        host_type* host_data;
        dev_array* dev_data_array;
        dev_type* temp_data;
        unsigned nbytes;
        bool owns_dev_array;
        bool current_data_on_device; //this is not yet working as intended and therefore deactivated
        bool current_data_on_host;
        bool is_continues;
        bool pinned;

        public:
        cCudaData(host_type* host_data, dev_array* dev_data_array, unsigned dim_x, unsigned dim_y=0, unsigned dim_z=0, bool is_pinned=false);
        cCudaData(host_type* host_data, unsigned dim_x, unsigned dim_y=0, unsigned dim_z=0, bool is_pinned=false);
        ~cCudaData();
        void* dev_data() {if(dev_data_array!=NULL) return dev_data_array->dev_data; else return NULL;};
        void set_dev_data(void* adev_data) {dev_data_array->dev_data=adev_data;};
        void set_dev_array(dev_array* adev_array) {dev_data_array=adev_array;};
        void set_host_data(host_type* host_data);
        void* get_host_data() { return host_data;};
        void set_buffer(void** buffer,int* buf_size,bool ais_continues);
        unsigned int* get_dim() {return dev_data_array->dim;};
        // if you want to upload data to the gpu, which will not change there, then set will_be_changed=false
        // if you want to upload data to the gpu and update it there, then set will_be_changed=true (default)
        void upload(bool will_be_changed=true);
        void uploadAsync(int stream, bool will_be_changed=true );
        // if you want to download data just to have a look at it, then set will_be_changed=false
        // if you are going to modify the downloaded data, then set will_be_changed=true (default)
        void download(bool will_be_changed=true);
        void downloadAsync(int stream);
        void memset_device(int value);
        void device_data_has_changed() {current_data_on_device=false;}
        void host_data_has_changed() {current_data_on_host=false;}
        int dev_size() {
                int size = dev_data_array->dim[0]*sizeof(dev_type);
                if(dev_data_array->dim[1]) size*=dev_data_array->dim[1];
                if(dev_data_array->dim[2]) size*=dev_data_array->dim[2];
                return size;}
};


template <typename host_type, typename dev_type, copy_mode mode>
cCudaData<host_type, dev_type, mode>
::cCudaData(host_type* host_data, dev_array* dev_data_array, unsigned dim_x, unsigned dim_y, unsigned dim_z, bool is_pinned)
{
        pinned=is_pinned;
        owns_dev_array = false;
        current_data_on_device = false;
        current_data_on_host = false;
        is_continues = false;
        this->host_data = host_data;
        this->dev_data_array = dev_data_array;
        unsigned ndev;
        if((mode == x)||(mode==xx))
        {
                ndev = dim_x;
                dev_data_array->dim[0] = dim_x;
                dev_data_array->dim[1] = 0;
                dev_data_array->dim[2] = 0;
        }
        else if(mode == xy || mode == yx )
        {
                ndev = dim_x * dim_y;
                dev_data_array->dim[0] = dim_x;
                dev_data_array->dim[1] = dim_y;
                dev_data_array->dim[2] = 0;
        }
        else
        {
                ndev = dim_x * dim_y * dim_z;
                dev_data_array->dim[0] = dim_x;
                dev_data_array->dim[1] = dim_y;
                dev_data_array->dim[2] = dim_z;
        }
        nbytes = ndev * sizeof(dev_type);
        if(nbytes<=0)
        {
                host_data=NULL;
                temp_data=NULL;
                dev_data_array->dev_data=NULL;
                return;
        }

        dev_data_array->dev_data = CudaWrapper_AllocCudaData(nbytes);
        if(((mode!=x)&&(mode!=xx)) || typeid(host_type) != typeid(dev_type))
        {
                if(not pinned)
                temp_data = new dev_type[ndev];
                else
                {
                        temp_data = (dev_type*) CudaWrapper_AllocPinnedHostData(ndev*sizeof(dev_type));
                }
        }
}

template <typename host_type, typename dev_type, copy_mode mode>
cCudaData<host_type, dev_type, mode>
::cCudaData(host_type* host_data, unsigned dim_x, unsigned dim_y, unsigned dim_z, bool is_pinned)
{
        pinned=is_pinned;
        this->dev_data_array = new dev_array;
        this->owns_dev_array = true;
        current_data_on_device = false;
        current_data_on_host = false;
        is_continues = false;
        this->host_data = host_data;
        unsigned ndev;
        if((mode == x)||(mode==xx))
        {
                ndev = dim_x;
                dev_data_array->dim[0] = dim_x;
                dev_data_array->dim[1] = 0;
                dev_data_array->dim[2] = 0;
        }
        else if(mode == xy || mode == yx )
        {
                ndev = dim_x * dim_y;
                dev_data_array->dim[0] = dim_x;
                dev_data_array->dim[1] = dim_y;
                dev_data_array->dim[2] = 0;
        }
        else
        {
                ndev = dim_x * dim_y * dim_z;
                dev_data_array->dim[0] = dim_x;
                dev_data_array->dim[1] = dim_y;
                dev_data_array->dim[2] = dim_z;
        }
        nbytes = ndev * sizeof(dev_type);
        if(nbytes<=0)
        {
                host_data=NULL;
                temp_data=NULL;
                dev_data_array->dev_data=NULL;
                return;
        }

        dev_data_array->dev_data = CudaWrapper_AllocCudaData(nbytes);
        if(((mode!=x)&&(mode!=xx)) || (typeid(host_type) != typeid(dev_type)))
        {
                if(not pinned)
                temp_data = new dev_type[ndev];
                else
                {
                        temp_data = (dev_type*) CudaWrapper_AllocPinnedHostData(ndev*sizeof(dev_type));
                }
        }
}

template <typename host_type, typename dev_type, copy_mode mode>
cCudaData<host_type, dev_type, mode>
::~cCudaData()
{
        if(((mode!=x)&&(mode!=xx)) || typeid(host_type) != typeid(dev_type))
        {
                if(not pinned)
                delete [] temp_data;
                else
                {
                        CudaWrapper_FreePinnedHostData((void*)temp_data);
                }
        }
        if((dev_data_array->dev_data)&&(nbytes>0))
        CudaWrapper_FreeCudaData(dev_data_array->dev_data,nbytes);
        if(owns_dev_array) delete dev_data_array;
}

template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::set_host_data(host_type* host_data)
{
        this->host_data = host_data;
}

template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::upload(bool will_be_changed)
{
        // if current data is already up, do not re-upload it
//        if(current_data_on_device) return;
    if(buffer&&is_continues)
    {
           printf("Actual Buffer: %p %i\n",*buffer,*buf_size);
            if(typeid(host_type)==typeid(double))
            {
              if(typeid(dev_type)==typeid(double))
              {
                      CudaData_Upload_DoubleDouble((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
              else if(typeid(dev_type)==typeid(float))
              {
                      CudaData_Upload_DoubleFloat((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
            }
            else if(typeid(host_type)==typeid(float))
            {
              if(typeid(dev_type)==typeid(double))
              {
                      CudaData_Upload_FloatDouble((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
              else if(typeid(dev_type)==typeid(float))
              {
                      CudaData_Upload_FloatFloat((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
            }
            else if(typeid(host_type)==typeid(int))
            {
              if(typeid(dev_type)==typeid(int))
              {
                      CudaData_Upload_IntInt((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
            }
    }
        switch(mode)
        {
                case x:
                {
                        if(typeid(host_type) == typeid(dev_type))
                                CudaWrapper_UploadCudaData(host_data, dev_data_array->dev_data, nbytes);
                        else
                        {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                          for(unsigned i=0; i<dev_data_array->dim[0]; ++i) temp_data[i] = static_cast<dev_type>(host_data[i]);
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                          CudaWrapper_UploadCudaData(temp_data, dev_data_array->dev_data, nbytes);
                        }
                        break;
                }

                case xx:
                {
                        if(typeid(host_type) == typeid(dev_type))
                                CudaWrapper_UploadCudaData(host_data, dev_data_array->dev_data, nbytes);
                        else
                        {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i) temp_data[i] = static_cast<dev_type>(host_data[i]);
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                                CudaWrapper_UploadCudaData(temp_data, dev_data_array->dev_data, nbytes);
                        }
                        break;
                }

                case xy:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        {
                                dev_type* temp = &temp_data[i * dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        temp[j] = static_cast<dev_type>((reinterpret_cast<host_type**>(host_data))[i][j]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaData(temp_data, dev_data_array->dev_data, nbytes);
                        break;
                }

                case yx:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[j*dev_data_array->dim[0]];
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                                {
                                        temp[i] = static_cast<dev_type>(reinterpret_cast<host_type**>(host_data)[i][j]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaData(temp_data, dev_data_array->dev_data, nbytes);
                        break;
                }
                case xyz:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[(i*dev_data_array->dim[1]+j)*dev_data_array->dim[2]];
                                for(unsigned k=0; k<dev_data_array->dim[2]; ++k)
                                {
                                        temp[k] = static_cast<dev_type>(reinterpret_cast<host_type***>(host_data)[i][j][k]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaData(temp_data, dev_data_array->dev_data, nbytes);
                        break;
                }

                case xzy:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        for(unsigned k=0; k<dev_data_array->dim[2]; ++k)
                        {
                                dev_type* temp = &temp_data[(i*dev_data_array->dim[2]+k)*dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        temp[j] = static_cast<dev_type>(reinterpret_cast<host_type***>(host_data)[i][j][k]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaData(temp_data, dev_data_array->dev_data, nbytes);
                        break;
                }
        }
        // we have uploaded the data to the device, i.e.:
        current_data_on_device = true;
        // the data is going to change on the device, making the host data out-dated
        if(will_be_changed) current_data_on_host = false;
}

template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::uploadAsync(int stream,bool will_be_changed)
{
        // if current data is already up, do not re-upload it
//        if(current_data_on_device) return;
    if(buffer&&is_continues)
    {
           printf("Actual Buffer: %p %i\n",*buffer,*buf_size);
            if(typeid(host_type)==typeid(double))
            {
              if(typeid(dev_type)==typeid(double))
              {
                      CudaData_Upload_DoubleDouble((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
              else if(typeid(dev_type)==typeid(float))
              {
                      CudaData_Upload_DoubleFloat((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
            }
            else if(typeid(host_type)==typeid(float))
            {
              if(typeid(dev_type)==typeid(double))
              {
                      CudaData_Upload_FloatDouble((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
              else if(typeid(dev_type)==typeid(float))
              {
                      CudaData_Upload_FloatFloat((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
            }
            else if(typeid(host_type)==typeid(int))
            {
              if(typeid(dev_type)==typeid(int))
              {
                      CudaData_Upload_IntInt((void*) host_data,dev_data_array->dev_data,
                                                                                              dev_data_array->dim,mode,*buffer);
                        current_data_on_device = true;
                        if(will_be_changed) current_data_on_host = false;
                        return;
              }
            }
    }
        switch(mode)
        {
                case x:
                {
                        if(typeid(host_type) == typeid(dev_type))
                                CudaWrapper_UploadCudaDataAsync(host_data, dev_data_array->dev_data, nbytes,stream);
                        else
                        {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                          for(unsigned i=0; i<dev_data_array->dim[0]; ++i) temp_data[i] = static_cast<dev_type>(host_data[i]);
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                          CudaWrapper_UploadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes,stream);
                        }
                        break;
                }

                case xx:
                {
                        if(typeid(host_type) == typeid(dev_type))
                                CudaWrapper_UploadCudaDataAsync(host_data, dev_data_array->dev_data, nbytes,stream);
                        else
                        {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i) temp_data[i] = static_cast<dev_type>(host_data[i]);
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                                CudaWrapper_UploadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes,stream);
                        }
                        break;
                }

                case xy:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        {
                                dev_type* temp = &temp_data[i * dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        temp[j] = static_cast<dev_type>((reinterpret_cast<host_type**>(host_data))[i][j]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes,stream);
                        break;
                }

                case yx:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[j*dev_data_array->dim[0]];
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                                {
                                        temp[i] = static_cast<dev_type>(reinterpret_cast<host_type**>(host_data)[i][j]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes,stream);
                        break;
                }
                case xyz:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[(i*dev_data_array->dim[1]+j)*dev_data_array->dim[2]];
                                for(unsigned k=0; k<dev_data_array->dim[2]; ++k)
                                {
                                        temp[k] = static_cast<dev_type>(reinterpret_cast<host_type***>(host_data)[i][j][k]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes,stream);
                        break;
                }

                case xzy:
                {
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        for(unsigned k=0; k<dev_data_array->dim[2]; ++k)
                        {
                                dev_type* temp = &temp_data[(i*dev_data_array->dim[2]+k)*dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        temp[j] = static_cast<dev_type>(reinterpret_cast<host_type***>(host_data)[i][j][k]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufUploadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        CudaWrapper_UploadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes,stream);
                        break;
                }
        }
        // we have uploaded the data to the device, i.e.:
        current_data_on_device = true;
        // the data is going to change on the device, making the host data out-dated
        if(will_be_changed) current_data_on_host = false;
}

template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::download(bool will_be_changed)
{
        // if current data is already down, do not re-download it
//        if(current_data_on_host) return;
        switch(mode)
        {
                case x:
                {
                        if(typeid(host_type) == typeid(dev_type))
                                CudaWrapper_DownloadCudaData(host_data, dev_data_array->dev_data, nbytes);
                        else
                        {
                                CudaWrapper_DownloadCudaData(temp_data, dev_data_array->dev_data, nbytes);
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i) host_data[i] = static_cast<host_type>(temp_data[i]);
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufDownloadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        }
                        break;
                }

                case xx:
                {
                        if(typeid(host_type) == typeid(dev_type))
                                CudaWrapper_DownloadCudaData(host_data, dev_data_array->dev_data, nbytes);
                        else
                        {
                                CudaWrapper_DownloadCudaData(temp_data, dev_data_array->dev_data, nbytes);
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i) host_data[i] = static_cast<host_type>(temp_data[i]);
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufDownloadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        }
                        break;
                }

                case xy:
                {
                        CudaWrapper_DownloadCudaData(temp_data, dev_data_array->dev_data, nbytes);
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        {
                                dev_type* temp = &temp_data[i * dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        reinterpret_cast<host_type**>(host_data)[i][j] = static_cast<host_type>(temp[j]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufDownloadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        break;
                }

                case yx:
                {
                        CudaWrapper_DownloadCudaData(temp_data, dev_data_array->dev_data, nbytes);
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[j*dev_data_array->dim[0]];
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                                {
                                        reinterpret_cast<host_type**>(host_data)[i][j] = static_cast<host_type>(temp[i]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufDownloadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        break;
                }

                case xyz:
                {
                        CudaWrapper_DownloadCudaData(temp_data, dev_data_array->dev_data, nbytes);
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[(i * dev_data_array->dim[1]+j)*dev_data_array->dim[2]];
                                for(unsigned k=0; k<dev_data_array->dim[2]; ++k)
                                {
                                        reinterpret_cast<host_type***>(host_data)[i][j][k] = static_cast<host_type>(temp[k]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufDownloadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        break;
                }

                case xzy:
                {
                        CudaWrapper_DownloadCudaData(temp_data, dev_data_array->dev_data, nbytes);
    timespec time1,time2;
    my_gettime(CLOCK_REALTIME,&time1);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        for(unsigned k=0; k<dev_data_array->dim[2]; ++k)
                        {
                                dev_type* temp = &temp_data[(i * dev_data_array->dim[2]+k)*dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        reinterpret_cast<host_type***>(host_data)[i][j][k] = static_cast<host_type>(temp[j]);
                                }
                        }
        my_gettime(CLOCK_REALTIME,&time2);
        CudaWrapper_AddCPUBufDownloadTime(
        time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000);
                        break;
                }
        }
        // we have downloaded the data to the host, i.e.:
        current_data_on_host = true;
        // the data is going to change on the host, making the device data out-dated
        if(will_be_changed) current_data_on_device = false;
}

template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::downloadAsync(int stream)
{
        switch(mode)
        {
                case x:
                {
                        if(typeid(host_type) == typeid(dev_type))
                        {
                                CudaWrapper_DownloadCudaDataAsync(host_data, dev_data_array->dev_data, nbytes, stream);
                                CudaWrapper_SyncStream(stream);
                        }
                        else
                        {
                                CudaWrapper_DownloadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes, stream);
                                CudaWrapper_SyncStream(stream);
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i) host_data[i] = static_cast<host_type>(temp_data[i]);
                        }
                        break;
                }

                case xx:
                {
                        if(typeid(host_type) == typeid(dev_type))
                        {
                                CudaWrapper_DownloadCudaDataAsync(host_data, dev_data_array->dev_data, nbytes, stream);
                            CudaWrapper_SyncStream(stream);
                        }
                        else
                        {
                                CudaWrapper_DownloadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes, stream);
                             CudaWrapper_SyncStream(stream);
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i) host_data[i] = static_cast<host_type>(temp_data[i]);
                        }
                        break;
                }

                case xy:
                {
                        CudaWrapper_DownloadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes, stream);
                        CudaWrapper_SyncStream(stream);
                        for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                        {
                                dev_type* temp = &temp_data[i * dev_data_array->dim[1]];
                                for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                                {
                                        reinterpret_cast<host_type**>(host_data)[i][j] = static_cast<host_type>(temp[j]);
                                }
                        }
                        break;
                }

                case yx:
                {
                        CudaWrapper_DownloadCudaDataAsync(temp_data, dev_data_array->dev_data, nbytes, stream);
                        CudaWrapper_SyncStream(stream);
                        for(unsigned j=0; j<dev_data_array->dim[1]; ++j)
                        {
                                dev_type* temp = &temp_data[j*dev_data_array->dim[0]];
                                for(unsigned i=0; i<dev_data_array->dim[0]; ++i)
                                {
                                        reinterpret_cast<host_type**>(host_data)[i][j] = static_cast<host_type>(temp[i]);
                                }
                        }
                        break;
                }
        }
}


template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::memset_device(int value)
{
   CudaWrapper_Memset(dev_data_array->dev_data,value, nbytes);
}

template <typename host_type, typename dev_type, copy_mode mode>
void cCudaData<host_type, dev_type, mode>
::set_buffer(void** abuffer,int* abuf_size,bool ais_continues)
{
   buffer = abuffer;
   buf_size = abuf_size;
   unsigned nbytes_buf=(nbytes/sizeof(dev_type))*sizeof(host_type);
   if(buffer!=NULL)
   if(not((typeid(host_type) == typeid(dev_type))&&(mode == x || mode == xx)))
   {
           printf("Allocate Buffer: %p %i\n",*buffer,*buf_size);
            if(((*buffer)!=NULL)&&(*buf_size<nbytes_buf))
            CudaWrapper_FreeCudaData(*buffer,*buf_size);
            if(*buf_size<nbytes_buf)
            {*buffer=CudaWrapper_AllocCudaData(nbytes_buf);*buf_size=nbytes_buf;}
           printf("Allocate Buffer2: %p %i\n",*buffer,*buf_size);

   }
   is_continues=ais_continues;
}
#endif // _CUDA_DATA_H_
