// -------------------------------------------------------------
// CUDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
//  $Revision: 5632 $
//  $Date: 2009-07-01 14:36:01 +1000 (Wed, 01 Jul 2009) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt in
// the root directory of this source distribution.
// ------------------------------------------------------------- 

/**
 * @file
 * vector_kernel.cu
 * 
 * @brief CUDA kernel methods for basic operations on vectors.  
 * 
 * CUDA kernel methods for basic operations on vectors.  
 * 
 * Examples: 
 * - vectorAddConstant(): d_vector + constant 
 * - vectorAddUniform():  d_vector + uniform (per-block constants)
 * - vectorAddVectorVector(): d_vector + d_vector
 */

// MJH: these functions assume there are 2N elements for N threads.  
// Is this always going to be a good idea?  There may be cases where
// we have as many threads as elements, but for large problems
// we are probably limited by max CTA size for simple kernels like 
// this so we should process multiple elements per thread.
// we may want to extend these with looping versions that process 
// many elements per thread.

#include "cudpp_util.h"
#include "sharedmem.h"
#include "cudpp.h"

/** \addtogroup cudpp_kernel
  * @{
  */

/** @name Vector Functions
 * CUDA kernel methods for basic operations on vectors.  
 * @{
 */

/** @brief Adds a constant value to all values in the input d_vector
 *  
 * Each thread adds two pairs of elements.
 * @todo Test this function -- it is currently not yet used.
 *
 * @param[in,out] d_vector The array of elements to be modified
 * @param[in] constant The constant value to be added to elements of 
 * \a d_vector
 * @param[in] n The number of elements in the d_vector to be modified
 * @param[in] baseIndex An optional offset to the beginning of the
 * elements in the input array to be processed
 */
template <class T>
__global__  void vectorAddConstant(T   *d_vector, 
                                   T   constant, 
                                   int n, 
                                   int baseIndex)
{
    // Compute this thread's output address
    unsigned int address = baseIndex + threadIdx.x +
        __mul24(blockIdx.x, (blockDim.x << 1)); 

    // note two adds per thread: one in first half of the block, one in last
    d_vector[address]              += constant;
    d_vector[address + blockDim.x] += (threadIdx.x + blockDim.x < n) * constant;
}

 /** @brief Add a uniform value to each data element of an array
  *
  * This function reads one value per CTA from \a d_uniforms into shared
  * memory and adds that value to all values "owned" by the CTA in \a
  * d_vector.  Each thread adds two pairs of values.
  *
  * @param[out] d_vector The d_vector whose values will have the uniform added
  * @param[in] d_uniforms The array of uniform values (one per CTA)
  * @param[in] numElements The number of elements in \a d_vector to process
  * @param[in] blockOffset an optional offset to the beginning of this block's
  * data.
  * @param[in] baseIndex an optional offset to the beginning of the array 
  * within \a d_vector.
  */
template <class T>
__global__ void vectorAddUniform(T       *d_vector, 
                                 const T *d_uniforms, 
                                 int     numElements, 
                                 int     blockOffset, 
                                 int     baseIndex)
{
    __shared__ T uni;
    // Get this block's uniform value from the uniform array in device memory
    // We store it in shared memory so that the hardware's shared memory 
    // broadcast capability can be used to share among all threads in each warp
    // in a single cycle
    if (threadIdx.x == 0)
    {
        uni = d_uniforms[blockIdx.x + __mul24(gridDim.x, blockIdx.y) + blockOffset];
    }

    // Compute this thread's output address
    int width = __mul24(gridDim.x,(blockDim.x << 1));

    unsigned int address = baseIndex + __mul24(width, blockIdx.y)
        + threadIdx.x + __mul24(blockIdx.x, (blockDim.x << 1)); 

    __syncthreads();

    // note two adds per thread: one in first half of the block, one in last
    d_vector[address]              += uni;
    if (threadIdx.x + blockDim.x < numElements) d_vector[address + blockDim.x] += uni;
}


/** @brief Add a uniform value to each data element of an array (vec4 version)
  *
  * This function reads one value per CTA from \a d_uniforms into shared
  * memory and adds that value to all values "owned" by the CTA in \a d_vector.  
  * Each thread adds the uniform value to eight values in \a d_vector.
  *
  * @param[out] d_vector The d_vector whose values will have the uniform added
  * @param[in] d_uniforms The array of uniform values (one per CTA)
  * @param[in] numElements The number of elements in \a d_vector to process
  * @param[in] vectorRowPitch For 2D arrays, the pitch (in elements) of the 
  * rows of \a d_vector.
  * @param[in] uniformRowPitch For 2D arrays, the pitch (in elements) of the 
  * rows of \a d_uniforms.
  * @param[in] blockOffset an optional offset to the beginning of this block's
  * data.
  * @param[in] baseIndex an optional offset to the beginning of the array 
  * within \a d_vector.
  */
template <class T, CUDPPOperator op, int elementsPerThread>
__global__ void vectorAddUniform4(T       *d_vector, 
                                  const T *d_uniforms, 
                                  int      numElements,             
                                  int      vectorRowPitch,     // width of input array in elements
                                  int      uniformRowPitch,    // width of uniform array in elements
                                  int      blockOffset, 
                                  int      baseIndex)
{
    __shared__ T uni;
    // Get this block's uniform value from the uniform array in device memory
    // We store it in shared memory so that the hardware's shared memory 
    // broadcast capability can be used to share among all threads in each warp
    // in a single cycle
    if (threadIdx.x == 0)
    {
        uni = d_uniforms[blockIdx.x + __umul24(uniformRowPitch, blockIdx.y) + blockOffset];
    }

    // Compute this thread's output address
    //int width = __mul24(gridDim.x,(blockDim.x << 1));
   
    unsigned int address = baseIndex + __umul24(vectorRowPitch, blockIdx.y)
        + threadIdx.x + __umul24(blockIdx.x, (blockDim.x * elementsPerThread)); 
    numElements += __umul24(vectorRowPitch, blockIdx.y);

    __syncthreads();

    switch (op)
    {
    case CUDPP_ADD:
        for (int i = 0; i < elementsPerThread && address < numElements; i++)
        {
            d_vector[address] += uni;
            address += blockDim.x;
        }
        break;

    case CUDPP_MULTIPLY:
        for (int i = 0; i < elementsPerThread && address < numElements; i++)
        {
            d_vector[address] *= uni;
            address += blockDim.x;
        }
        break;

    case CUDPP_MAX:
        for (int i = 0; i < elementsPerThread && address < numElements; i++)
        {
            d_vector[address] = max(d_vector[address], uni);
            address += blockDim.x;
        }
        break;

    case CUDPP_MIN:
        for (int i = 0; i < elementsPerThread && address < numElements; i++)
        {
            d_vector[address] = min(d_vector[address], uni);
            address += blockDim.x;
        }
        break;
    default:
        break;
    }    
}

/** @brief Adds together two vectors
 *  
 * Each thread adds two pairs of elements.
 * @todo Test this function -- it is currently not yet used.
 *
 * @param[out] d_vectorA The left operand array and the result
 * @param[in] d_vectorB The right operand array
 * @param[in] numElements The number of elements in the vectors to be added.
 * @param[in] baseIndex An optional offset to the beginning of the
 * elements in the input arrays to be processed
 */
template <class T>
__global__ void vectorAddVector(T       *d_vectorA,        // A += B
                                const T *d_vectorB,
                                int     numElements,
                                int     baseIndex)
{
    // Compute this thread's output address
    unsigned int address = baseIndex + threadIdx.x +
        __mul24(blockIdx.x, (blockDim.x << 1)); 

    // note two adds per thread: one in first half of the block, one in last
    d_vectorA[address]              += d_vectorB[address];
    d_vectorA[address + blockDim.x] += 
        (threadIdx.x + blockDim.x < numElements) * d_vectorB[address];
}

/** @brief Add a uniform value to data elements of an array (vec4 version)
  *
  * This function reads one value per CTA from \a d_uniforms into shared
  * memory and adds that value to values "owned" by the CTA in \a d_vector.
  * The uniform value is added to only those values "owned" by the CTA which
  * have an index less than d_maxIndex. If d_maxIndex for that CTA is UINT_MAX
  * it adds the uniform to all values "owned" by the CTA.
  * Each thread adds the uniform value to eight values in \a d_vector.
  *
  * @param[out] d_vector The d_vector whose values will have the uniform added
  * @param[in] d_uniforms The array of uniform values (one per CTA)
  * @param[in] d_maxIndices The array of maximum indices (one per CTA). This is
  *            index upto which the uniform would be added. If this is UINT_MAX
  *            the uniform is added to all elements of the CTA. This index is
  *            1-based.
  * @param[in] numElements The number of elements in \a d_vector to process
  * @param[in] blockOffset an optional offset to the beginning of this block's
  * data.
  * @param[in] baseIndex an optional offset to the beginning of the array 
  * within \a d_vector.
  */
template <class T, CUDPPOperator oper, bool isLastBlockFull>
__global__ void vectorSegmentedAddUniform4(T                  *d_vector, 
                                           const T            *d_uniforms, 
                                           const unsigned int *d_maxIndices,
                                           unsigned int       numElements,
                                           int                blockOffset, 
                                           int                baseIndex)
{
    __shared__ T uni[2];

    unsigned int blockAddress = 
        blockIdx.x + __mul24(gridDim.x, blockIdx.y) + blockOffset;

    // Get this block's uniform value from the uniform array in device memory
    // We store it in shared memory so that the hardware's shared memory 
    // broadcast capability can be used to share among all threads in each warp
    // in a single cycle
    
    if (threadIdx.x == 0)
    {
        if (blockAddress > 0)
            uni[0] = d_uniforms[blockAddress-1];
        else
            uni[0] = Operator<T, oper>::identity(); 
        
        // Tacit assumption that T is four-byte wide
        uni[1] = (T)(d_maxIndices[blockAddress]);
    }

    // Compute this thread's output address
    int width = __mul24(gridDim.x,(blockDim.x << 1));

    unsigned int address = baseIndex + __mul24(width, blockIdx.y)
                           + threadIdx.x + __mul24(blockIdx.x, (blockDim.x << 3)); 

    __syncthreads();

    unsigned int maxIndex = (unsigned int)(uni[1]);

    bool isLastBlock = (blockIdx.x == (gridDim.x-1));
     
    if (maxIndex < UINT_MAX)
    {
        // Since maxIndex is a 1 based index
        --maxIndex;
        bool leftLess = address < maxIndex;
        bool rightLess = (address + 7 * blockDim.x) < maxIndex;

        if (leftLess)
        {
            if (rightLess)
            {
                for (unsigned int i = 0; i < 8; ++i)
                    d_vector[address + i * blockDim.x] = 
                        Operator<T, oper>::op(d_vector[address + i * blockDim.x], uni[0]);
            }
            else
            {
                for (unsigned int i=0; i < 8; ++i)
                {
                    if (address < maxIndex)
                        d_vector[address] = 
                            Operator<T, oper>::op(d_vector[address], uni[0]);

                    address += blockDim.x;
                }
            }
        }
    }
    else
    {
        if (!isLastBlockFull && isLastBlock)
        {
            for (unsigned int i = 0; i < 8; ++i)
            {
                if (address < numElements)
                    d_vector[address] = 
                        Operator<T, oper>::op(d_vector[address], uni[0]);
                
                address += blockDim.x;
            }
        }
        else
        {
            for (unsigned int i=0; i<8; ++i)
            {
                d_vector[address] = 
                    Operator<T, oper>::op(d_vector[address], uni[0]);
                
                address += blockDim.x;
            }            
        }
    }
}

/** @brief Add a uniform value to data elements of an array (vec4 version)
  *
  * This function reads one value per CTA from \a d_uniforms into shared
  * memory and adds that value to values "owned" by the CTA in \a d_vector.
  * The uniform value is added to only those values "owned" by the CTA which
  * have an index greater than d_minIndex. If d_minIndex for that CTA is 0
  * it adds the uniform to all values "owned" by the CTA.
  * Each thread adds the uniform value to eight values in \a d_vector.
  *
  * @param[out] d_vector The d_vector whose values will have the uniform added
  * @param[in] d_uniforms The array of uniform values (one per CTA)
  * @param[in] d_minIndices The array of minimum indices (one per CTA). The
  *            uniform is added to the right of this index (that is, to every index
  *            that is greater than this index). If this is 0, the uniform is 
  *            added to all elements of the CTA. This index is 1-based to
  *            prevent overloading of what 0 means. In our case it means
  *            absence of a flag. But if the first element of a CTA has
  *            flag the index will also be 0. Hence we use 1-based indices
  *            so the index is 1 in the latter case.
  * @param[in] numElements The number of elements in \a d_vector to process
  * @param[in] blockOffset an optional offset to the beginning of this block's
  * data.
  * @param[in] baseIndex an optional offset to the beginning of the array 
  * within \a d_vector.
  *
  */
template <class T, CUDPPOperator oper, bool isLastBlockFull>
__global__ void vectorSegmentedAddUniformToRight4(T                  *d_vector, 
                                                  const T            *d_uniforms, 
                                                  const unsigned int *d_minIndices,
                                                  unsigned int       numElements,
                                                  int                blockOffset, 
                                                  int                baseIndex)
{
    __shared__ T uni[2];

    unsigned int blockAddress = 
        blockIdx.x + __mul24(gridDim.x, blockIdx.y) + blockOffset;

    // Get this block's uniform value from the uniform array in device memory
    // We store it in shared memory so that the hardware's shared memory 
    // broadcast capability can be used to share among all threads in each warp
    // in a single cycle
    
    if (threadIdx.x == 0)
    {
        // FIXME - blockAddress test here is incompatible with how it is calculated
        // above
        if (blockAddress < (gridDim.x-1))
            uni[0] = d_uniforms[blockAddress+1];
        else
            uni[0] = Operator<T, oper>::identity(); 
        
        // Tacit assumption that T is four-byte wide
        uni[1] = (T)(d_minIndices[blockAddress]);
    }

    // Compute this thread's output address
    int width = __mul24(gridDim.x,(blockDim.x << 1));

    unsigned int address = baseIndex + __mul24(width, blockIdx.y)
                           + threadIdx.x + __mul24(blockIdx.x, (blockDim.x << 3)); 

    __syncthreads();

    unsigned int minIndex = (unsigned int)(uni[1]);

    bool isLastBlock = (blockIdx.x == (gridDim.x-1));
     
    if (minIndex > 0)
    {
        // Since minIndex is a 1 based index
        --minIndex;
        bool leftInRange = address > minIndex;
        bool rightInRange = (address + 7 * blockDim.x) > minIndex;

        if (rightInRange)
        {
            if (leftInRange)
            {
                for (unsigned int i = 0; i < 8; ++i)
                    d_vector[address + i * blockDim.x] = 
                        Operator<T, oper>::op(d_vector[address + i * blockDim.x], uni[0]);
            }
            else
            {
                for (unsigned int i=0; i < 8; ++i)
                {
                    if (address > minIndex)
                        d_vector[address] = 
                            Operator<T, oper>::op(d_vector[address], uni[0]);

                    address += blockDim.x;
                }
            }
        }
    }
    else
    {
        if (!isLastBlockFull && isLastBlock)
        {
            for (unsigned int i = 0; i < 8; ++i)
            {
                if (address < numElements)
                    d_vector[address] = 
                        Operator<T, oper>::op(d_vector[address], uni[0]);
                
                address += blockDim.x;
            }
        }
        else
        {
            for (unsigned int i=0; i<8; ++i)
            {
                d_vector[address] = 
                    Operator<T, oper>::op(d_vector[address], uni[0]);
                
                address += blockDim.x;
            }            
        }
    }
}

/** @} */ // end d_vector functions
/** @} */ // end cudpp_kernel
