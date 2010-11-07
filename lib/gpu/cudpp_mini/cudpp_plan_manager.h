// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 3572$
// $Date: 2007-11-19 13:58:06 +0000 (Mon, 19 Nov 2007) $
// ------------------------------------------------------------- 
// This source code is distributed under the terms of license.txt
// in the root directory of this source distribution.
// ------------------------------------------------------------- 
#ifndef __CUDPP_PLAN_MANAGER_H__
#define __CUDPP_PLAN_MANAGER_H__

#include <map>

class CUDPPPlan;
typedef void* KernelPointer;

/** @brief Singleton manager class for CUDPPPlan objects
  * 
  * This class manages all active plans in CUDPP.  It is a singleton class,
  * meaning that only one instance can exist.  It is created automatically the
  * first time AddPlan() is called, and destroyed when the last plan is removed
  * using RemovePlan().
  */
class CUDPPPlanManager
{
public:
    static CUDPPHandle AddPlan(CUDPPPlan* plan);
    static bool        RemovePlan(CUDPPHandle handle);
    static CUDPPPlan*  GetPlan(CUDPPHandle handle);
    
    static size_t      numCTAs(KernelPointer kernel);
    static void        computeNumCTAs(KernelPointer kernel, 
                                      size_t bytesDynamicSharedMem, 
                                      size_t threadsPerBlock);
    
protected:
    static CUDPPPlanManager* m_instance;
    std::map<CUDPPHandle, CUDPPPlan*> plans;
    std::map<void*, size_t> numCTAsTable;

private:
    

    //! @internal Instantiate the plan manager singleton object
    static void Instantiate();
    //! @internal Destroy the plan manager singleton object
    static void Destroy();

private:
    CUDPPPlanManager() {}
    CUDPPPlanManager(const CUDPPPlanManager&) {}
    ~CUDPPPlanManager();
};

#endif // __CUDPP_PLAN_MANAGER_H__
