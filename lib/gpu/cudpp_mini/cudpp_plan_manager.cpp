// -------------------------------------------------------------
// cuDPP -- CUDA Data Parallel Primitives library
// -------------------------------------------------------------
// $Revision: 3572$
// $Date: 2007-11-19 13:58:06 +0000 (Mon, 19 Nov 2007) $
// -------------------------------------------------------------
// This source code is distributed under the terms of license.txt
// in the root directory of this source distribution.
// -------------------------------------------------------------
#include "cudpp.h"
#include "cudpp_plan.h"
#include "cudpp_plan_manager.h"
#include "cudpp_maximal_launch.h"

typedef void* KernelPointer;

extern "C" size_t getNumCTAs(KernelPointer kernel)
{
    return CUDPPPlanManager::numCTAs(kernel);
}
extern "C" void compNumCTAs(KernelPointer kernel, size_t bytesDynamicSharedMem, size_t threadsPerBlock)
{
    CUDPPPlanManager::computeNumCTAs(kernel, bytesDynamicSharedMem, threadsPerBlock);
}

//! @internal Instantiate the plan manager singleton object
void CUDPPPlanManager::Instantiate()
{
    if (nullptr == m_instance)
        m_instance = new CUDPPPlanManager;
}

//! @internal Destroy the plan manager singleton object
void CUDPPPlanManager::Destroy()
{
    if (nullptr != m_instance)
    {
        delete m_instance;
        m_instance = nullptr;
    }
}

/** @brief Plan Manager destructor
* Destroys all plans as well as the plan manager.
*/
CUDPPPlanManager::~CUDPPPlanManager()
{
    std::map<CUDPPHandle,CUDPPPlan*>::iterator it;

    for (it = m_instance->plans.begin(); it != m_instance->plans.end(); it++)
    {
        CUDPPPlan* plan = it->second;
        delete plan;
        plan = nullptr;
    }
    m_instance->plans.clear();

    m_instance->numCTAsTable.clear();
}

/** @brief Add a plan to the plan manager
*
* @returns a valid CUDPPHandle if the plan was successfully added, or
* CUDPP_INVALID_HANDLE otherwise
* @param[in] plan The plan to add
*/
CUDPPHandle CUDPPPlanManager::AddPlan(CUDPPPlan* plan)
{
    Instantiate();

    std::pair<std::map<CUDPPHandle, CUDPPPlan*>::iterator, bool> ret;

    CUDPPHandle handle = (CUDPPHandle)m_instance->plans.size();
    ret = m_instance->plans.insert(std::pair<CUDPPHandle,CUDPPPlan*>(handle, plan));
    if (ret.second == true)
        return handle;
    else
        return CUDPP_INVALID_HANDLE;
}

/** @brief Remove a plan from the plan manager
*
* @returns true if the plan was successfully removed, false otherwise
* @param[in] handle The handle to the plan to remove
*/
bool CUDPPPlanManager::RemovePlan(CUDPPHandle handle)
{
    if (m_instance == nullptr)
    {
        return false;
    }

    std::map<CUDPPHandle,CUDPPPlan*>::iterator it;
    it = m_instance->plans.find(handle);

    if (it != m_instance->plans.end())
    {
        CUDPPPlan* plan = it->second;
        delete plan;
        plan = nullptr;
        m_instance->plans.erase(it);

        if (0 == m_instance->plans.size())
        {
            Destroy();
        }

        return true;
    }
    else
    {
        return false;
    }
}

/** @brief Get a plan from the plan manager by handle
*
* @returns A pointer to the plan if found, or nullptr otherwise
* @param handle The handle to the requested plan
*/
CUDPPPlan* CUDPPPlanManager::GetPlan(CUDPPHandle handle)
{
    if (m_instance == nullptr)
    {
        return nullptr;
    }

    std::map<CUDPPHandle, CUDPPPlan*>::iterator it;
    it = m_instance->plans.find(handle);
    if (it != m_instance->plans.end())
    {
        return it->second;
    }
    else
    {
        return nullptr;
    }
}

size_t CUDPPPlanManager::numCTAs(KernelPointer kernel)
{
    if (m_instance == nullptr)
    {
        return 0;
    }

    return m_instance->numCTAsTable[kernel];
}

void CUDPPPlanManager::computeNumCTAs(KernelPointer kernel, size_t bytesDynamicSharedMem, size_t threadsPerBlock)
{
    Instantiate();

    m_instance->numCTAsTable[kernel] = maxBlocks(kernel, bytesDynamicSharedMem, threadsPerBlock);
}
