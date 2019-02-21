/*Copyright (c) 2016 PM Larsen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <algorithm>
#include "ptm_constants.h"
#include "ptm_initialize_data.h"

namespace ptm {

#define NUM_ALLOY_TYPES 3
static uint32_t typedata[NUM_ALLOY_TYPES][3] = {
        {PTM_MATCH_FCC, PTM_ALLOY_L10,    0x00000db6},
        {PTM_MATCH_FCC, PTM_ALLOY_L12_CU, 0x00000492},
        {PTM_MATCH_FCC, PTM_ALLOY_L12_AU, 0x00001ffe},
};

static bool test_pure(int num_nbrs, int32_t* numbers)
{
        for (int i=1;i<num_nbrs + 1;i++)
                if (numbers[i] != numbers[0])
                        return false;
        return true;
}

static bool test_binary(int num_nbrs, int32_t* numbers)
{
        int a = numbers[0], b = -1;
        for (int i=1;i<num_nbrs + 1;i++)
        {
                if (numbers[i] != a)
                {
                        if (b == -1)
                                b = numbers[i];
                        else if (numbers[i] != b)
                                return false;
                }
        }

        return true;
}

static bool test_shell_structure(const refdata_t* ref, int8_t* mapping, int32_t* numbers, int num_inner)
{
        int8_t binary[PTM_MAX_POINTS];
        for (int i=0;i<ref->num_nbrs+1;i++)
                binary[i] = numbers[mapping[i]] == numbers[0] ? 0 : 1;

        for (int i=1;i<num_inner + 1;i++)
                if (binary[i] == binary[0])
                        return false;

        for (int i=num_inner+1;i<ref->num_nbrs+1;i++)
                if (binary[i] != binary[0])
                        return false;

        return true;
}

static int32_t canonical_alloy_representation(const refdata_t* ref, int8_t* mapping, int32_t* numbers)
{
        int8_t binary[PTM_MAX_POINTS];
        for (int i=0;i<ref->num_nbrs+1;i++)
                binary[i] = numbers[mapping[i]] == numbers[0] ? 0 : 1;

        int8_t temp[PTM_MAX_POINTS];
        uint32_t best = 0xFFFFFFFF;
        for (int j=0;j<ref->num_mappings;j++)
        {
                for (int i=0;i<ref->num_nbrs+1;i++)
                        temp[ref->mapping[j][i]] = binary[i];

                uint32_t code = 0;
                for (int i=0;i<ref->num_nbrs+1;i++)
                        code |= (temp[i] << i);

                best = std::min(best, code);
        }

        return best;
}

int32_t find_alloy_type(const refdata_t* ref, int8_t* mapping, int32_t* numbers)
{
        for (int i=0;i<ref->num_nbrs+1;i++)
                if (numbers[i] == -1)
                        return PTM_ALLOY_NONE;

        if (test_pure(ref->num_nbrs, numbers))
                return PTM_ALLOY_PURE;

        if (!test_binary(ref->num_nbrs, numbers))
                return PTM_ALLOY_NONE;

        uint32_t code = canonical_alloy_representation(ref, mapping, numbers);
        for (int i=0;i<NUM_ALLOY_TYPES;i++)
                if ((uint32_t)ref->type == typedata[i][0] && code == typedata[i][2])
                        return typedata[i][1];

        if (ref->type == PTM_MATCH_BCC)
                if (test_shell_structure(ref, mapping, numbers, 8))
                        return PTM_ALLOY_B2;

        if (ref->type == PTM_MATCH_DCUB || ref->type == PTM_MATCH_DHEX)
                if (test_shell_structure(ref, mapping, numbers, 4))
                        return PTM_ALLOY_SIC;


        if (ref->type == PTM_MATCH_GRAPHENE)
                if (test_shell_structure(ref, mapping, numbers, 3))
                        return PTM_ALLOY_BN;

        return PTM_ALLOY_NONE;
}

}

