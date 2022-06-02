#ifndef LEPTON_VECTOR_EXPRESSION_H_
#define LEPTON_VECTOR_EXPRESSION_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2022 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ExpressionTreeNode.h"
#include "windowsIncludes.h"
#include <array>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#ifdef LEPTON_USE_JIT
#if defined(__ARM__) || defined(__ARM64__)
#include "asmjit/a64.h"
#else
#include "asmjit/x86.h"
#endif
#endif

namespace Lepton {

class Operation;
class ParsedExpression;

/**
 * A CompiledVectorExpression is a highly optimized representation of an expression for cases when you want to evaluate
 * it many times as quickly as possible.  It is similar to CompiledExpression, with the extra feature that it uses the CPU's
 * vector unit (AVX on x86, NEON on ARM) to evaluate the expression for multiple sets of arguments at once.  It also differs
 * from CompiledExpression and ParsedExpression in using single precision rather than double precision to evaluate the expression.
 * You should treat it as an opaque object; none of the internal representation is visible.
 *
 * A CompiledVectorExpression is created by calling createCompiledVectorExpression() on a ParsedExpression.  When you create
 * it, you must specify the width of the vectors on which to compute the expression.  The allowed widths depend on the type of
 * CPU it is running on.  4 is always allowed, and 8 is allowed on x86 processors with AVX.  Call getAllowedWidths() to query
 * the allowed values.
 *
 * WARNING: CompiledVectorExpression is NOT thread safe.  You should never access a CompiledVectorExpression from two threads at
 * the same time.
 */

class LEPTON_EXPORT CompiledVectorExpression {
public:
    CompiledVectorExpression();
    CompiledVectorExpression(const CompiledVectorExpression& expression);
    ~CompiledVectorExpression();
    CompiledVectorExpression& operator=(const CompiledVectorExpression& expression);
    /**
     * Get the width of the vectors on which the expression is computed.
     */
    int getWidth() const;
    /**
     * Get the names of all variables used by this expression.
     */
    const std::set<std::string>& getVariables() const;
    /**
     * Get a pointer to the memory location where the value of a particular variable is stored.  This can be used
     * to set the value of the variable before calling evaluate().
     *
     * @param name    the name of the variable to query
     * @return a pointer to N floating point values, where N is the vector width
     */
    float* getVariablePointer(const std::string& name);
    /**
     * You can optionally specify the memory locations from which the values of variables should be read.
     * This is useful, for example, when several expressions all use the same variable.  You can then set
     * the value of that variable in one place, and it will be seen by all of them.  The location should
     * be a pointer to N floating point values, where N is the vector width.
     */
    void setVariableLocations(std::map<std::string, float*>& variableLocations);
    /**
     * Evaluate the expression.  The values of all variables should have been set before calling this.
     *
     * @return a pointer to N floating point values, where N is the vector width
     */
    const float* evaluate() const;
    /**
     * Get the list of vector widths that are supported on the current processor.
     */
    static const std::vector<int>& getAllowedWidths();
private:
    friend class ParsedExpression;
    CompiledVectorExpression(const ParsedExpression& expression, int width);
    void compileExpression(const ExpressionTreeNode& node, std::vector<std::pair<ExpressionTreeNode, int> >& temps, int& workspaceSize);
    int findTempIndex(const ExpressionTreeNode& node, std::vector<std::pair<ExpressionTreeNode, int> >& temps);
    int width;
    std::map<std::string, float*> variablePointers;
    std::vector<std::pair<float*, float*> > variablesToCopy;
    std::vector<std::vector<int> > arguments;
    std::vector<int> target;
    std::vector<Operation*> operation;
    std::map<std::string, int> variableIndices;
    std::set<std::string> variableNames;
    mutable std::vector<float> workspace;
    mutable std::vector<double> argValues;
    std::map<std::string, double> dummyVariables;
    void (*jitCode)();
#ifdef LEPTON_USE_JIT
    void findPowerGroups(std::vector<std::vector<int> >& groups, std::vector<std::vector<int> >& groupPowers, std::vector<int>& stepGroup);
    void generateJitCode();
#if defined(__ARM__) || defined(__ARM64__)
    void generateSingleArgCall(asmjit::a64::Compiler& c, asmjit::arm::Vec& dest, asmjit::arm::Vec& arg, float (*function)(float));
    void generateTwoArgCall(asmjit::a64::Compiler& c, asmjit::arm::Vec& dest, asmjit::arm::Vec& arg1, asmjit::arm::Vec& arg2, float (*function)(float, float));
#else
    void generateSingleArgCall(asmjit::x86::Compiler& c, asmjit::x86::Ymm& dest, asmjit::x86::Ymm& arg, float (*function)(float));
    void generateTwoArgCall(asmjit::x86::Compiler& c, asmjit::x86::Ymm& dest, asmjit::x86::Ymm& arg1, asmjit::x86::Ymm& arg2, float (*function)(float, float));
#endif
    std::vector<float> constants;
    asmjit::JitRuntime runtime;
#endif
};

} // namespace Lepton

#endif /*LEPTON_VECTOR_EXPRESSION_H_*/
