#ifndef LEPTON_COMPILED_EXPRESSION_H_
#define LEPTON_COMPILED_EXPRESSION_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2016 Stanford University and the Authors.      *
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
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#ifdef LEPTON_USE_JIT
    #include "asmjit.h"
#endif

namespace Lepton {

class Operation;
class ParsedExpression;

/**
 * A CompiledExpression is a highly optimized representation of an expression for cases when you want to evaluate
 * it many times as quickly as possible.  You should treat it as an opaque object; none of the internal representation
 * is visible.
 * 
 * A CompiledExpression is created by calling createCompiledExpression() on a ParsedExpression.
 * 
 * WARNING: CompiledExpression is NOT thread safe.  You should never access a CompiledExpression from two threads at
 * the same time.
 */

class LEPTON_EXPORT CompiledExpression {
public:
    CompiledExpression();
    CompiledExpression(const CompiledExpression& expression);
    ~CompiledExpression();
    CompiledExpression& operator=(const CompiledExpression& expression);
    /**
     * Get the names of all variables used by this expression.
     */
    const std::set<std::string>& getVariables() const;
    /**
     * Get a reference to the memory location where the value of a particular variable is stored.  This can be used
     * to set the value of the variable before calling evaluate().
     */
    double& getVariableReference(const std::string& name);
    /**
     * You can optionally specify the memory locations from which the values of variables should be read.
     * This is useful, for example, when several expressions all use the same variable.  You can then set
     * the value of that variable in one place, and it will be seen by all of them.
     */
    void setVariableLocations(std::map<std::string, double*>& variableLocations);
    /**
     * Evaluate the expression.  The values of all variables should have been set before calling this.
     */
    double evaluate() const;
private:
    friend class ParsedExpression;
    CompiledExpression(const ParsedExpression& expression);
    void compileExpression(const ExpressionTreeNode& node, std::vector<std::pair<ExpressionTreeNode, int> >& temps);
    int findTempIndex(const ExpressionTreeNode& node, std::vector<std::pair<ExpressionTreeNode, int> >& temps);
    std::map<std::string, double*> variablePointers;
    std::vector<std::pair<double*, double*> > variablesToCopy;
    std::vector<std::vector<int> > arguments;
    std::vector<int> target;
    std::vector<Operation*> operation;
    std::map<std::string, int> variableIndices;
    std::set<std::string> variableNames;
    mutable std::vector<double> workspace;
    mutable std::vector<double> argValues;
    std::map<std::string, double> dummyVariables;
    void* jitCode;
#ifdef LEPTON_USE_JIT
    void generateJitCode();
    void generateSingleArgCall(asmjit::X86Compiler& c, asmjit::X86XmmVar& dest, asmjit::X86XmmVar& arg, double (*function)(double));
    std::vector<double> constants;
    asmjit::JitRuntime runtime;
#endif
};

} // namespace Lepton

#endif /*LEPTON_COMPILED_EXPRESSION_H_*/
