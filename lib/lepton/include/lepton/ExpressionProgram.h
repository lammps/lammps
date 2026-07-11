#ifndef LEPTON_EXPRESSION_PROGRAM_H_
#define LEPTON_EXPRESSION_PROGRAM_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2018 Stanford University and the Authors.      *
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
#include <string>
#include <vector>

namespace Lepton {

class ParsedExpression;

/**
 * An ExpressionProgram is a linear sequence of Operations for evaluating an expression.  The evaluation
 * is done with a stack.  The arguments to each Operation are first taken off the stack in order, then it is
 * evaluated and the result is pushed back onto the stack.  At the end, the stack contains a single value,
 * which is the value of the expression.
 *
 * An ExpressionProgram is created by calling createProgram() on a ParsedExpression.
 */

class LEPTON_EXPORT ExpressionProgram {
public:
    ExpressionProgram();
    ExpressionProgram(const ExpressionProgram& program);
    ~ExpressionProgram();
    ExpressionProgram& operator=(const ExpressionProgram& program);
    /**
     * Get the number of Operations that make up this program.
     */
    int getNumOperations() const;
    /**
     * Get an Operation in this program.
     */
    const Operation& getOperation(int index) const;
    /**
     * Change an Operation in this program.
     *
     * The Operation must have been allocated on the heap with the "new" operator.
     * The ExpressionProgram assumes ownership of it and will delete it when it
     * is no longer needed.
     */
    void setOperation(int index, Operation* operation);
    /**
     * Get the size of the stack needed to execute this program.  This is the largest number of elements present
     * on the stack at any point during evaluation.
     */
    int getStackSize() const;
    /**
     * Evaluate the expression.  If the expression involves any variables, this method will throw an exception.
     */
    double evaluate() const;
    /**
     * Evaluate the expression.
     *
     * @param variables    a map specifying the values of all variables that appear in the expression.  If any
     *                     variable appears in the expression but is not included in this map, an exception
     *                     will be thrown.
     */
    double evaluate(const std::map<std::string, double>& variables) const;
private:
    friend class ParsedExpression;
    ExpressionProgram(const ParsedExpression& expression);
    void buildProgram(const ExpressionTreeNode& node);
    std::vector<Operation*> operations;
    int maxArgs, stackSize;
};

} // namespace Lepton

#endif /*LEPTON_EXPRESSION_PROGRAM_H_*/
