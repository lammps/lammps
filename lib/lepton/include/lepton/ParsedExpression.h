#ifndef LEPTON_PARSED_EXPRESSION_H_
#define LEPTON_PARSED_EXPRESSION_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2022 Stanford University and the Authors.      *
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

namespace Lepton {

class CompiledExpression;
class ExpressionProgram;
class CompiledVectorExpression;

/**
 * This class represents the result of parsing an expression.  It provides methods for working with the
 * expression in various ways, such as evaluating it, getting the tree representation of the expresson, etc.
 */

class LEPTON_EXPORT ParsedExpression {
public:
    /**
     * Create an uninitialized ParsedExpression.  This exists so that ParsedExpressions can be put in STL containers.
     * Doing anything with it will produce an exception.
     */
    ParsedExpression();
    /**
     * Create a ParsedExpression.  Normally you will not call this directly.  Instead, use the Parser class
     * to parse expression.
     */
    ParsedExpression(const ExpressionTreeNode& rootNode);
    /**
     * Get the root node of the expression's abstract syntax tree.
     */
    const ExpressionTreeNode& getRootNode() const;
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
    /**
     * Create a new ParsedExpression which produces the same result as this one, but is faster to evaluate.
     */
    ParsedExpression optimize() const;
    /**
     * Create a new ParsedExpression which produces the same result as this one, but is faster to evaluate.
     *
     * @param variables    a map specifying values for a subset of variables that appear in the expression.
     *                     All occurrences of these variables in the expression are replaced with the values
     *                     specified.
     */
    ParsedExpression optimize(const std::map<std::string, double>& variables) const;
    /**
     * Create a new ParsedExpression which is the analytic derivative of this expression with respect to a
     * particular variable.
     *
     * @param variable     the variable with respect to which the derivate should be taken
     */
    ParsedExpression differentiate(const std::string& variable) const;
    /**
     * Create an ExpressionProgram that represents the same calculation as this expression.
     */
    ExpressionProgram createProgram() const;
    /**
     * Create a CompiledExpression that represents the same calculation as this expression.
     */
    CompiledExpression createCompiledExpression() const;
    /**
     * Create a CompiledVectorExpression that allows the expression to be evaluated efficiently
     * using the CPU's vector unit.
     *
     * @param width    the width of the vectors to evaluate it on.  The allowed values
     *                 depend on the CPU.  4 is always allowed, and 8 is allowed on
     *                 x86 processors with AVX.  Call CompiledVectorExpression::getAllowedWidths()
     *                 to query the allowed widths on the current processor.
     */
    CompiledVectorExpression createCompiledVectorExpression(int width) const;
    /**
     * Create a new ParsedExpression which is identical to this one, except that the names of some
     * variables have been changed.
     *
     * @param replacements    a map whose keys are the names of variables, and whose values are the
     *                        new names to replace them with
     */
    ParsedExpression renameVariables(const std::map<std::string, std::string>& replacements) const;
private:
    static double evaluate(const ExpressionTreeNode& node, const std::map<std::string, double>& variables);
    static ExpressionTreeNode preevaluateVariables(const ExpressionTreeNode& node, const std::map<std::string, double>& variables);
    static ExpressionTreeNode precalculateConstantSubexpressions(const ExpressionTreeNode& node, std::map<int, ExpressionTreeNode>& nodeCache);
    static ExpressionTreeNode substituteSimplerExpression(const ExpressionTreeNode& node, std::map<int, ExpressionTreeNode>& nodeCache);
    static ExpressionTreeNode differentiate(const ExpressionTreeNode& node, const std::string& variable, std::map<int, ExpressionTreeNode>& nodeCache);
    static bool isConstant(const ExpressionTreeNode& node);
    static double getConstantValue(const ExpressionTreeNode& node);
    static ExpressionTreeNode renameNodeVariables(const ExpressionTreeNode& node, const std::map<std::string, std::string>& replacements);
    ExpressionTreeNode rootNode;
};

LEPTON_EXPORT std::ostream& operator<<(std::ostream& out, const ExpressionTreeNode& node);

LEPTON_EXPORT std::ostream& operator<<(std::ostream& out, const ParsedExpression& exp);

} // namespace Lepton

#endif /*LEPTON_PARSED_EXPRESSION_H_*/
