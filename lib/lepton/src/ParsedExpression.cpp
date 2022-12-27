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

#include "lepton/ParsedExpression.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CompiledVectorExpression.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/Operation.h"
#include <limits>
#include <vector>

using namespace Lepton;
using namespace std;

ParsedExpression::ParsedExpression() : rootNode(ExpressionTreeNode()) {
}

ParsedExpression::ParsedExpression(const ExpressionTreeNode& rootNode) : rootNode(rootNode) {
}

const ExpressionTreeNode& ParsedExpression::getRootNode() const {
    if (&rootNode.getOperation() == NULL)
        throw Exception("Illegal call to an initialized ParsedExpression");
    return rootNode;
}

double ParsedExpression::evaluate() const {
    return evaluate(getRootNode(), map<string, double>());
}

double ParsedExpression::evaluate(const map<string, double>& variables) const {
    return evaluate(getRootNode(), variables);
}

double ParsedExpression::evaluate(const ExpressionTreeNode& node, const map<string, double>& variables) {
    int numArgs = (int) node.getChildren().size();
    vector<double> args(max(numArgs, 1));
    for (int i = 0; i < numArgs; i++)
        args[i] = evaluate(node.getChildren()[i], variables);
    return node.getOperation().evaluate(&args[0], variables);
}

ParsedExpression ParsedExpression::optimize() const {
    ExpressionTreeNode result = getRootNode();
    vector<const ExpressionTreeNode*> examples;
    result.assignTags(examples);
    map<int, ExpressionTreeNode> nodeCache;
    result = precalculateConstantSubexpressions(result, nodeCache);
    while (true) {
        examples.clear();
        result.assignTags(examples);
        nodeCache.clear();
        ExpressionTreeNode simplified = substituteSimplerExpression(result, nodeCache);
        if (simplified == result)
            break;
        result = simplified;
    }
    return ParsedExpression(result);
}

ParsedExpression ParsedExpression::optimize(const map<string, double>& variables) const {
    ExpressionTreeNode result = preevaluateVariables(getRootNode(), variables);
    vector<const ExpressionTreeNode*> examples;
    result.assignTags(examples);
    map<int, ExpressionTreeNode> nodeCache;
    result = precalculateConstantSubexpressions(result, nodeCache);
    while (true) {
        examples.clear();
        result.assignTags(examples);
        nodeCache.clear();
        ExpressionTreeNode simplified = substituteSimplerExpression(result, nodeCache);
        if (simplified == result)
            break;
        result = simplified;
    }
    return ParsedExpression(result);
}

ExpressionTreeNode ParsedExpression::preevaluateVariables(const ExpressionTreeNode& node, const map<string, double>& variables) {
    if (node.getOperation().getId() == Operation::VARIABLE) {
        const Operation::Variable& var = dynamic_cast<const Operation::Variable&>(node.getOperation());
        map<string, double>::const_iterator iter = variables.find(var.getName());
        if (iter == variables.end())
            return node;
        return ExpressionTreeNode(new Operation::Constant(iter->second));
    }
    vector<ExpressionTreeNode> children(node.getChildren().size());
    for (int i = 0; i < (int) children.size(); i++)
        children[i] = preevaluateVariables(node.getChildren()[i], variables);
    return ExpressionTreeNode(node.getOperation().clone(), children);
}

ExpressionTreeNode ParsedExpression::precalculateConstantSubexpressions(const ExpressionTreeNode& node, map<int, ExpressionTreeNode>& nodeCache) {
    auto cached = nodeCache.find(node.tag);
    if (cached != nodeCache.end())
        return cached->second;
    vector<ExpressionTreeNode> children(node.getChildren().size());
    for (int i = 0; i < (int) children.size(); i++)
        children[i] = precalculateConstantSubexpressions(node.getChildren()[i], nodeCache);
    ExpressionTreeNode result = ExpressionTreeNode(node.getOperation().clone(), children);
    if (node.getOperation().getId() == Operation::VARIABLE || node.getOperation().getId() == Operation::CUSTOM) {
        nodeCache[node.tag] = result;
        return result;
    }
    for (int i = 0; i < (int) children.size(); i++)
        if (children[i].getOperation().getId() != Operation::CONSTANT) {
            nodeCache[node.tag] = result;
            return result;
        }
    result = ExpressionTreeNode(new Operation::Constant(evaluate(result, map<string, double>())));
    nodeCache[node.tag] = result;
    return result;
}

ExpressionTreeNode ParsedExpression::substituteSimplerExpression(const ExpressionTreeNode& node, map<int, ExpressionTreeNode>& nodeCache) {
    vector<ExpressionTreeNode> children(node.getChildren().size());
    for (int i = 0; i < (int) children.size(); i++) {
        const ExpressionTreeNode& child = node.getChildren()[i];
        auto cached = nodeCache.find(child.tag);
        if (cached == nodeCache.end()) {
            children[i] = substituteSimplerExpression(child, nodeCache);
            nodeCache[child.tag] = children[i];
        }
        else
            children[i] = cached->second;
    }

    // Collect some info on constant expressions in children
    bool first_const = children.size() > 0 && isConstant(children[0]); // is first child constant?
    bool second_const = children.size() > 1 && isConstant(children[1]); // is second child constant?
    double first, second; // if yes, value of first and second child
    first = second = 0.0;
    if (first_const)
        first = getConstantValue(children[0]);
    if (second_const)
        second = getConstantValue(children[1]);

    switch (node.getOperation().getId()) {
        case Operation::ADD:
        {
            if (first_const) {
                if (first == 0.0) { // Add 0
                    return children[1];
                } else { // Add a constant
                    return ExpressionTreeNode(new Operation::AddConstant(first), children[1]);
                }
            }
            if (second_const) {
                if (second == 0.0) { // Add 0
                    return children[0];
                } else { // Add a constant
                    return ExpressionTreeNode(new Operation::AddConstant(second), children[0]);
                }
            }
            if (children[1].getOperation().getId() == Operation::NEGATE) // a+(-b) = a-b
                return ExpressionTreeNode(new Operation::Subtract(), children[0], children[1].getChildren()[0]);
            if (children[0].getOperation().getId() == Operation::NEGATE) // (-a)+b = b-a
                return ExpressionTreeNode(new Operation::Subtract(), children[1], children[0].getChildren()[0]);
            break;
        }
        case Operation::SUBTRACT:
        {
            if (children[0] == children[1])
                return ExpressionTreeNode(new Operation::Constant(0.0)); // Subtracting anything from itself is 0
            if (first_const) {
                if (first == 0.0) // Subtract from 0
                    return ExpressionTreeNode(new Operation::Negate(), children[1]);
            }
            if (second_const) {
                if (second == 0.0) { // Subtract 0
                    return children[0];
                } else { // Subtract a constant
                    return ExpressionTreeNode(new Operation::AddConstant(-second), children[0]);
                }
            }
            if (children[1].getOperation().getId() == Operation::NEGATE) // a-(-b) = a+b
                return ExpressionTreeNode(new Operation::Add(), children[0], children[1].getChildren()[0]);
            break;
        }
        case Operation::MULTIPLY:
        {
            if ((first_const && first == 0.0) || (second_const && second == 0.0)) // Multiply by 0
                return ExpressionTreeNode(new Operation::Constant(0.0));
            if (first_const && first == 1.0) // Multiply by 1
                return children[1];
            if (second_const && second == 1.0) // Multiply by 1
                return children[0];
            if (first_const) { // Multiply by a constant
                if (children[1].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Combine two multiplies into a single one
                    return ExpressionTreeNode(new Operation::MultiplyConstant(first*dynamic_cast<const Operation::MultiplyConstant*>(&children[1].getOperation())->getValue()), children[1].getChildren()[0]);
                return ExpressionTreeNode(new Operation::MultiplyConstant(first), children[1]);
            }
            if (second_const) { // Multiply by a constant
                if (children[0].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Combine two multiplies into a single one
                    return ExpressionTreeNode(new Operation::MultiplyConstant(second*dynamic_cast<const Operation::MultiplyConstant*>(&children[0].getOperation())->getValue()), children[0].getChildren()[0]);
                return ExpressionTreeNode(new Operation::MultiplyConstant(second), children[0]);
            }
            if (children[0].getOperation().getId() == Operation::NEGATE && children[1].getOperation().getId() == Operation::NEGATE) // The two negations cancel
                return ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[0], children[1].getChildren()[0]);
            if (children[0].getOperation().getId() == Operation::NEGATE && children[1].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Negate the constant
                return ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[0], ExpressionTreeNode(new Operation::MultiplyConstant(-dynamic_cast<const Operation::MultiplyConstant*>(&children[1].getOperation())->getValue()), children[1].getChildren()[0]));
            if (children[1].getOperation().getId() == Operation::NEGATE && children[0].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Negate the constant
                return ExpressionTreeNode(new Operation::Multiply(), ExpressionTreeNode(new Operation::MultiplyConstant(-dynamic_cast<const Operation::MultiplyConstant*>(&children[0].getOperation())->getValue()), children[0].getChildren()[0]), children[1].getChildren()[0]);
            if (children[0].getOperation().getId() == Operation::NEGATE) // Pull the negation out so it can possibly be optimized further
                return ExpressionTreeNode(new Operation::Negate(), ExpressionTreeNode(new Operation::Multiply(), children[0].getChildren()[0], children[1]));
            if (children[1].getOperation().getId() == Operation::NEGATE) // Pull the negation out so it can possibly be optimized further
                return ExpressionTreeNode(new Operation::Negate(), ExpressionTreeNode(new Operation::Multiply(), children[0], children[1].getChildren()[0]));
            if (children[1].getOperation().getId() == Operation::RECIPROCAL) // a*(1/b) = a/b
                return ExpressionTreeNode(new Operation::Divide(), children[0], children[1].getChildren()[0]);
            if (children[0].getOperation().getId() == Operation::RECIPROCAL) // (1/a)*b = b/a
                return ExpressionTreeNode(new Operation::Divide(), children[1], children[0].getChildren()[0]);
            if (children[0] == children[1])
                return ExpressionTreeNode(new Operation::Square(), children[0]); // x*x = square(x)
            if (children[0].getOperation().getId() == Operation::SQUARE && children[0].getChildren()[0] == children[1])
                return ExpressionTreeNode(new Operation::Cube(), children[1]); // x*x*x = cube(x)
            if (children[1].getOperation().getId() == Operation::SQUARE && children[1].getChildren()[0] == children[0])
                return ExpressionTreeNode(new Operation::Cube(), children[0]); // x*x*x = cube(x)
            break;
        }
        case Operation::DIVIDE:
        {
            if (children[0] == children[1])
                return ExpressionTreeNode(new Operation::Constant(1.0)); // Dividing anything from itself is 0
            if (first_const && first == 0.0) // 0 divided by something
                return ExpressionTreeNode(new Operation::Constant(0.0));
            if (first_const && first == 1.0) // 1 divided by something
                return ExpressionTreeNode(new Operation::Reciprocal(), children[1]);
            if (second_const && second == 1.0) // Divide by 1
                return children[0];
            if (second_const) {
                if (children[0].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Combine a multiply and a divide into one multiply
                    return ExpressionTreeNode(new Operation::MultiplyConstant(dynamic_cast<const Operation::MultiplyConstant*>(&children[0].getOperation())->getValue()/second), children[0].getChildren()[0]);
                return ExpressionTreeNode(new Operation::MultiplyConstant(1.0/second), children[0]); // Replace a divide with a multiply
            }
            if (children[0].getOperation().getId() == Operation::NEGATE && children[1].getOperation().getId() == Operation::NEGATE) // The two negations cancel
                return ExpressionTreeNode(new Operation::Divide(), children[0].getChildren()[0], children[1].getChildren()[0]);
            if (children[1].getOperation().getId() == Operation::NEGATE && children[0].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Negate the constant
                return ExpressionTreeNode(new Operation::Divide(), ExpressionTreeNode(new Operation::MultiplyConstant(-dynamic_cast<const Operation::MultiplyConstant*>(&children[0].getOperation())->getValue()), children[0].getChildren()[0]), children[1].getChildren()[0]);
            if (children[0].getOperation().getId() == Operation::NEGATE) // Pull the negation out so it can possibly be optimized further
                return ExpressionTreeNode(new Operation::Negate(), ExpressionTreeNode(new Operation::Divide(), children[0].getChildren()[0], children[1]));
            if (children[1].getOperation().getId() == Operation::NEGATE) // Pull the negation out so it can possibly be optimized further
                return ExpressionTreeNode(new Operation::Negate(), ExpressionTreeNode(new Operation::Divide(), children[0], children[1].getChildren()[0]));
            if (children[1].getOperation().getId() == Operation::RECIPROCAL) // a/(1/b) = a*b
                return ExpressionTreeNode(new Operation::Multiply(), children[0], children[1].getChildren()[0]);
            break;
        }
        case Operation::POWER:
        {
            if (first_const && first == 0.0) // 0 to any power is 0
                return ExpressionTreeNode(new Operation::Constant(0.0));
            if (first_const && first == 1.0) // 1 to any power is 1
                return ExpressionTreeNode(new Operation::Constant(1.0));
            if (second_const) { // Constant exponent
                if (second == 0.0) // x^0 = 1
                    return ExpressionTreeNode(new Operation::Constant(1.0));
                if (second == 1.0) // x^1 = x
                    return children[0];
                if (second == -1.0) // x^-1 = recip(x)
                    return ExpressionTreeNode(new Operation::Reciprocal(), children[0]);
                if (second == 2.0) // x^2 = square(x)
                    return ExpressionTreeNode(new Operation::Square(), children[0]);
                if (second == 3.0) // x^3 = cube(x)
                    return ExpressionTreeNode(new Operation::Cube(), children[0]);
                if (second == 0.5) // x^0.5 = sqrt(x)
                    return ExpressionTreeNode(new Operation::Sqrt(), children[0]);
                // Constant power
                return ExpressionTreeNode(new Operation::PowerConstant(second), children[0]);
            }
            break;
        }
        case Operation::NEGATE:
        {
            if (children[0].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Combine a multiply and a negate into a single multiply
                return ExpressionTreeNode(new Operation::MultiplyConstant(-dynamic_cast<const Operation::MultiplyConstant*>(&children[0].getOperation())->getValue()), children[0].getChildren()[0]);
            if (first_const) // Negate a constant
                return ExpressionTreeNode(new Operation::Constant(-first));
            if (children[0].getOperation().getId() == Operation::NEGATE) // The two negations cancel
                return children[0].getChildren()[0];
            break;
        }
        case Operation::MULTIPLY_CONSTANT:
        {
            if (children[0].getOperation().getId() == Operation::MULTIPLY_CONSTANT) // Combine two multiplies into a single one
                return ExpressionTreeNode(new Operation::MultiplyConstant(dynamic_cast<const Operation::MultiplyConstant*>(&node.getOperation())->getValue()*dynamic_cast<const Operation::MultiplyConstant*>(&children[0].getOperation())->getValue()), children[0].getChildren()[0]);
            if (first_const) // Multiply two constants
                return ExpressionTreeNode(new Operation::Constant(dynamic_cast<const Operation::MultiplyConstant*>(&node.getOperation())->getValue()*getConstantValue(children[0])));
            if (children[0].getOperation().getId() == Operation::NEGATE) // Combine a multiply and a negate into a single multiply
                return ExpressionTreeNode(new Operation::MultiplyConstant(-dynamic_cast<const Operation::MultiplyConstant*>(&node.getOperation())->getValue()), children[0].getChildren()[0]);
            break;
        }
        case Operation::SQRT:
        {
            if (children[0].getOperation().getId() == Operation::SQUARE) // sqrt(square(x)) = abs(x)
                return ExpressionTreeNode(new Operation::Abs(), children[0].getChildren()[0]);
            break;
        }
        case Operation::SQUARE:
        {
            if (children[0].getOperation().getId() == Operation::SQRT) // square(sqrt(x)) = x
                return children[0].getChildren()[0];
            break;
        }
        case Operation::SELECT:
        {
            if (children[1] == children[2]) // Select between two identical values
                return children[1];
            break;
        }
        default:
        {
            // If operation ID is not one of the above,
            // we don't substitute a simpler expression.
            break;
        }

    }
    return ExpressionTreeNode(node.getOperation().clone(), children);
}

ParsedExpression ParsedExpression::differentiate(const string& variable) const {
    vector<const ExpressionTreeNode*> examples;
    getRootNode().assignTags(examples);
    map<int, ExpressionTreeNode> nodeCache;
    return differentiate(getRootNode(), variable, nodeCache);
}

ExpressionTreeNode ParsedExpression::differentiate(const ExpressionTreeNode& node, const string& variable, map<int, ExpressionTreeNode>& nodeCache) {
    auto cached = nodeCache.find(node.tag);
    if (cached != nodeCache.end())
        return cached->second;
    vector<ExpressionTreeNode> childDerivs(node.getChildren().size());
    for (int i = 0; i < (int) childDerivs.size(); i++)
        childDerivs[i] = differentiate(node.getChildren()[i], variable, nodeCache);
    ExpressionTreeNode result = node.getOperation().differentiate(node.getChildren(), childDerivs, variable);
    nodeCache[node.tag] = result;
    return result;
}

bool ParsedExpression::isConstant(const ExpressionTreeNode& node) {
    return (node.getOperation().getId() == Operation::CONSTANT);
}

double ParsedExpression::getConstantValue(const ExpressionTreeNode& node) {
    if (node.getOperation().getId() != Operation::CONSTANT) {
        throw Exception("getConstantValue called on a non-constant ExpressionNode");
    }
    return dynamic_cast<const Operation::Constant&>(node.getOperation()).getValue();
}

ExpressionProgram ParsedExpression::createProgram() const {
    return ExpressionProgram(*this);
}

CompiledExpression ParsedExpression::createCompiledExpression() const {
    return CompiledExpression(*this);
}

CompiledVectorExpression ParsedExpression::createCompiledVectorExpression(int width) const {
    return CompiledVectorExpression(*this, width);
}

ParsedExpression ParsedExpression::renameVariables(const map<string, string>& replacements) const {
    return ParsedExpression(renameNodeVariables(getRootNode(), replacements));
}

ExpressionTreeNode ParsedExpression::renameNodeVariables(const ExpressionTreeNode& node, const map<string, string>& replacements) {
    if (node.getOperation().getId() == Operation::VARIABLE) {
        map<string, string>::const_iterator replace = replacements.find(node.getOperation().getName());
        if (replace != replacements.end())
            return ExpressionTreeNode(new Operation::Variable(replace->second));
    }
    vector<ExpressionTreeNode> children;
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        children.push_back(renameNodeVariables(node.getChildren()[i], replacements));
    return ExpressionTreeNode(node.getOperation().clone(), children);
}

ostream& Lepton::operator<<(ostream& out, const ExpressionTreeNode& node) {
    if (node.getOperation().isInfixOperator() && node.getChildren().size() == 2) {
        out << "(" << node.getChildren()[0] << ")" << node.getOperation().getName() << "(" << node.getChildren()[1] << ")";
    }
    else if (node.getOperation().isInfixOperator() && node.getChildren().size() == 1) {
        out << "(" << node.getChildren()[0] << ")" << node.getOperation().getName();
    }
    else {
        out << node.getOperation().getName();
        if (node.getChildren().size() > 0) {
            out << "(";
            for (int i = 0; i < (int) node.getChildren().size(); i++) {
                if (i > 0)
                    out << ", ";
                out << node.getChildren()[i];
            }
            out << ")";
        }
    }
    return out;
}

ostream& Lepton::operator<<(ostream& out, const ParsedExpression& exp) {
    out << exp.getRootNode();
    return out;
}
