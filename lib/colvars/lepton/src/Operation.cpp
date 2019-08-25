
/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2015 Stanford University and the Authors.      *
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

#include "lepton/Operation.h"
#include "lepton/ExpressionTreeNode.h"
#include "MSVC_erfc.h"

using namespace Lepton;
using namespace std;

double Operation::Erf::evaluate(double* args, const map<string, double>& variables) const {
    return erf(args[0]);
}

double Operation::Erfc::evaluate(double* args, const map<string, double>& variables) const {
    return erfc(args[0]);
}

ExpressionTreeNode Operation::Constant::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Constant(0.0));
}

ExpressionTreeNode Operation::Variable::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    if (variable == name)
        return ExpressionTreeNode(new Operation::Constant(1.0));
    return ExpressionTreeNode(new Operation::Constant(0.0));
}

ExpressionTreeNode Operation::Custom::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    if (function->getNumArguments() == 0)
        return ExpressionTreeNode(new Operation::Constant(0.0));
    ExpressionTreeNode result = ExpressionTreeNode(new Operation::Multiply(), ExpressionTreeNode(new Operation::Custom(*this, 0), children), childDerivs[0]);
    for (int i = 1; i < getNumArguments(); i++) {
        result = ExpressionTreeNode(new Operation::Add(),
                                    result,
                                    ExpressionTreeNode(new Operation::Multiply(), ExpressionTreeNode(new Operation::Custom(*this, i), children), childDerivs[i]));
    }
    return result;
}

ExpressionTreeNode Operation::Add::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Add(), childDerivs[0], childDerivs[1]);
}

ExpressionTreeNode Operation::Subtract::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Subtract(), childDerivs[0], childDerivs[1]);
}

ExpressionTreeNode Operation::Multiply::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Add(),
                              ExpressionTreeNode(new Operation::Multiply(), children[0], childDerivs[1]),
                              ExpressionTreeNode(new Operation::Multiply(), children[1], childDerivs[0]));
}

ExpressionTreeNode Operation::Divide::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Divide(),
                              ExpressionTreeNode(new Operation::Subtract(),
                                                 ExpressionTreeNode(new Operation::Multiply(), children[1], childDerivs[0]),
                                                 ExpressionTreeNode(new Operation::Multiply(), children[0], childDerivs[1])),
                              ExpressionTreeNode(new Operation::Square(), children[1]));
}

ExpressionTreeNode Operation::Power::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Add(),
                              ExpressionTreeNode(new Operation::Multiply(),
                                                 ExpressionTreeNode(new Operation::Multiply(),
                                                                    children[1],
                                                                    ExpressionTreeNode(new Operation::Power(),
                                                                                       children[0], ExpressionTreeNode(new Operation::AddConstant(-1.0), children[1]))),
                                                 childDerivs[0]),
                              ExpressionTreeNode(new Operation::Multiply(),
                                                 ExpressionTreeNode(new Operation::Multiply(),
                                                                    ExpressionTreeNode(new Operation::Log(), children[0]),
                                                                    ExpressionTreeNode(new Operation::Power(), children[0], children[1])),
                                                 childDerivs[1]));
}

ExpressionTreeNode Operation::Negate::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Negate(), childDerivs[0]);
}

ExpressionTreeNode Operation::Sqrt::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::MultiplyConstant(0.5),
                                                 ExpressionTreeNode(new Operation::Reciprocal(),
                                                                    ExpressionTreeNode(new Operation::Sqrt(), children[0]))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Exp::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Exp(), children[0]),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Log::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Reciprocal(), children[0]),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Sin::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Cos(), children[0]),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Cos::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Negate(),
                                                 ExpressionTreeNode(new Operation::Sin(), children[0])),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Sec::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Multiply(),
                                                 ExpressionTreeNode(new Operation::Sec(), children[0]),
                                                 ExpressionTreeNode(new Operation::Tan(), children[0])),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Csc::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Negate(),
                                                 ExpressionTreeNode(new Operation::Multiply(),
                                                                    ExpressionTreeNode(new Operation::Csc(), children[0]),
                                                                    ExpressionTreeNode(new Operation::Cot(), children[0]))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Tan::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Square(),
                                                 ExpressionTreeNode(new Operation::Sec(), children[0])),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Cot::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Negate(),
                                                 ExpressionTreeNode(new Operation::Square(),
                                                                    ExpressionTreeNode(new Operation::Csc(), children[0]))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Asin::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Reciprocal(),
                                                 ExpressionTreeNode(new Operation::Sqrt(),
                                                                    ExpressionTreeNode(new Operation::Subtract(),
                                                                                       ExpressionTreeNode(new Operation::Constant(1.0)),
                                                                                       ExpressionTreeNode(new Operation::Square(), children[0])))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Acos::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Negate(),
                                                 ExpressionTreeNode(new Operation::Reciprocal(),
                                                                    ExpressionTreeNode(new Operation::Sqrt(),
                                                                                       ExpressionTreeNode(new Operation::Subtract(),
                                                                                                          ExpressionTreeNode(new Operation::Constant(1.0)),
                                                                                                          ExpressionTreeNode(new Operation::Square(), children[0]))))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Atan::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Reciprocal(),
                                                 ExpressionTreeNode(new Operation::AddConstant(1.0),
                                                                    ExpressionTreeNode(new Operation::Square(), children[0]))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Sinh::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Cosh(),
                                                 children[0]),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Cosh::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Sinh(),
                                                 children[0]),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Tanh::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Subtract(),
                                                 ExpressionTreeNode(new Operation::Constant(1.0)),
                                                 ExpressionTreeNode(new Operation::Square(),
                                                                    ExpressionTreeNode(new Operation::Tanh(), children[0]))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Erf::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Multiply(),
                                                 ExpressionTreeNode(new Operation::Constant(2.0/sqrt(M_PI))),
                                                 ExpressionTreeNode(new Operation::Exp(),
                                                                    ExpressionTreeNode(new Operation::Negate(),
                                                                                       ExpressionTreeNode(new Operation::Square(), children[0])))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Erfc::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Multiply(),
                                                 ExpressionTreeNode(new Operation::Constant(-2.0/sqrt(M_PI))),
                                                 ExpressionTreeNode(new Operation::Exp(),
                                                                    ExpressionTreeNode(new Operation::Negate(),
                                                                                       ExpressionTreeNode(new Operation::Square(), children[0])))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Step::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Constant(0.0));
}

ExpressionTreeNode Operation::Delta::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Constant(0.0));
}

ExpressionTreeNode Operation::Square::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::MultiplyConstant(2.0),
                                                 children[0]),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Cube::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::MultiplyConstant(3.0),
                                                 ExpressionTreeNode(new Operation::Square(), children[0])),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Reciprocal::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::Negate(),
                                                 ExpressionTreeNode(new Operation::Reciprocal(),
                                                                    ExpressionTreeNode(new Operation::Square(), children[0]))),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::AddConstant::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return childDerivs[0];
}

ExpressionTreeNode Operation::MultiplyConstant::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::MultiplyConstant(value),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::PowerConstant::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Multiply(),
                              ExpressionTreeNode(new Operation::MultiplyConstant(value),
                                                 ExpressionTreeNode(new Operation::PowerConstant(value-1),
                                                                    children[0])),
                              childDerivs[0]);
}

ExpressionTreeNode Operation::Min::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    ExpressionTreeNode step(new Operation::Step(),
                            ExpressionTreeNode(new Operation::Subtract(), children[0], children[1]));
    return ExpressionTreeNode(new Operation::Subtract(),
                              ExpressionTreeNode(new Operation::Multiply(), childDerivs[1], step),
                              ExpressionTreeNode(new Operation::Multiply(), childDerivs[0],
                                                 ExpressionTreeNode(new Operation::AddConstant(-1), step)));
}

ExpressionTreeNode Operation::Max::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    ExpressionTreeNode step(new Operation::Step(),
                            ExpressionTreeNode(new Operation::Subtract(), children[0], children[1]));
    return ExpressionTreeNode(new Operation::Subtract(),
                              ExpressionTreeNode(new Operation::Multiply(), childDerivs[0], step),
                              ExpressionTreeNode(new Operation::Multiply(), childDerivs[1],
                                                 ExpressionTreeNode(new Operation::AddConstant(-1), step)));
}

ExpressionTreeNode Operation::Abs::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    ExpressionTreeNode step(new Operation::Step(), children[0]);
    return ExpressionTreeNode(new Operation::Multiply(),
                              childDerivs[0],
                              ExpressionTreeNode(new Operation::AddConstant(-1),
                                                 ExpressionTreeNode(new Operation::MultiplyConstant(2), step)));
}

ExpressionTreeNode Operation::Floor::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Constant(0.0));
}

ExpressionTreeNode Operation::Ceil::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    return ExpressionTreeNode(new Operation::Constant(0.0));
}

ExpressionTreeNode Operation::Select::differentiate(const std::vector<ExpressionTreeNode>& children, const std::vector<ExpressionTreeNode>& childDerivs, const std::string& variable) const {
    vector<ExpressionTreeNode> derivChildren;
    derivChildren.push_back(children[0]);
    derivChildren.push_back(childDerivs[1]);
    derivChildren.push_back(childDerivs[2]);
    return ExpressionTreeNode(new Operation::Select(), derivChildren);
}
