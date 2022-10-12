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

#include "lepton/ExpressionProgram.h"
#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"

using namespace Lepton;
using namespace std;

ExpressionProgram::ExpressionProgram() : maxArgs(0), stackSize(0) {
}

ExpressionProgram::ExpressionProgram(const ParsedExpression& expression) : maxArgs(0), stackSize(0) {
    buildProgram(expression.getRootNode());
    int currentStackSize = 0;
    for (int i = 0; i < (int) operations.size(); i++) {
        int args = operations[i]->getNumArguments();
        if (args > maxArgs)
            maxArgs = args;
        currentStackSize += 1-args;
        if (currentStackSize > stackSize)
            stackSize = currentStackSize;
    }
}

ExpressionProgram::~ExpressionProgram() {
    for (int i = 0; i < (int) operations.size(); i++)
        delete operations[i];
}

ExpressionProgram::ExpressionProgram(const ExpressionProgram& program) {
    *this = program;
}

ExpressionProgram& ExpressionProgram::operator=(const ExpressionProgram& program) {
    maxArgs = program.maxArgs;
    stackSize = program.stackSize;
    operations.resize(program.operations.size());
    for (int i = 0; i < (int) operations.size(); i++)
        operations[i] = program.operations[i]->clone();
    return *this;
}

void ExpressionProgram::buildProgram(const ExpressionTreeNode& node) {
    for (int i = (int) node.getChildren().size()-1; i >= 0; i--)
        buildProgram(node.getChildren()[i]);
    operations.push_back(node.getOperation().clone());
}

int ExpressionProgram::getNumOperations() const {
    return (int) operations.size();
}

const Operation& ExpressionProgram::getOperation(int index) const {
    return *operations[index];
}

void ExpressionProgram::setOperation(int index, Operation* operation) {
    delete operations[index];
    operations[index] = operation;
}

int ExpressionProgram::getStackSize() const {
    return stackSize;
}

double ExpressionProgram::evaluate() const {
    return evaluate(map<string, double>());
}

double ExpressionProgram::evaluate(const std::map<std::string, double>& variables) const {
    vector<double> stack(stackSize+1);
    int stackPointer = stackSize;
    for (int i = 0; i < (int) operations.size(); i++) {
        int numArgs = operations[i]->getNumArguments();
        double result = operations[i]->evaluate(&stack[stackPointer], variables);
        stackPointer += numArgs-1;
        stack[stackPointer] = result;
    }
    return stack[stackSize-1];
}
