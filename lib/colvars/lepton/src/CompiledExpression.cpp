/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2019 Stanford University and the Authors.      *
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

#include "lepton/CompiledExpression.h"
#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"
#include <utility>

using namespace Lepton;
using namespace std;
#ifdef LEPTON_USE_JIT
    using namespace asmjit;
#endif

CompiledExpression::CompiledExpression() : jitCode(NULL) {
}

CompiledExpression::CompiledExpression(const ParsedExpression& expression) : jitCode(NULL) {
    ParsedExpression expr = expression.optimize(); // Just in case it wasn't already optimized.
    vector<pair<ExpressionTreeNode, int> > temps;
    compileExpression(expr.getRootNode(), temps);
    int maxArguments = 1;
    for (int i = 0; i < (int) operation.size(); i++)
        if (operation[i]->getNumArguments() > maxArguments)
            maxArguments = operation[i]->getNumArguments();
    argValues.resize(maxArguments);
#ifdef LEPTON_USE_JIT
    generateJitCode();
#endif
}

CompiledExpression::~CompiledExpression() {
    for (int i = 0; i < (int) operation.size(); i++)
        if (operation[i] != NULL)
            delete operation[i];
}

CompiledExpression::CompiledExpression(const CompiledExpression& expression) : jitCode(NULL) {
    *this = expression;
}

CompiledExpression& CompiledExpression::operator=(const CompiledExpression& expression) {
    arguments = expression.arguments;
    target = expression.target;
    variableIndices = expression.variableIndices;
    variableNames = expression.variableNames;
    workspace.resize(expression.workspace.size());
    argValues.resize(expression.argValues.size());
    operation.resize(expression.operation.size());
    for (int i = 0; i < (int) operation.size(); i++)
        operation[i] = expression.operation[i]->clone();
    setVariableLocations(variablePointers);
    return *this;
}

void CompiledExpression::compileExpression(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, int> >& temps) {
    if (findTempIndex(node, temps) != -1)
        return; // We have already processed a node identical to this one.

    // Process the child nodes.

    vector<int> args;
    for (int i = 0; i < node.getChildren().size(); i++) {
        compileExpression(node.getChildren()[i], temps);
        args.push_back(findTempIndex(node.getChildren()[i], temps));
    }

    // Process this node.

    if (node.getOperation().getId() == Operation::VARIABLE) {
        variableIndices[node.getOperation().getName()] = (int) workspace.size();
        variableNames.insert(node.getOperation().getName());
    }
    else {
        int stepIndex = (int) arguments.size();
        arguments.push_back(vector<int>());
        target.push_back((int) workspace.size());
        operation.push_back(node.getOperation().clone());
        if (args.size() == 0)
            arguments[stepIndex].push_back(0); // The value won't actually be used.  We just need something there.
        else {
            // If the arguments are sequential, we can just pass a pointer to the first one.

            bool sequential = true;
            for (int i = 1; i < args.size(); i++)
                if (args[i] != args[i-1]+1)
                    sequential = false;
            if (sequential)
                arguments[stepIndex].push_back(args[0]);
            else
                arguments[stepIndex] = args;
        }
    }
    temps.push_back(make_pair(node, (int) workspace.size()));
    workspace.push_back(0.0);
}

int CompiledExpression::findTempIndex(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, int> >& temps) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return i;
    return -1;
}

const set<string>& CompiledExpression::getVariables() const {
    return variableNames;
}

double& CompiledExpression::getVariableReference(const string& name) {
    map<string, double*>::iterator pointer = variablePointers.find(name);
    if (pointer != variablePointers.end())
        return *pointer->second;
    map<string, int>::iterator index = variableIndices.find(name);
    if (index == variableIndices.end())
        throw Exception("getVariableReference: Unknown variable '"+name+"'");
    return workspace[index->second];
}

void CompiledExpression::setVariableLocations(map<string, double*>& variableLocations) {
    variablePointers = variableLocations;
#ifdef LEPTON_USE_JIT
    // Rebuild the JIT code.

    if (workspace.size() > 0)
        generateJitCode();
#else
    // Make a list of all variables we will need to copy before evaluating the expression.

    variablesToCopy.clear();
    for (map<string, int>::const_iterator iter = variableIndices.begin(); iter != variableIndices.end(); ++iter) {
        map<string, double*>::iterator pointer = variablePointers.find(iter->first);
        if (pointer != variablePointers.end())
            variablesToCopy.push_back(make_pair(&workspace[iter->second], pointer->second));
    }
#endif
}

double CompiledExpression::evaluate() const {
#ifdef LEPTON_USE_JIT
    return jitCode();
#else
    for (int i = 0; i < variablesToCopy.size(); i++)
        *variablesToCopy[i].first = *variablesToCopy[i].second;

    // Loop over the operations and evaluate each one.

    for (int step = 0; step < operation.size(); step++) {
        const vector<int>& args = arguments[step];
        if (args.size() == 1)
            workspace[target[step]] = operation[step]->evaluate(&workspace[args[0]], dummyVariables);
        else {
            for (int i = 0; i < args.size(); i++)
                argValues[i] = workspace[args[i]];
            workspace[target[step]] = operation[step]->evaluate(&argValues[0], dummyVariables);
        }
    }
    return workspace[workspace.size()-1];
#endif
}

#ifdef LEPTON_USE_JIT
static double evaluateOperation(Operation* op, double* args) {
    static map<string, double> dummyVariables;
    return op->evaluate(args, dummyVariables);
}

void CompiledExpression::generateJitCode() {
    CodeHolder code;
    code.init(runtime.getCodeInfo());
    X86Compiler c(&code);
    c.addFunc(FuncSignature0<double>());
    vector<X86Xmm> workspaceVar(workspace.size());
    for (int i = 0; i < (int) workspaceVar.size(); i++)
        workspaceVar[i] = c.newXmmSd();
    X86Gp argsPointer = c.newIntPtr();
    c.mov(argsPointer, imm_ptr(&argValues[0]));

    // Load the arguments into variables.

    for (set<string>::const_iterator iter = variableNames.begin(); iter != variableNames.end(); ++iter) {
        map<string, int>::iterator index = variableIndices.find(*iter);
        X86Gp variablePointer = c.newIntPtr();
        c.mov(variablePointer, imm_ptr(&getVariableReference(index->first)));
        c.movsd(workspaceVar[index->second], x86::ptr(variablePointer, 0, 0));
    }

    // Make a list of all constants that will be needed for evaluation.

    vector<int> operationConstantIndex(operation.size(), -1);
    for (int step = 0; step < (int) operation.size(); step++) {
        // Find the constant value (if any) used by this operation.

        Operation& op = *operation[step];
        double value;
        if (op.getId() == Operation::CONSTANT)
            value = dynamic_cast<Operation::Constant&>(op).getValue();
        else if (op.getId() == Operation::ADD_CONSTANT)
            value = dynamic_cast<Operation::AddConstant&>(op).getValue();
        else if (op.getId() == Operation::MULTIPLY_CONSTANT)
            value = dynamic_cast<Operation::MultiplyConstant&>(op).getValue();
        else if (op.getId() == Operation::RECIPROCAL)
            value = 1.0;
        else if (op.getId() == Operation::STEP)
            value = 1.0;
        else if (op.getId() == Operation::DELTA)
            value = 1.0;
        else
            continue;

        // See if we already have a variable for this constant.

        for (int i = 0; i < (int) constants.size(); i++)
            if (value == constants[i]) {
                operationConstantIndex[step] = i;
                break;
            }
        if (operationConstantIndex[step] == -1) {
            operationConstantIndex[step] = constants.size();
            constants.push_back(value);
        }
    }

    // Load constants into variables.

    vector<X86Xmm> constantVar(constants.size());
    if (constants.size() > 0) {
        X86Gp constantsPointer = c.newIntPtr();
        c.mov(constantsPointer, imm_ptr(&constants[0]));
        for (int i = 0; i < (int) constants.size(); i++) {
            constantVar[i] = c.newXmmSd();
            c.movsd(constantVar[i], x86::ptr(constantsPointer, 8*i, 0));
        }
    }

    // Evaluate the operations.

    for (int step = 0; step < (int) operation.size(); step++) {
        Operation& op = *operation[step];
        vector<int> args = arguments[step];
        if (args.size() == 1) {
            // One or more sequential arguments.  Fill out the list.

            for (int i = 1; i < op.getNumArguments(); i++)
                args.push_back(args[0]+i);
        }

        // Generate instructions to execute this operation.

        switch (op.getId()) {
            case Operation::CONSTANT:
                c.movsd(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::ADD:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.addsd(workspaceVar[target[step]], workspaceVar[args[1]]);
                break;
            case Operation::SUBTRACT:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.subsd(workspaceVar[target[step]], workspaceVar[args[1]]);
                break;
            case Operation::MULTIPLY:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.mulsd(workspaceVar[target[step]], workspaceVar[args[1]]);
                break;
            case Operation::DIVIDE:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.divsd(workspaceVar[target[step]], workspaceVar[args[1]]);
                break;
            case Operation::POWER:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], pow);
                break;
            case Operation::NEGATE:
                c.xorps(workspaceVar[target[step]], workspaceVar[target[step]]);
                c.subsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::SQRT:
                c.sqrtsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::EXP:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], exp);
                break;
            case Operation::LOG:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], log);
                break;
            case Operation::SIN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], sin);
                break;
            case Operation::COS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], cos);
                break;
            case Operation::TAN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], tan);
                break;
            case Operation::ASIN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], asin);
                break;
            case Operation::ACOS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], acos);
                break;
            case Operation::ATAN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], atan);
                break;
            case Operation::ATAN2:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], atan2);
                break;
            case Operation::SINH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], sinh);
                break;
            case Operation::COSH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], cosh);
                break;
            case Operation::TANH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], tanh);
                break;
            case Operation::STEP:
                c.xorps(workspaceVar[target[step]], workspaceVar[target[step]]);
                c.cmpsd(workspaceVar[target[step]], workspaceVar[args[0]], imm(18)); // Comparison mode is _CMP_LE_OQ = 18
                c.andps(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::DELTA:
                c.xorps(workspaceVar[target[step]], workspaceVar[target[step]]);
                c.cmpsd(workspaceVar[target[step]], workspaceVar[args[0]], imm(16)); // Comparison mode is _CMP_EQ_OS = 16
                c.andps(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::SQUARE:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.mulsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::CUBE:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.mulsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.mulsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::RECIPROCAL:
                c.movsd(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                c.divsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::ADD_CONSTANT:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.addsd(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::MULTIPLY_CONSTANT:
                c.movsd(workspaceVar[target[step]], workspaceVar[args[0]]);
                c.mulsd(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::ABS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], fabs);
                break;
            case Operation::FLOOR:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], floor);
                break;
            case Operation::CEIL:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], ceil);
                break;
            default:
                // Just invoke evaluateOperation().

                for (int i = 0; i < (int) args.size(); i++)
                    c.movsd(x86::ptr(argsPointer, 8*i, 0), workspaceVar[args[i]]);
                X86Gp fn = c.newIntPtr();
                c.mov(fn, imm_ptr((void*) evaluateOperation));
                CCFuncCall* call = c.call(fn, FuncSignature2<double, Operation*, double*>());
                call->setArg(0, imm_ptr(&op));
                call->setArg(1, imm_ptr(&argValues[0]));
                call->setRet(0, workspaceVar[target[step]]);
        }
    }
    c.ret(workspaceVar[workspace.size()-1]);
    c.endFunc();
    c.finalize();
    runtime.add(&jitCode, &code);
}

void CompiledExpression::generateSingleArgCall(X86Compiler& c, X86Xmm& dest, X86Xmm& arg, double (*function)(double)) {
    X86Gp fn = c.newIntPtr();
    c.mov(fn, imm_ptr((void*) function));
    CCFuncCall* call = c.call(fn, FuncSignature1<double, double>());
    call->setArg(0, arg);
    call->setRet(0, dest);
}

void CompiledExpression::generateTwoArgCall(X86Compiler& c, X86Xmm& dest, X86Xmm& arg1, X86Xmm& arg2, double (*function)(double, double)) {
    X86Gp fn = c.newIntPtr();
    c.mov(fn, imm_ptr((void*) function));
    CCFuncCall* call = c.call(fn, FuncSignature2<double, double, double>());
    call->setArg(0, arg1);
    call->setArg(1, arg2);
    call->setRet(0, dest);
}
#endif
