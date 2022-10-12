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
#endif
    // Make a list of all variables we will need to copy before evaluating the expression.

    variablesToCopy.clear();
    for (map<string, int>::const_iterator iter = variableIndices.begin(); iter != variableIndices.end(); ++iter) {
        map<string, double*>::iterator pointer = variablePointers.find(iter->first);
        if (pointer != variablePointers.end())
            variablesToCopy.push_back(make_pair(&workspace[iter->second], pointer->second));
    }
}

double CompiledExpression::evaluate() const {
    if (jitCode)
        return jitCode();
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
}

#ifdef LEPTON_USE_JIT
static double evaluateOperation(Operation* op, double* args) {
    static map<string, double> dummyVariables;
    return op->evaluate(args, dummyVariables);
}

void CompiledExpression::findPowerGroups(vector<vector<int> >& groups, vector<vector<int> >& groupPowers, vector<int>& stepGroup) {
    // Identify every step that raises an argument to an integer power.

    vector<int> stepPower(operation.size(), 0);
    vector<int> stepArg(operation.size(), -1);
    for (int step = 0; step < operation.size(); step++) {
        Operation& op = *operation[step];
        int power = 0;
        if (op.getId() == Operation::SQUARE)
            power = 2;
        else if (op.getId() == Operation::CUBE)
            power = 3;
        else if (op.getId() == Operation::POWER_CONSTANT) {
            double realPower = dynamic_cast<const Operation::PowerConstant*>(&op)->getValue();
            if (realPower == (int) realPower)
                power = (int) realPower;
        }
        if (power != 0) {
            stepPower[step] = power;
            stepArg[step] = arguments[step][0];
        }
    }

    // Find groups that operate on the same argument and whose powers have the same sign.

    stepGroup.resize(operation.size(), -1);
    for (int i = 0; i < operation.size(); i++) {
        if (stepGroup[i] != -1)
            continue;
        vector<int> group, power;
        for (int j = i; j < operation.size(); j++) {
            if (stepArg[i] == stepArg[j] && stepPower[i]*stepPower[j] > 0) {
                stepGroup[j] = groups.size();
                group.push_back(j);
                power.push_back(stepPower[j]);
            }
        }
        groups.push_back(group);
        groupPowers.push_back(power);
    }
}

#if defined(__ARM__) || defined(__ARM64__)
void CompiledExpression::generateJitCode() {
    CodeHolder code;
    code.init(runtime.environment());
    a64::Compiler c(&code);
    c.addFunc(FuncSignatureT<double>());
    vector<arm::Vec> workspaceVar(workspace.size());
    for (int i = 0; i < (int) workspaceVar.size(); i++)
        workspaceVar[i] = c.newVecD();
    arm::Gp argsPointer = c.newIntPtr();
    c.mov(argsPointer, imm(&argValues[0]));
    vector<vector<int> > groups, groupPowers;
    vector<int> stepGroup;
    findPowerGroups(groups, groupPowers, stepGroup);

    // Load the arguments into variables.

    for (set<string>::const_iterator iter = variableNames.begin(); iter != variableNames.end(); ++iter) {
        map<string, int>::iterator index = variableIndices.find(*iter);
        arm::Gp variablePointer = c.newIntPtr();
        c.mov(variablePointer, imm(&getVariableReference(index->first)));
        c.ldr(workspaceVar[index->second], arm::ptr(variablePointer, 0));
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
        else if (op.getId() == Operation::POWER_CONSTANT) {
            if (stepGroup[step] == -1)
                value = dynamic_cast<Operation::PowerConstant&>(op).getValue();
            else
                value = 1.0;
        }
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

    vector<arm::Vec> constantVar(constants.size());
    if (constants.size() > 0) {
        arm::Gp constantsPointer = c.newIntPtr();
        c.mov(constantsPointer, imm(&constants[0]));
        for (int i = 0; i < (int) constants.size(); i++) {
            constantVar[i] = c.newVecD();
            c.ldr(constantVar[i], arm::ptr(constantsPointer, 8*i));
        }
    }

    // Evaluate the operations.

    vector<bool> hasComputedPower(operation.size(), false);
    for (int step = 0; step < (int) operation.size(); step++) {
        if (hasComputedPower[step])
            continue;

        // When one or more steps involve raising the same argument to multiple integer
        // powers, we can compute them all together for efficiency.

        if (stepGroup[step] != -1) {
            vector<int>& group = groups[stepGroup[step]];
            vector<int>& powers = groupPowers[stepGroup[step]];
            arm::Vec multiplier = c.newVecD();
            if (powers[0] > 0)
                c.fmov(multiplier, workspaceVar[arguments[step][0]]);
            else {
                c.fdiv(multiplier, constantVar[operationConstantIndex[step]], workspaceVar[arguments[step][0]]);
                for (int i = 0; i < powers.size(); i++)
                    powers[i] = -powers[i];
            }
            vector<bool> hasAssigned(group.size(), false);
            bool done = false;
            while (!done) {
                done = true;
                for (int i = 0; i < group.size(); i++) {
                    if (powers[i]%2 == 1) {
                        if (!hasAssigned[i])
                            c.fmov(workspaceVar[target[group[i]]], multiplier);
                        else
                            c.fmul(workspaceVar[target[group[i]]], workspaceVar[target[group[i]]], multiplier);
                        hasAssigned[i] = true;
                    }
                    powers[i] >>= 1;
                    if (powers[i] != 0)
                        done = false;
                }
                if (!done)
                    c.fmul(multiplier, multiplier, multiplier);
            }
            for (int step : group)
                hasComputedPower[step] = true;
            continue;
        }

        // Evaluate the step.

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
                c.fmov(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::ADD:
                c.fadd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::SUBTRACT:
                c.fsub(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::MULTIPLY:
                c.fmul(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::DIVIDE:
                c.fdiv(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::POWER:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], pow);
                break;
            case Operation::NEGATE:
                c.fneg(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::SQRT:
                c.fsqrt(workspaceVar[target[step]], workspaceVar[args[0]]);
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
                c.cmge(workspaceVar[target[step]], workspaceVar[args[0]], imm(0));
                c.and_(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::DELTA:
                c.cmeq(workspaceVar[target[step]], workspaceVar[args[0]], imm(0));
                c.and_(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::SQUARE:
                c.fmul(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
                break;
            case Operation::CUBE:
                c.fmul(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
                c.fmul(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::RECIPROCAL:
                c.fdiv(workspaceVar[target[step]], constantVar[operationConstantIndex[step]], workspaceVar[args[0]]);
                break;
            case Operation::ADD_CONSTANT:
                c.fadd(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::MULTIPLY_CONSTANT:
                c.fmul(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::POWER_CONSTANT:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]], pow);
                break;
            case Operation::MIN:
                c.fmin(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::MAX:
                c.fmax(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::ABS:
                c.fabs(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::FLOOR:
                c.frintm(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::CEIL:
                c.frintp(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::SELECT:
                c.fcmeq(workspaceVar[target[step]], workspaceVar[args[0]], imm(0));
                c.bsl(workspaceVar[target[step]], workspaceVar[args[2]], workspaceVar[args[1]]);
                break;
            default:
                // Just invoke evaluateOperation().

                for (int i = 0; i < (int) args.size(); i++)
                    c.str(workspaceVar[args[i]], arm::ptr(argsPointer, 8*i));
                arm::Gp fn = c.newIntPtr();
                c.mov(fn, imm((void*) evaluateOperation));
                InvokeNode* invoke;
                c.invoke(&invoke, fn, FuncSignatureT<double, Operation*, double*>());
                invoke->setArg(0, imm(&op));
                invoke->setArg(1, imm(&argValues[0]));
                invoke->setRet(0, workspaceVar[target[step]]);
        }
    }
    c.ret(workspaceVar[workspace.size()-1]);
    c.endFunc();
    c.finalize();
    runtime.add(&jitCode, &code);
}

void CompiledExpression::generateSingleArgCall(a64::Compiler& c, arm::Vec& dest, arm::Vec& arg, double (*function)(double)) {
    arm::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    InvokeNode* invoke;
    c.invoke(&invoke, fn, FuncSignatureT<double, double>());
    invoke->setArg(0, arg);
    invoke->setRet(0, dest);
}

void CompiledExpression::generateTwoArgCall(a64::Compiler& c, arm::Vec& dest, arm::Vec& arg1, arm::Vec& arg2, double (*function)(double, double)) {
    arm::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    InvokeNode* invoke;
    c.invoke(&invoke, fn, FuncSignatureT<double, double, double>());
    invoke->setArg(0, arg1);
    invoke->setArg(1, arg2);
    invoke->setRet(0, dest);
}
#else
void CompiledExpression::generateJitCode() {
    const CpuInfo& cpu = CpuInfo::host();
    if (!cpu.hasFeature(CpuFeatures::X86::kAVX))
        return;
    CodeHolder code;
    code.init(runtime.environment());
    x86::Compiler c(&code);
    FuncNode* funcNode = c.addFunc(FuncSignatureT<double>());
    funcNode->frame().setAvxEnabled();
    vector<x86::Xmm> workspaceVar(workspace.size());
    for (int i = 0; i < (int) workspaceVar.size(); i++)
        workspaceVar[i] = c.newXmmSd();
    x86::Gp argsPointer = c.newIntPtr();
    c.mov(argsPointer, imm(&argValues[0]));
    vector<vector<int> > groups, groupPowers;
    vector<int> stepGroup;
    findPowerGroups(groups, groupPowers, stepGroup);

    // Load the arguments into variables.

    x86::Gp variablePointer = c.newIntPtr();
    for (set<string>::const_iterator iter = variableNames.begin(); iter != variableNames.end(); ++iter) {
        map<string, int>::iterator index = variableIndices.find(*iter);
        c.mov(variablePointer, imm(&getVariableReference(index->first)));
        c.vmovsd(workspaceVar[index->second], x86::ptr(variablePointer, 0, 0));
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
        else if (op.getId() == Operation::ABS) {
            long long mask = 0x7FFFFFFFFFFFFFFF;
            value = *reinterpret_cast<double*>(&mask);
        }
        else if (op.getId() == Operation::POWER_CONSTANT) {
            if (stepGroup[step] == -1)
                value = dynamic_cast<Operation::PowerConstant&>(op).getValue();
            else
                value = 1.0;
        }
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

    vector<x86::Xmm> constantVar(constants.size());
    if (constants.size() > 0) {
        x86::Gp constantsPointer = c.newIntPtr();
        c.mov(constantsPointer, imm(&constants[0]));
        for (int i = 0; i < (int) constants.size(); i++) {
            constantVar[i] = c.newXmmSd();
            c.vmovsd(constantVar[i], x86::ptr(constantsPointer, 8*i, 0));
        }
    }

    // Evaluate the operations.

    vector<bool> hasComputedPower(operation.size(), false);
    for (int step = 0; step < (int) operation.size(); step++) {
        if (hasComputedPower[step])
            continue;

        // When one or more steps involve raising the same argument to multiple integer
        // powers, we can compute them all together for efficiency.

        if (stepGroup[step] != -1) {
            vector<int>& group = groups[stepGroup[step]];
            vector<int>& powers = groupPowers[stepGroup[step]];
            x86::Xmm multiplier = c.newXmmSd();
            if (powers[0] > 0)
                c.vmovsd(multiplier, workspaceVar[arguments[step][0]], workspaceVar[arguments[step][0]]);
            else {
                c.vdivsd(multiplier, constantVar[operationConstantIndex[step]], workspaceVar[arguments[step][0]]);
                for (int i = 0; i < powers.size(); i++)
                    powers[i] = -powers[i];
            }
            vector<bool> hasAssigned(group.size(), false);
            bool done = false;
            while (!done) {
                done = true;
                for (int i = 0; i < group.size(); i++) {
                    if (powers[i]%2 == 1) {
                        if (!hasAssigned[i])
                            c.vmovsd(workspaceVar[target[group[i]]], multiplier, multiplier);
                        else
                            c.vmulsd(workspaceVar[target[group[i]]], workspaceVar[target[group[i]]], multiplier);
                        hasAssigned[i] = true;
                    }
                    powers[i] >>= 1;
                    if (powers[i] != 0)
                        done = false;
                }
                if (!done)
                    c.vmulsd(multiplier, multiplier, multiplier);
            }
            for (int step : group)
                hasComputedPower[step] = true;
            continue;
        }

        // Evaluate the step.

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
                c.vmovsd(workspaceVar[target[step]], constantVar[operationConstantIndex[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::ADD:
                c.vaddsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::SUBTRACT:
                c.vsubsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::MULTIPLY:
                c.vmulsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::DIVIDE:
                c.vdivsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::POWER:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], pow);
                break;
            case Operation::NEGATE:
                c.vxorps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[target[step]]);
                c.vsubsd(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::SQRT:
                c.vsqrtsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
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
                c.vxorps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[target[step]]);
                c.vcmpsd(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]], imm(18)); // Comparison mode is _CMP_LE_OQ = 18
                c.vandps(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::DELTA:
                c.vxorps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[target[step]]);
                c.vcmpsd(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]], imm(16)); // Comparison mode is _CMP_EQ_OS = 16
                c.vandps(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::SQUARE:
                c.vmulsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
                break;
            case Operation::CUBE:
                c.vmulsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
                c.vmulsd(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::RECIPROCAL:
                c.vdivsd(workspaceVar[target[step]], constantVar[operationConstantIndex[step]], workspaceVar[args[0]]);
                break;
            case Operation::ADD_CONSTANT:
                c.vaddsd(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::MULTIPLY_CONSTANT:
                c.vmulsd(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::POWER_CONSTANT:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]], pow);
                break;
            case Operation::MIN:
                c.vminsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::MAX:
                c.vmaxsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::ABS:
                c.vandpd(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::FLOOR:
                c.vroundsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]], imm(1));
                break;
            case Operation::CEIL:
                c.vroundsd(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]], imm(2));
                break;
            case Operation::SELECT:
            {
                x86::Xmm mask = c.newXmmSd();
                c.vxorps(mask, mask, mask);
                c.vcmpsd(mask, mask, workspaceVar[args[0]], imm(0)); // Comparison mode is _CMP_EQ_OQ = 0
                c.vblendvps(workspaceVar[target[step]], workspaceVar[args[1]], workspaceVar[args[2]], mask);
                break;
            }
            default:
                // Just invoke evaluateOperation().

                for (int i = 0; i < (int) args.size(); i++)
                    c.vmovsd(x86::ptr(argsPointer, 8*i, 0), workspaceVar[args[i]]);
                x86::Gp fn = c.newIntPtr();
                c.mov(fn, imm((void*) evaluateOperation));
                InvokeNode* invoke;
                c.invoke(&invoke, fn, FuncSignatureT<double, Operation*, double*>());
                invoke->setArg(0, imm(&op));
                invoke->setArg(1, imm(&argValues[0]));
                invoke->setRet(0, workspaceVar[target[step]]);
        }
    }
    c.ret(workspaceVar[workspace.size()-1]);
    c.endFunc();
    c.finalize();
    runtime.add(&jitCode, &code);
}

void CompiledExpression::generateSingleArgCall(x86::Compiler& c, x86::Xmm& dest, x86::Xmm& arg, double (*function)(double)) {
    x86::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    InvokeNode* invoke;
    c.invoke(&invoke, fn, FuncSignatureT<double, double>());
    invoke->setArg(0, arg);
    invoke->setRet(0, dest);
}

void CompiledExpression::generateTwoArgCall(x86::Compiler& c, x86::Xmm& dest, x86::Xmm& arg1, x86::Xmm& arg2, double (*function)(double, double)) {
    x86::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    InvokeNode* invoke;
    c.invoke(&invoke, fn, FuncSignatureT<double, double, double>());
    invoke->setArg(0, arg1);
    invoke->setArg(1, arg2);
    invoke->setRet(0, dest);
}
#endif
#endif
