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

#include "lepton/CompiledVectorExpression.h"
#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"
#include <algorithm>
#include <utility>

using namespace Lepton;
using namespace std;
#ifdef LEPTON_USE_JIT
using namespace asmjit;
#endif

CompiledVectorExpression::CompiledVectorExpression() : jitCode(NULL) {
}

CompiledVectorExpression::CompiledVectorExpression(const ParsedExpression& expression, int width) : width(width), jitCode(NULL) {
    const vector<int> allowedWidths = getAllowedWidths();
    if (find(allowedWidths.begin(), allowedWidths.end(), width) == allowedWidths.end())
        throw Exception("Unsupported width for vector expression: "+to_string(width));
    ParsedExpression expr = expression.optimize(); // Just in case it wasn't already optimized.
    vector<pair<ExpressionTreeNode, int> > temps;
    int workspaceSize = 0;
    compileExpression(expr.getRootNode(), temps, workspaceSize);
    workspace.resize(workspaceSize*width);
    int maxArguments = 1;
    for (int i = 0; i < (int) operation.size(); i++)
        if (operation[i]->getNumArguments() > maxArguments)
            maxArguments = operation[i]->getNumArguments();
    argValues.resize(maxArguments);
#ifdef LEPTON_USE_JIT
    generateJitCode();
#endif
}

CompiledVectorExpression::~CompiledVectorExpression() {
    for (int i = 0; i < (int) operation.size(); i++)
        if (operation[i] != NULL)
            delete operation[i];
}

CompiledVectorExpression::CompiledVectorExpression(const CompiledVectorExpression& expression) : jitCode(NULL) {
    *this = expression;
}

CompiledVectorExpression& CompiledVectorExpression::operator=(const CompiledVectorExpression& expression) {
    arguments = expression.arguments;
    width = expression.width;
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

const vector<int>& CompiledVectorExpression::getAllowedWidths() {
    static vector<int> widths;
    if (widths.size() == 0) {
        widths.push_back(4);
#ifdef LEPTON_USE_JIT
        const CpuInfo& cpu = CpuInfo::host();
        if (cpu.hasFeature(CpuFeatures::X86::kAVX))
            widths.push_back(8);
#endif
    }
    return widths;
}

void CompiledVectorExpression::compileExpression(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, int> >& temps, int& workspaceSize) {
    if (findTempIndex(node, temps) != -1)
        return; // We have already processed a node identical to this one.

    // Process the child nodes.

    vector<int> args;
    for (int i = 0; i < (int)node.getChildren().size(); i++) {
        compileExpression(node.getChildren()[i], temps, workspaceSize);
        args.push_back(findTempIndex(node.getChildren()[i], temps));
    }

    // Process this node.

    if (node.getOperation().getId() == Operation::VARIABLE) {
        variableIndices[node.getOperation().getName()] = workspaceSize;
        variableNames.insert(node.getOperation().getName());
    }
    else {
        int stepIndex = (int) arguments.size();
        arguments.push_back(vector<int>());
        target.push_back(workspaceSize);
        operation.push_back(node.getOperation().clone());
        if (args.size() == 0)
            arguments[stepIndex].push_back(0); // The value won't actually be used.  We just need something there.
        else {
            // If the arguments are sequential, we can just pass a pointer to the first one.

            bool sequential = true;
            for (int i = 1; i < (int)args.size(); i++)
                if (args[i] != args[i - 1] + 1)
                    sequential = false;
            if (sequential)
                arguments[stepIndex].push_back(args[0]);
            else
                arguments[stepIndex] = args;
        }
    }
    temps.push_back(make_pair(node, workspaceSize));
    workspaceSize++;
}

int CompiledVectorExpression::findTempIndex(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, int> >& temps) {
    for (int i = 0; i < (int) temps.size(); i++)
        if (temps[i].first == node)
            return i;
    return -1;
}

int CompiledVectorExpression::getWidth() const {
    return width;
}

const set<string>& CompiledVectorExpression::getVariables() const {
    return variableNames;
}

float* CompiledVectorExpression::getVariablePointer(const string& name) {
    map<string, float*>::iterator pointer = variablePointers.find(name);
    if (pointer != variablePointers.end())
        return pointer->second;
    map<string, int>::iterator index = variableIndices.find(name);
    if (index == variableIndices.end())
        throw Exception("getVariableReference: Unknown variable '" + name + "'");
    return &workspace[index->second*width];
}

void CompiledVectorExpression::setVariableLocations(map<string, float*>& variableLocations) {
    variablePointers = variableLocations;
#ifdef LEPTON_USE_JIT
    // Rebuild the JIT code.

    if (workspace.size() > 0)
        generateJitCode();
#endif
    // Make a list of all variables we will need to copy before evaluating the expression.

    variablesToCopy.clear();
    for (map<string, int>::const_iterator iter = variableIndices.begin(); iter != variableIndices.end(); ++iter) {
        map<string, float*>::iterator pointer = variablePointers.find(iter->first);
        if (pointer != variablePointers.end())
            variablesToCopy.push_back(make_pair(&workspace[iter->second*width], pointer->second));
    }
}

const float* CompiledVectorExpression::evaluate() const {
    if (jitCode) {
        jitCode();
        return &workspace[workspace.size()-width];
    }
    for (int i = 0; i < (int)variablesToCopy.size(); i++)
        for (int j = 0; j < width; j++)
            variablesToCopy[i].first[j] = variablesToCopy[i].second[j];

    // Loop over the operations and evaluate each one.

    for (int step = 0; step < (int)operation.size(); step++) {
        const vector<int>& args = arguments[step];
        if (args.size() == 1) {
            for (int j = 0; j < width; j++) {
                for (int i = 0; i < operation[step]->getNumArguments(); i++)
                    argValues[i] = workspace[(args[0]+i)*width+j];
                workspace[target[step]*width+j] = operation[step]->evaluate(&argValues[0], dummyVariables);
            }
        } else {
            for (int j = 0; j < width; j++) {
              for (int i = 0; i < (int)args.size(); i++)
                    argValues[i] = workspace[args[i]*width+j];
                workspace[target[step]*width+j] = operation[step]->evaluate(&argValues[0], dummyVariables);
            }
        }
    }
    return &workspace[workspace.size()-width];
}

#ifdef LEPTON_USE_JIT

static double evaluateOperation(Operation* op, double* args) {
    static map<string, double> dummyVariables;
    return op->evaluate(args, dummyVariables);
}

void CompiledVectorExpression::findPowerGroups(vector<vector<int> >& groups, vector<vector<int> >& groupPowers, vector<int>& stepGroup) {
    // Identify every step that raises an argument to an integer power.

    vector<int> stepPower(operation.size(), 0);
    vector<int> stepArg(operation.size(), -1);
    for (int step = 0; step < (int)operation.size(); step++) {
        Operation& op = *operation[step];
        int power = 0;
        if (op.getId() == Operation::SQUARE)
            power = 2;
        else if (op.getId() == Operation::CUBE)
            power = 3;
        else if (op.getId() == Operation::POWER_CONSTANT) {
            double realPower = dynamic_cast<const Operation::PowerConstant*> (&op)->getValue();
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
    for (int i = 0; i < (int)operation.size(); i++) {
        if (stepGroup[i] != -1)
            continue;
        vector<int> group, power;
        for (int j = i; j < (int)operation.size(); j++) {
            if (stepArg[i] == stepArg[j] && stepPower[i] * stepPower[j] > 0) {
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

void CompiledVectorExpression::generateJitCode() {
    CodeHolder code;
    code.init(runtime.environment());
    a64::Compiler c(&code);
    c.addFunc(FuncSignatureT<void>());
    vector<arm::Vec> workspaceVar(workspace.size()/width);
    for (int i = 0; i < (int) workspaceVar.size(); i++)
        workspaceVar[i] = c.newVecQ();
    arm::Gp argsPointer = c.newIntPtr();
    c.mov(argsPointer, imm(&argValues[0]));
    vector<vector<int> > groups, groupPowers;
    vector<int> stepGroup;
    findPowerGroups(groups, groupPowers, stepGroup);

    // Load the arguments into variables.

    arm::Gp variablePointer = c.newIntPtr();
    for (set<string>::const_iterator iter = variableNames.begin(); iter != variableNames.end(); ++iter) {
        map<string, int>::iterator index = variableIndices.find(*iter);
        c.mov(variablePointer, imm(getVariablePointer(index->first)));
        c.ldr(workspaceVar[index->second].s4(), arm::ptr(variablePointer, 0));
    }

    // Make a list of all constants that will be needed for evaluation.

    vector<int> operationConstantIndex(operation.size(), -1);
    for (int step = 0; step < (int) operation.size(); step++) {
        // Find the constant value (if any) used by this operation.

        Operation& op = *operation[step];
        float value;
        if (op.getId() == Operation::CONSTANT)
            value = dynamic_cast<Operation::Constant&> (op).getValue();
        else if (op.getId() == Operation::ADD_CONSTANT)
            value = dynamic_cast<Operation::AddConstant&> (op).getValue();
        else if (op.getId() == Operation::MULTIPLY_CONSTANT)
            value = dynamic_cast<Operation::MultiplyConstant&> (op).getValue();
        else if (op.getId() == Operation::RECIPROCAL)
            value = 1.0;
        else if (op.getId() == Operation::STEP)
            value = 1.0;
        else if (op.getId() == Operation::DELTA)
            value = 1.0;
        else if (op.getId() == Operation::POWER_CONSTANT) {
            if (stepGroup[step] == -1)
                value = dynamic_cast<Operation::PowerConstant&> (op).getValue();
            else
                value = 1.0;
        } else
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
        for (int i = 0; i < (int) constants.size(); i++) {
            c.mov(constantsPointer, imm(&constants[i]));
            constantVar[i] = c.newVecQ();
            c.ld1r(constantVar[i].s4(), arm::ptr(constantsPointer));
        }
    }

    // Evaluate the operations.

    vector<bool> hasComputedPower(operation.size(), false);
    arm::Vec argReg = c.newVecS();
    arm::Vec doubleArgReg = c.newVecD();
    arm::Vec doubleResultReg = c.newVecD();
    for (int step = 0; step < (int) operation.size(); step++) {
        if (hasComputedPower[step])
            continue;

        // When one or more steps involve raising the same argument to multiple integer
        // powers, we can compute them all together for efficiency.

        if (stepGroup[step] != -1) {
            vector<int>& group = groups[stepGroup[step]];
            vector<int>& powers = groupPowers[stepGroup[step]];
            arm::Vec multiplier = c.newVecQ();
            if (powers[0] > 0)
                c.mov(multiplier.s4(), workspaceVar[arguments[step][0]].s4());
            else {
                c.fdiv(multiplier.s4(), constantVar[operationConstantIndex[step]].s4(), workspaceVar[arguments[step][0]].s4());
                for (int i = 0; i < powers.size(); i++)
                    powers[i] = -powers[i];
            }
            vector<bool> hasAssigned(group.size(), false);
            bool done = false;
            while (!done) {
                done = true;
                for (int i = 0; i < group.size(); i++) {
                    if (powers[i] % 2 == 1) {
                        if (!hasAssigned[i])
                            c.mov(workspaceVar[target[group[i]]].s4(), multiplier.s4());
                        else
                            c.fmul(workspaceVar[target[group[i]]].s4(), workspaceVar[target[group[i]]].s4(), multiplier.s4());
                        hasAssigned[i] = true;
                    }
                    powers[i] >>= 1;
                    if (powers[i] != 0)
                        done = false;
                }
                if (!done)
                    c.fmul(multiplier.s4(), multiplier.s4(), multiplier.s4());
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
                args.push_back(args[0] + i);
        }

        // Generate instructions to execute this operation.

        switch (op.getId()) {
            case Operation::CONSTANT:
                c.mov(workspaceVar[target[step]].s4(), constantVar[operationConstantIndex[step]].s4());
                break;
            case Operation::ADD:
                c.fadd(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[1]].s4());
                break;
            case Operation::SUBTRACT:
                c.fsub(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[1]].s4());
                break;
            case Operation::MULTIPLY:
                c.fmul(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[1]].s4());
                break;
            case Operation::DIVIDE:
                c.fdiv(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[1]].s4());
                break;
            case Operation::POWER:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], powf);
                break;
            case Operation::NEGATE:
                c.fneg(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::SQRT:
                c.fsqrt(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::EXP:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], expf);
                break;
            case Operation::LOG:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], logf);
                break;
            case Operation::SIN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], sinf);
                break;
            case Operation::COS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], cosf);
                break;
            case Operation::TAN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], tanf);
                break;
            case Operation::ASIN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], asinf);
                break;
            case Operation::ACOS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], acosf);
                break;
            case Operation::ATAN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], atanf);
                break;
            case Operation::ATAN2:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], atan2f);
                break;
            case Operation::SINH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], sinhf);
                break;
            case Operation::COSH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], coshf);
                break;
            case Operation::TANH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], tanhf);
                break;
            case Operation::STEP:
                c.cmge(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), imm(0));
                c.and_(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::DELTA:
                c.cmeq(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), imm(0));
                c.and_(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::SQUARE:
                c.fmul(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::CUBE:
                c.fmul(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[0]].s4());
                c.fmul(workspaceVar[target[step]].s4(), workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::RECIPROCAL:
                c.fdiv(workspaceVar[target[step]].s4(), constantVar[operationConstantIndex[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::ADD_CONSTANT:
                c.fadd(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), constantVar[operationConstantIndex[step]].s4());
                break;
            case Operation::MULTIPLY_CONSTANT:
                c.fmul(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), constantVar[operationConstantIndex[step]].s4());
                break;
            case Operation::POWER_CONSTANT:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]], powf);
                break;
            case Operation::MIN:
                c.fmin(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[1]].s4());
                break;
            case Operation::MAX:
                c.fmax(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), workspaceVar[args[1]].s4());
                break;
            case Operation::ABS:
                c.fabs(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::FLOOR:
                c.frintm(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::CEIL:
                c.frintp(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4());
                break;
            case Operation::SELECT:
                c.fcmeq(workspaceVar[target[step]].s4(), workspaceVar[args[0]].s4(), imm(0));
                c.bsl(workspaceVar[target[step]], workspaceVar[args[2]], workspaceVar[args[1]]);
                break;
            default:
                // Just invoke evaluateOperation().
                for (int element = 0; element < width; element++) {
                    for (int i = 0; i < (int) args.size(); i++) {
                        c.ins(argReg.s(0), workspaceVar[args[i]].s(element));
                        c.fcvt(doubleArgReg, argReg);
                        c.str(doubleArgReg, arm::ptr(argsPointer, 8*i));
                    }
                    arm::Gp fn = c.newIntPtr();
                    c.mov(fn, imm((void*) evaluateOperation));
                    InvokeNode* invoke;
                    c.invoke(&invoke, fn, FuncSignatureT<double, Operation*, double*>());
                    invoke->setArg(0, imm(&op));
                    invoke->setArg(1, imm(&argValues[0]));
                    invoke->setRet(0, doubleResultReg);
                    c.fcvt(argReg, doubleResultReg);
                    c.ins(workspaceVar[target[step]].s(element), argReg.s(0));
                }
        }
    }
    arm::Gp resultPointer = c.newIntPtr();
    c.mov(resultPointer, imm(&workspace[workspace.size()-width]));
    c.str(workspaceVar.back().s4(), arm::ptr(resultPointer, 0));
    c.endFunc();
    c.finalize();
    runtime.add(&jitCode, &code);
}

void CompiledVectorExpression::generateSingleArgCall(a64::Compiler& c, arm::Vec& dest, arm::Vec& arg, float (*function)(float)) {
    arm::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    arm::Vec a = c.newVecS();
    arm::Vec d = c.newVecS();
    for (int element = 0; element < width; element++) {
        c.ins(a.s(0), arg.s(element));
        InvokeNode* invoke;
        c.invoke(&invoke, fn, FuncSignatureT<float, float>());
        invoke->setArg(0, a);
        invoke->setRet(0, d);
        c.ins(dest.s(element), d.s(0));
    }
}

void CompiledVectorExpression::generateTwoArgCall(a64::Compiler& c, arm::Vec& dest, arm::Vec& arg1, arm::Vec& arg2, float (*function)(float, float)) {
    arm::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    arm::Vec a1 = c.newVecS();
    arm::Vec a2 = c.newVecS();
    arm::Vec d = c.newVecS();
    for (int element = 0; element < width; element++) {
        c.ins(a1.s(0), arg1.s(element));
        c.ins(a2.s(0), arg2.s(element));
        InvokeNode* invoke;
        c.invoke(&invoke, fn, FuncSignatureT<float, float, float>());
        invoke->setArg(0, a1);
        invoke->setArg(1, a2);
        invoke->setRet(0, d);
        c.ins(dest.s(element), d.s(0));
    }
}
#else

union int_to_float {
  int_to_float(const int &_i) { i = _i; }
  int i;
  float  f;
};

void CompiledVectorExpression::generateJitCode() {
    const CpuInfo& cpu = CpuInfo::host();
    if (!cpu.hasFeature(CpuFeatures::X86::kAVX))
        return;
    CodeHolder code;
    code.init(runtime.environment());
    x86::Compiler c(&code);
    FuncNode* funcNode = c.addFunc(FuncSignatureT<void>());
    funcNode->frame().setAvxEnabled();
    vector<x86::Ymm> workspaceVar(workspace.size()/width);
    for (int i = 0; i < (int) workspaceVar.size(); i++)
        workspaceVar[i] = c.newYmmPs();
    x86::Gp argsPointer = c.newIntPtr();
    c.mov(argsPointer, imm(&argValues[0]));
    vector<vector<int> > groups, groupPowers;
    vector<int> stepGroup;
    findPowerGroups(groups, groupPowers, stepGroup);

    // Load the arguments into variables.

    for (set<string>::const_iterator iter = variableNames.begin(); iter != variableNames.end(); ++iter) {
        map<string, int>::iterator index = variableIndices.find(*iter);
        x86::Gp variablePointer = c.newIntPtr();
        c.mov(variablePointer, imm(getVariablePointer(index->first)));
        if (width == 4)
            c.vmovdqu(workspaceVar[index->second].xmm(), x86::ptr(variablePointer, 0, 0));
        else
            c.vmovdqu(workspaceVar[index->second], x86::ptr(variablePointer, 0, 0));
    }

    // Make a list of all constants that will be needed for evaluation.

    vector<int> operationConstantIndex(operation.size(), -1);
    for (int step = 0; step < (int) operation.size(); step++) {
        // Find the constant value (if any) used by this operation.

        Operation& op = *operation[step];
        double value;
        if (op.getId() == Operation::CONSTANT)
            value = dynamic_cast<Operation::Constant&> (op).getValue();
        else if (op.getId() == Operation::ADD_CONSTANT)
            value = dynamic_cast<Operation::AddConstant&> (op).getValue();
        else if (op.getId() == Operation::MULTIPLY_CONSTANT)
            value = dynamic_cast<Operation::MultiplyConstant&> (op).getValue();
        else if (op.getId() == Operation::RECIPROCAL)
            value = 1.0;
        else if (op.getId() == Operation::STEP)
            value = 1.0;
        else if (op.getId() == Operation::DELTA)
            value = 1.0;
        else if (op.getId() == Operation::ABS)
            value = int_to_float(0x7FFFFFFF).f;
        else if (op.getId() == Operation::POWER_CONSTANT) {
            if (stepGroup[step] == -1)
                value = dynamic_cast<Operation::PowerConstant&> (op).getValue();
            else
                value = 1.0;
        } else
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

    vector<x86::Ymm> constantVar(constants.size());
    if (constants.size() > 0) {
        x86::Gp constantsPointer = c.newIntPtr();
        c.mov(constantsPointer, imm(&constants[0]));
        for (int i = 0; i < (int) constants.size(); i++) {
            constantVar[i] = c.newYmmPs();
            c.vbroadcastss(constantVar[i], x86::ptr(constantsPointer, 4*i, 0));
        }
    }

    // Evaluate the operations.

    vector<bool> hasComputedPower(operation.size(), false);
    x86::Ymm argReg = c.newYmm();
    x86::Ymm doubleArgReg = c.newYmm();
    x86::Ymm doubleResultReg = c.newYmm();
    for (int step = 0; step < (int) operation.size(); step++) {
        if (hasComputedPower[step])
            continue;

        // When one or more steps involve raising the same argument to multiple integer
        // powers, we can compute them all together for efficiency.

        if (stepGroup[step] != -1) {
            vector<int>& group = groups[stepGroup[step]];
            vector<int>& powers = groupPowers[stepGroup[step]];
            x86::Ymm multiplier = c.newYmmPs();
            if (powers[0] > 0)
                c.vmovdqu(multiplier, workspaceVar[arguments[step][0]]);
            else {
                c.vdivps(multiplier, constantVar[operationConstantIndex[step]], workspaceVar[arguments[step][0]]);
                for (int i = 0; i < (int)powers.size(); i++)
                    powers[i] = -powers[i];
            }
            vector<bool> hasAssigned(group.size(), false);
            bool done = false;
            while (!done) {
                done = true;
                for (int i = 0; i < (int)group.size(); i++) {
                    if (powers[i] % 2 == 1) {
                        if (!hasAssigned[i])
                            c.vmovdqu(workspaceVar[target[group[i]]], multiplier);
                        else
                            c.vmulps(workspaceVar[target[group[i]]], workspaceVar[target[group[i]]], multiplier);
                        hasAssigned[i] = true;
                    }
                    powers[i] >>= 1;
                    if (powers[i] != 0)
                        done = false;
                }
                if (!done)
                    c.vmulps(multiplier, multiplier, multiplier);
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
                args.push_back(args[0] + i);
        }

        // Generate instructions to execute this operation.

        switch (op.getId()) {
            case Operation::CONSTANT:
                c.vmovdqu(workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::ADD:
                c.vaddps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::SUBTRACT:
                c.vsubps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::MULTIPLY:
                c.vmulps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::DIVIDE:
                c.vdivps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::POWER:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], powf);
                break;
            case Operation::NEGATE:
                c.vxorps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[target[step]]);
                c.vsubps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::SQRT:
                c.vsqrtps(workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::EXP:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], expf);
                break;
            case Operation::LOG:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], logf);
                break;
            case Operation::SIN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], sinf);
                break;
            case Operation::COS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], cosf);
                break;
            case Operation::TAN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], tanf);
                break;
            case Operation::ASIN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], asinf);
                break;
            case Operation::ACOS:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], acosf);
                break;
            case Operation::ATAN:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], atanf);
                break;
            case Operation::ATAN2:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]], atan2f);
                break;
            case Operation::SINH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], sinhf);
                break;
            case Operation::COSH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], coshf);
                break;
            case Operation::TANH:
                generateSingleArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], tanhf);
                break;
            case Operation::STEP:
                c.vxorps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[target[step]]);
                c.vcmpps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]], imm(18)); // Comparison mode is _CMP_LE_OQ = 18
                c.vandps(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::DELTA:
                c.vxorps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[target[step]]);
                c.vcmpps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]], imm(16)); // Comparison mode is _CMP_EQ_OQ = 0
                c.vandps(workspaceVar[target[step]], workspaceVar[target[step]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::SQUARE:
                c.vmulps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
                break;
            case Operation::CUBE:
                c.vmulps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[0]]);
                c.vmulps(workspaceVar[target[step]], workspaceVar[target[step]], workspaceVar[args[0]]);
                break;
            case Operation::RECIPROCAL:
                c.vdivps(workspaceVar[target[step]], constantVar[operationConstantIndex[step]], workspaceVar[args[0]]);
                break;
            case Operation::ADD_CONSTANT:
                c.vaddps(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::MULTIPLY_CONSTANT:
                c.vmulps(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::POWER_CONSTANT:
                generateTwoArgCall(c, workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]], powf);
                break;
            case Operation::MIN:
                c.vminps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::MAX:
                c.vmaxps(workspaceVar[target[step]], workspaceVar[args[0]], workspaceVar[args[1]]);
                break;
            case Operation::ABS:
                c.vandps(workspaceVar[target[step]], workspaceVar[args[0]], constantVar[operationConstantIndex[step]]);
                break;
            case Operation::FLOOR:
                c.vroundps(workspaceVar[target[step]], workspaceVar[args[0]], imm(1));
                break;
            case Operation::CEIL:
                c.vroundps(workspaceVar[target[step]], workspaceVar[args[0]], imm(2));
                break;
            case Operation::SELECT:
            {
                x86::Ymm mask = c.newYmmPs();
                c.vxorps(mask, mask, mask);
                c.vcmpps(mask, mask, workspaceVar[args[0]], imm(0)); // Comparison mode is _CMP_EQ_OQ = 0
                c.vblendvps(workspaceVar[target[step]], workspaceVar[args[1]], workspaceVar[args[2]], mask);
                break;
            }
            default:
                // Just invoke evaluateOperation().

                for (int element = 0; element < width; element++) {
                    for (int i = 0; i < (int) args.size(); i++) {
                        if (element < 4)
                            c.vshufps(argReg, workspaceVar[args[i]], workspaceVar[args[i]], imm(element));
                        else {
                            c.vperm2f128(argReg, workspaceVar[args[i]], workspaceVar[args[i]], imm(1));
                            c.vshufps(argReg, argReg, argReg, imm(element-4));
                        }
                        c.vcvtss2sd(doubleArgReg.xmm(), doubleArgReg.xmm(), argReg.xmm());
                        c.vmovsd(x86::ptr(argsPointer, 8*i, 0), doubleArgReg.xmm());
                    }
                    x86::Gp fn = c.newIntPtr();
                    c.mov(fn, imm((void*) evaluateOperation));
                    InvokeNode* invoke;
                    c.invoke(&invoke, fn, FuncSignatureT<double, Operation*, double*>());
                    invoke->setArg(0, imm(&op));
                    invoke->setArg(1, imm(&argValues[0]));
                    invoke->setRet(0, doubleResultReg);
                    c.vcvtsd2ss(argReg.xmm(), argReg.xmm(), doubleResultReg.xmm());
                    if (element > 3)
                        c.vperm2f128(argReg, argReg, argReg, imm(0));
                    if (element != 0)
                        c.vshufps(argReg, argReg, argReg, imm(0));
                    c.vblendps(workspaceVar[target[step]], workspaceVar[target[step]], argReg, 1<<element);
                }
        }
    }
    x86::Gp resultPointer = c.newIntPtr();
    c.mov(resultPointer, imm(&workspace[workspace.size()-width]));
    if (width == 4)
        c.vmovdqu(x86::ptr(resultPointer, 0, 0), workspaceVar.back().xmm());
    else
        c.vmovdqu(x86::ptr(resultPointer, 0, 0), workspaceVar.back());
    c.endFunc();
    c.finalize();
    runtime.add(&jitCode, &code);
}

void CompiledVectorExpression::generateSingleArgCall(x86::Compiler& c, x86::Ymm& dest, x86::Ymm& arg, float (*function)(float)) {
    x86::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    x86::Ymm a = c.newYmm();
    x86::Ymm d = c.newYmm();
    for (int element = 0; element < width; element++) {
        if (element < 4)
            c.vshufps(a, arg, arg, imm(element));
        else {
            c.vperm2f128(a, arg, arg, imm(1));
            c.vshufps(a, a, a, imm(element-4));
        }
        InvokeNode* invoke;
        c.invoke(&invoke, fn, FuncSignatureT<float, float>());
        invoke->setArg(0, a);
        invoke->setRet(0, d);
        if (element > 3)
            c.vperm2f128(d, d, d, imm(0));
        if (element != 0)
            c.vshufps(d, d, d, imm(0));
        c.vblendps(dest, dest, d, 1<<element);
    }
}

void CompiledVectorExpression::generateTwoArgCall(x86::Compiler& c, x86::Ymm& dest, x86::Ymm& arg1, x86::Ymm& arg2, float (*function)(float, float)) {
    x86::Gp fn = c.newIntPtr();
    c.mov(fn, imm((void*) function));
    x86::Ymm a1 = c.newYmm();
    x86::Ymm a2 = c.newYmm();
    x86::Ymm d = c.newYmm();
    for (int element = 0; element < width; element++) {
        if (element < 4) {
            c.vshufps(a1, arg1, arg1, imm(element));
            c.vshufps(a2, arg2, arg2, imm(element));
        }
        else {
            c.vperm2f128(a1, arg1, arg1, imm(1));
            c.vperm2f128(a2, arg2, arg2, imm(1));
            c.vshufps(a1, a1, a1, imm(element-4));
            c.vshufps(a2, a2, a2, imm(element-4));
        }
        InvokeNode* invoke;
        c.invoke(&invoke, fn, FuncSignatureT<float, float, float>());
        invoke->setArg(0, a1);
        invoke->setArg(1, a2);
        invoke->setRet(0, d);
        if (element > 3)
            c.vperm2f128(d, d, d, imm(0));
        if (element != 0)
            c.vshufps(d, d, d, imm(0));
        c.vblendps(dest, dest, d, 1<<element);
    }
}
#endif
#endif
