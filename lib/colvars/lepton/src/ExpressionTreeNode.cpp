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

#include "lepton/ExpressionTreeNode.h"
#include "lepton/Exception.h"
#include "lepton/Operation.h"

using namespace Lepton;
using namespace std;

ExpressionTreeNode::ExpressionTreeNode(Operation* operation, const vector<ExpressionTreeNode>& children) : operation(operation), children(children) {
    if (operation->getNumArguments() != children.size())
        throw Exception("wrong number of arguments to function: "+operation->getName());
}

ExpressionTreeNode::ExpressionTreeNode(Operation* operation, const ExpressionTreeNode& child1, const ExpressionTreeNode& child2) : operation(operation) {
    children.push_back(child1);
    children.push_back(child2);
    if (operation->getNumArguments() != children.size())
        throw Exception("wrong number of arguments to function: "+operation->getName());
}

ExpressionTreeNode::ExpressionTreeNode(Operation* operation, const ExpressionTreeNode& child) : operation(operation) {
    children.push_back(child);
    if (operation->getNumArguments() != children.size())
        throw Exception("wrong number of arguments to function: "+operation->getName());
}

ExpressionTreeNode::ExpressionTreeNode(Operation* operation) : operation(operation) {
    if (operation->getNumArguments() != children.size())
        throw Exception("wrong number of arguments to function: "+operation->getName());
}

ExpressionTreeNode::ExpressionTreeNode(const ExpressionTreeNode& node) : operation(node.operation == NULL ? NULL : node.operation->clone()), children(node.getChildren()) {
}

ExpressionTreeNode::ExpressionTreeNode() : operation(NULL) {
}

ExpressionTreeNode::~ExpressionTreeNode() {
    if (operation != NULL)
        delete operation;
}

bool ExpressionTreeNode::operator!=(const ExpressionTreeNode& node) const {
    if (node.getOperation() != getOperation())
        return true;
    if (getOperation().isSymmetric() && getChildren().size() == 2) {
        if (getChildren()[0] == node.getChildren()[0] && getChildren()[1] == node.getChildren()[1])
            return false;
        if (getChildren()[0] == node.getChildren()[1] && getChildren()[1] == node.getChildren()[0])
            return false;
        return true;
    }
    for (int i = 0; i < (int) getChildren().size(); i++)
        if (getChildren()[i] != node.getChildren()[i])
            return true;
    return false;
}

bool ExpressionTreeNode::operator==(const ExpressionTreeNode& node) const {
    return !(*this != node);
}

ExpressionTreeNode& ExpressionTreeNode::operator=(const ExpressionTreeNode& node) {
    if (operation != NULL)
        delete operation;
    operation = node.getOperation().clone();
    children = node.getChildren();
    return *this;
}

const Operation& ExpressionTreeNode::getOperation() const {
    return *operation;
}

const vector<ExpressionTreeNode>& ExpressionTreeNode::getChildren() const {
    return children;
}
