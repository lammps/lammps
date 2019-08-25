#ifndef LEPTON_EXPRESSION_TREE_NODE_H_
#define LEPTON_EXPRESSION_TREE_NODE_H_

/* -------------------------------------------------------------------------- *
 *                                   Lepton                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the Lepton expression parser originating from              *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "windowsIncludes.h"
#include <string>
#include <vector>

namespace Lepton {

class Operation;

/**
 * This class represents a node in the abstract syntax tree representation of an expression.
 * Each node is defined by an Operation and a set of children.  When the expression is
 * evaluated, each child is first evaluated in order, then the resulting values are passed
 * as the arguments to the Operation's evaluate() method.
 */

class LEPTON_EXPORT ExpressionTreeNode {
public:
    /**
     * Create a new ExpressionTreeNode.
     *
     * @param operation    the operation for this node.  The ExpressionTreeNode takes over ownership
     *                     of this object, and deletes it when the node is itself deleted.
     * @param children     the children of this node
     */
    ExpressionTreeNode(Operation* operation, const std::vector<ExpressionTreeNode>& children);
    /**
     * Create a new ExpressionTreeNode with two children.
     *
     * @param operation    the operation for this node.  The ExpressionTreeNode takes over ownership
     *                     of this object, and deletes it when the node is itself deleted.
     * @param child1       the first child of this node
     * @param child2       the second child of this node
     */
    ExpressionTreeNode(Operation* operation, const ExpressionTreeNode& child1, const ExpressionTreeNode& child2);
    /**
     * Create a new ExpressionTreeNode with one child.
     *
     * @param operation    the operation for this node.  The ExpressionTreeNode takes over ownership
     *                     of this object, and deletes it when the node is itself deleted.
     * @param child        the child of this node
     */
    ExpressionTreeNode(Operation* operation, const ExpressionTreeNode& child);
    /**
     * Create a new ExpressionTreeNode with no children.
     *
     * @param operation    the operation for this node.  The ExpressionTreeNode takes over ownership
     *                     of this object, and deletes it when the node is itself deleted.
     */
    ExpressionTreeNode(Operation* operation);
    ExpressionTreeNode(const ExpressionTreeNode& node);
    ExpressionTreeNode();
    ~ExpressionTreeNode();
    bool operator==(const ExpressionTreeNode& node) const;
    bool operator!=(const ExpressionTreeNode& node) const;
    ExpressionTreeNode& operator=(const ExpressionTreeNode& node);
    /**
     * Get the Operation performed by this node.
     */
    const Operation& getOperation() const;
    /**
     * Get this node's child nodes.
     */
    const std::vector<ExpressionTreeNode>& getChildren() const;
private:
    Operation* operation;
    std::vector<ExpressionTreeNode> children;
};

} // namespace Lepton

#endif /*LEPTON_EXPRESSION_TREE_NODE_H_*/
