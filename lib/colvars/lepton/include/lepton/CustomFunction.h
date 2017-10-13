#ifndef LEPTON_CUSTOM_FUNCTION_H_
#define LEPTON_CUSTOM_FUNCTION_H_

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

namespace Lepton {

/**
 * This class is the interface for defining your own function that may be included in expressions.
 * To use it, create a concrete subclass that implements all of the virtual methods for each new function
 * you want to define.  Then when you call Parser::parse() to parse an expression, pass a map of
 * function names to CustomFunction objects.
 */

class LEPTON_EXPORT CustomFunction {
public:
    virtual ~CustomFunction() {
    }
    /**
     * Get the number of arguments this function expects.
     */
    virtual int getNumArguments() const = 0;
    /**
     * Evaluate the function.
     *
     * @param arguments    the array of argument values
     */
    virtual double evaluate(const double* arguments) const = 0;
    /**
     * Evaluate a derivative of the function.
     *
     * @param arguments    the array of argument values
     * @param derivOrder   an array specifying the number of times the function has been differentiated
     *                     with respect to each of its arguments.  For example, the array {0, 2} indicates
     *                     a second derivative with respect to the second argument.
     */
    virtual double evaluateDerivative(const double* arguments, const int* derivOrder) const = 0;
    /**
     * Create a new duplicate of this object on the heap using the "new" operator.
     */
    virtual CustomFunction* clone() const = 0;
};

} // namespace Lepton

#endif /*LEPTON_CUSTOM_FUNCTION_H_*/
