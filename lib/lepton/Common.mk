# -*- makefile -*-
# common settings for Lepton library makefiles

SRC= \
    src/CompiledExpression.cpp \
    src/CompiledVectorExpression.cpp \
    src/ExpressionProgram.cpp \
    src/ExpressionTreeNode.cpp \
    src/Operation.cpp \
    src/ParsedExpression.cpp \
    src/Parser.cpp
OBJ=$(SRC:src/%.cpp=build/lepton.%.o)

JITARM= \
    asmjit/arm/a64assembler.cpp \
    asmjit/arm/a64builder.cpp \
    asmjit/arm/a64compiler.cpp \
    asmjit/arm/a64emithelper.cpp \
    asmjit/arm/a64formatter.cpp \
    asmjit/arm/a64func.cpp \
    asmjit/arm/a64instapi.cpp \
    asmjit/arm/a64instdb.cpp \
    asmjit/arm/a64operand.cpp \
    asmjit/arm/a64rapass.cpp \
    asmjit/arm/armformatter.cpp
JITX86 = \
    asmjit/x86/x86assembler.cpp \
    asmjit/x86/x86builder.cpp \
    asmjit/x86/x86compiler.cpp \
    asmjit/x86/x86emithelper.cpp \
    asmjit/x86/x86formatter.cpp \
    asmjit/x86/x86func.cpp \
    asmjit/x86/x86instapi.cpp \
    asmjit/x86/x86instdb.cpp \
    asmjit/x86/x86operand.cpp \
    asmjit/x86/x86rapass.cpp
JITCORE= \
    asmjit/core/archtraits.cpp \
    asmjit/core/assembler.cpp \
    asmjit/core/builder.cpp \
    asmjit/core/codeholder.cpp \
    asmjit/core/codewriter.cpp \
    asmjit/core/compiler.cpp \
    asmjit/core/constpool.cpp \
    asmjit/core/cpuinfo.cpp \
    asmjit/core/emithelper.cpp \
    asmjit/core/emitter.cpp \
    asmjit/core/emitterutils.cpp \
    asmjit/core/environment.cpp \
    asmjit/core/errorhandler.cpp \
    asmjit/core/formatter.cpp \
    asmjit/core/funcargscontext.cpp \
    asmjit/core/func.cpp \
    asmjit/core/globals.cpp \
    asmjit/core/inst.cpp \
    asmjit/core/jitallocator.cpp \
    asmjit/core/jitruntime.cpp \
    asmjit/core/logger.cpp \
    asmjit/core/operand.cpp \
    asmjit/core/osutils.cpp \
    asmjit/core/ralocal.cpp \
    asmjit/core/rapass.cpp \
    asmjit/core/rastack.cpp \
    asmjit/core/string.cpp \
    asmjit/core/support.cpp \
    asmjit/core/target.cpp \
    asmjit/core/type.cpp \
    asmjit/core/virtmem.cpp \
    asmjit/core/zone.cpp \
    asmjit/core/zonehash.cpp \
    asmjit/core/zonelist.cpp \
    asmjit/core/zonestack.cpp \
    asmjit/core/zonetree.cpp \
    asmjit/core/zonevector.cpp

JITOBJ=$(JITX86:asmjit/x86/%.cpp=build/x86.%.o) \
       $(JITARM:asmjit/arm/%.cpp=build/arm.%.o) \
       $(JIXCORE:asmjit/core/%.cpp=build/core.%.o)

LEPTON_DIR=.

include $(LEPTON_DIR)/Settings.mk

EXTRAMAKE=Makefile.lammps.empty
LIB=liblepton.a

ifeq ($(ENABLE_JIT),1)
OBJ += $(JITOBJ)
endif

INC += $(LEPTON_INC)
CXXFLAGS += $(LEPTON_DEF)

all: $(LIB) Makefile.lammps

build:
	mkdir -p build

build/lepton.%.o: src/%.cpp build
	$(CXX) $(INC) $(CXXFLAGS) -c $< -o $@

build/arm.%.o: asmjit/arm/%.cpp build
	$(CXX) $(INC) $(CXXFLAGS) -c $< -o $@

build/x86.%.o: asmjit/x86/%.cpp build
	$(CXX) $(INC) $(CXXFLAGS) -c $< -o $@

build/core.%.o: asmjit/core/%.cpp build
	$(CXX) $(INC) $(CXXFLAGS) -c $< -o $@

Makefile.lammps:
	cp $(EXTRAMAKE) $@
	sed -i -e 's,^.*lepton_SYSINC *=.*$$,lepton_SYSINC = $(DEF),' $@

.PHONY: all lib clean

$(LIB) : $(OBJ)
	$(AR) $(ARFLAGS) $@ $^

clean:
	rm -f build/*.o $(LIB) *~ Makefile.lammps


