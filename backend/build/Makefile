# ===============================
# Compiler
# ===============================
CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -pedantic-errors -g \
            -msse2 -march=native -funroll-loops \
            -fomit-frame-pointer -fstrict-aliasing
ifeq ($(OS),Windows_NT)
EXE_EXT := .exe
LDFLAGS := -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive,-Bdynamic

else
LDFLAGS :=
endif

# ===============================
# Source directory (unchanged)
# ===============================
SDIR := ../src
SRCS := $(wildcard $(SDIR)/*.cpp)

# ===============================
# Object files in PWD
# ===============================
OBJS := $(notdir $(SRCS:.cpp=.o))

TARGET := sse_sim_test$(EXE_EXT)

# ===============================
# Default target
# ===============================
all: $(TARGET)
	@echo Compilation Successful

# ===============================
# Link
# ===============================
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# ===============================
# Compile EACH source separately
# ===============================
%.o: $(SDIR)/%.cpp $(wildcard ../include/*.hpp)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ===============================
# Clean (PWD only)
# ===============================
clean:
ifeq ($(OS),Windows_NT)
	@cmd /c "del /Q *.o 2>NUL"
	@cmd /c "del /Q $(TARGET) 2>NUL"
else
	@rm -f *.o $(TARGET)
endif

# ===============================
# Debug target — enables SSE_DEBUG
# Recompiles everything from scratch with -DSSE_DEBUG so the diagnostic
# output (init check, op-string dump, linked list, loop BFS trace) is active.
# Usage:  make debug     (from backend/build/)
# ===============================
.PHONY: debug
debug:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DSSE_DEBUG"

# ===============================
# Verbose target — enables VERBOSE print output
# Recompiles with -DVERBOSE so parameter dumps and [SSE-CHECK] lines are shown.
# Usage:  make verbose   (from backend/build/)
# ===============================
.PHONY: verbose
verbose:
	$(MAKE) clean
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -DVERBOSE"

# ===============================
# Debug binary — separate binary with SSE_DEBUG enabled
# Does NOT overwrite sse_sim_test.exe.
# Usage:  make debugbin   (from backend/build/)
# ===============================
DBG_TARGET := sse_sim_debug$(EXE_EXT)
.PHONY: debugbin
debugbin: $(DBG_TARGET)

$(DBG_TARGET): $(SRCS) $(wildcard ../include/*.hpp)
	$(CXX) $(CXXFLAGS) -DSSE_DEBUG $(SRCS) $(LDFLAGS) -o $(DBG_TARGET)
	@echo Debug binary built: $(DBG_TARGET)

.PHONY: all clean

.PHONY: benchmark
benchmark: sse_benchmark$(EXE_EXT)

.PHONY: layered
layered: sse_layered$(EXE_EXT)

sse_layered$(EXE_EXT): ../benchmarks/layered.cpp $(filter-out ../src/main.cpp,$(SRCS)) $(wildcard ../include/*.hpp)
	$(CXX) -O2 -std=c++17 -Wall ../benchmarks/layered.cpp ../src/updates.cpp ../src/variables.cpp ../src/observables.cpp $(LDFLAGS) -o $@

.PHONY: layered-contract
layered-contract: layered_contract$(EXE_EXT)

layered_contract$(EXE_EXT): ../benchmarks/layered_contract.cpp ../benchmarks/layered.cpp ../src/variables.cpp $(wildcard ../include/*.hpp)
	$(CXX) -O2 -std=c++17 -Wall ../benchmarks/layered_contract.cpp ../src/variables.cpp $(LDFLAGS) -o $@

.PHONY: projector
projector: sse_projector$(EXE_EXT) emc_validation$(EXE_EXT)

emc_validation$(EXE_EXT): ../benchmarks/emc.cpp ../src/variables.cpp $(wildcard ../include/*.hpp)
	$(CXX) -O2 -std=c++17 -Wall ../benchmarks/emc.cpp ../src/variables.cpp $(LDFLAGS) -o $@

sse_projector$(EXE_EXT): ../benchmarks/projector.cpp ../src/updates.cpp ../src/variables.cpp ../src/observables.cpp $(wildcard ../include/*.hpp)
	$(CXX) -O2 -std=c++17 -Wall ../benchmarks/projector.cpp ../src/updates.cpp ../src/variables.cpp ../src/observables.cpp $(LDFLAGS) -o $@
sse_benchmark$(EXE_EXT): ../benchmarks/runner.cpp $(filter-out ../src/main.cpp,$(SRCS)) $(wildcard ../include/*.hpp)
	$(CXX) -O2 -std=c++17 -Wall ../benchmarks/runner.cpp ../src/updates.cpp ../src/variables.cpp ../src/observables.cpp $(LDFLAGS) -o $@
