BINDIR = bin
SRCDIR = src
CXX = g++
CXXFLAGS = -O3 -pedantic -Wall -std=c++11 -static

TARGETS = hotspot2_part1 hotspot2_part2 resolveOverlapsInSummit-CenteredPeaks findVarWidthPeaks
EXE = $(addprefix $(BINDIR)/,$(TARGETS))

default: $(EXE)

$(BINDIR)/% : $(SRCDIR)/%.cpp
	mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(EXE)
