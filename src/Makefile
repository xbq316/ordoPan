CXX = g++
CXXFLAGS = -O3 -fomit-frame-pointer -Wno-deprecated -std=c++17

PROGS = dagchainer clusterSort cluster finalAdjustment

all: $(PROGS)

dagchainer: dagchainer.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

cluster: cluster.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clusterSort: clusterSort.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

finalAdjustment: finalAdjustment.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
clean:
	rm -f $(PROGS)

.PHONY: all clean