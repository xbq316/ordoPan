#include <iostream>
#include <vector>
#include <unordered_map>
#include <limits>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <unordered_set>
#include <cmath>
#include <stdexcept>

class Graph {
public:
    Graph() {} // Default construction, nodes can be added dynamically.

    // Add a node and return the node number.
    int addNode() {
        adj.emplace_back();
        revAdj.emplace_back();
        maxNext.emplace_back(-1, -std::numeric_limits<float>::infinity());
        maxPrev.emplace_back(-1, -std::numeric_limits<float>::infinity());
        return adj.size() - 1;
    }

    std::vector<int> getNextNodesAboveThreshold(size_t u,const float& threshold) const {
        std::vector<int> result;
        if (u >= adj.size()) return result;
        for (const auto& [v, w] : adj[u]) {
            if (w > threshold) result.push_back(v);
        }
        return result;
    }

    // Return all previous nodes of v with weights greater than the threshold.
    std::vector<int> getPrevNodesAboveThreshold(size_t v,const float& threshold) const {
        std::vector<int> result;
        if (v >= revAdj.size()) return result;
        for (const auto& [u, w] : revAdj[v]) {
            if (w >= threshold) result.push_back(u);
        }
        return result;
    }

    // Add a directed edge u -> v with weight w
    void addEdge(size_t u, size_t v, float w) {
        ensureNode(u);
        ensureNode(v);
        adj[u][v] = w;
        revAdj[v][u] = w;

        // Update maximum outgoing edge
        if (w > maxNext[u].second) {
            maxNext[u] = {v, w};
        }

        // Update maximum inbound edge
        if (w > maxPrev[v].second) {
            maxPrev[v] = {u, w};
        }
    }

    // O(1) returns the weights of u -> v
    float getWeight(size_t u, size_t v) const {
        if (u >= adj.size()) return std::numeric_limits<float>::quiet_NaN();
        auto it = adj[u].find(v);
        if (it != adj[u].end()) return it->second;
        return std::numeric_limits<float>::quiet_NaN();
    }

    // O(1) returns the largest outgoing edge node of u.
    int getMaxWeightNextNode(size_t u) const {
        if (u >= maxNext.size()) return -1;
        return maxNext[u].first;
    }

    // O(1) returns the largest incoming edge node of v.
    int getMaxWeightPrevNode(size_t v) const {
        if (v >= maxPrev.size()) return -1;
        return maxPrev[v].first;
    }

private:
    std::vector<std::unordered_map<int, float>> adj;    // Out of bounds
    std::vector<std::unordered_map<int, float>> revAdj; // Entering edge
    std::vector<std::pair<int, float>> maxNext; // Maximum outgoing edges for each node
    std::vector<std::pair<int, float>> maxPrev; // Maximum inbound edge of each node

    // Ensure the node exists
    void ensureNode(size_t u) {
        while (u >= adj.size()) addNode();
    }
};

class Cluster {
public:
    std::vector<int> start;
    std::vector<int> end;
    size_t speciesNum;
    std::vector<std::vector<int>> content;


    // Constructor: Single-line initialization
    Cluster(const std::vector<int>& line) {
        start = line;
        end = line;
        speciesNum = line.size();
        content.push_back(line);
    }

    // Multi-line initialization
    static Cluster fromMultiline(const std::vector<std::vector<int>>& lines) {
        if (lines.empty())
            throw std::invalid_argument("Input lines cannot be empty");

        Cluster newCluster({});
        newCluster.content = lines;
        newCluster.start = newCluster.updateStart();
        newCluster.end = newCluster.updateEnd();
        newCluster.speciesNum = lines[0].size();
        return newCluster;
    }

    // renew start
    std::vector<int> updateStart() const {
        std::vector<int> newLine = content[0];
        for (const auto& line : content) {
            for (size_t i = 0; i < newLine.size(); ++i) {
                if (newLine[i] == -9) newLine[i] = line[i];
            }
        }
        return newLine;
    }

    // renew end
    std::vector<int> updateEnd() const {
        std::vector<int> newLine = content.back();
        for (auto it = content.rbegin(); it != content.rend(); ++it) {
            for (size_t i = 0; i < newLine.size(); ++i) {
                if (newLine[i] == -9) newLine[i] = (*it)[i];
            }
        }
        return newLine;
    }

    // Merge two Cluster
    Cluster operator+(const Cluster& nextCluster) const {
        std::vector<std::vector<int>> combined = content;
        combined.insert(combined.end(), nextCluster.content.begin(), nextCluster.content.end());
        return Cluster::fromMultiline(combined);
    }

    // Determine if they are adjacent
    bool isLinked(const Cluster& nextCluster) const {
        for (size_t i = 0; i < end.size(); ++i) {
            int lastEnd = end[i];
            int thisStart = nextCluster.start[i];
            if (lastEnd == -9 || thisStart == -9) continue;
            if (lastEnd + 1 == thisStart) return true;
        }
        return false;
    }

    // Number of supported species
    int supportedSpecies() const {
        int count = 0;
        for (int s : start) {
            if (s != -9) ++count;
        }
        return count;
    }

    float link_rate(const Cluster& nextCluster) const {
        //Calculate the adjacency rate:
        //Number of adjacent species / (number of species that can be compared))
        //Number of comparable species: If none are missing, increment by 1.
        size_t norm_num = 0;
        size_t adjacency = 0;

        for (size_t i = 0; i < this->speciesNum; i++){
            int nextStart = nextCluster.start[i];
            int thisEnd = this->end[i];
            if (thisEnd == -9 || nextStart == -9) continue;
            if (std::abs(thisEnd - nextStart) == 1) adjacency++;
            norm_num++;
        }
        if (norm_num==0) return 0;
        else return (float)adjacency/norm_num;
    }

    float linkScore(const Cluster& ortherCluster) const {
        float result = 0;
        for (size_t i = 0; i < this->speciesNum; i++){
            if (this->end[i] == -9 || ortherCluster.start[i] == -9) { continue; }
            int diff = ortherCluster.start[i] - this->end[i];
            if (diff < 0) { continue; }
            else{result += std::exp(1-diff);}
        }
        return result;
    }

};

bool isNormMatrix(const std::vector<std::vector<std::vector<int>>> &clusterVector) {
    size_t speciesNum;
    if (clusterVector.empty()){return false ;}
    if (clusterVector[0].empty()){return false ;}
    if (clusterVector[0][0].empty()){return false ;}
    speciesNum = clusterVector[0][0].size();
    
    for (size_t i = 0; i < clusterVector.size(); i++) {
        for (size_t j = 0; j < clusterVector[i].size(); j++) {
            if (clusterVector[i][j].size() != speciesNum) {
                return false;
            }
        }
    }
    return true;
}

std::vector<std::vector<std::vector<int>>> readCluster(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::vector<std::vector<std::vector<int>>> clusters;
    std::vector<std::vector<int>> current_cluster;
    std::string line;

    while (std::getline(file, line)) {
        // Skip blank lines
        if (line.empty()) continue;

        // Check if it is a new Cluster block (e.g., "#Cluster-1843:").
        if (line[0]=='#') {
            if (!current_cluster.empty()) {
                clusters.push_back(current_cluster);
                current_cluster.clear();
            }
        }
        else {
            // Parse rows of numbers (separated by tabs).
            std::vector<int> row;
            std::istringstream iss(line);
            std::string token;

            while (std::getline(iss, token, '\t')) {
                if (token == "-") {row.push_back(-9);}
                else {row.push_back(std::stoi(token));}
            }
            current_cluster.push_back(row);
        }
    }

    // Add the last cluster
    if (!current_cluster.empty()) {
        clusters.push_back(current_cluster);
    }
    return clusters;
}

std::string vec2Key(const std::vector<int>& v) {
    std::string s;
    s.resize(v.size() * sizeof(int)); // Fixed size
    std::memcpy(s.data(), v.data(), v.size() * sizeof(int));
    return s;
}


bool hasDuplicate(const std::vector<size_t>& v) {
    std::unordered_set<size_t> seen;
    for (auto x : v) {
        if (!seen.insert(x).second) {
            return true; // Insertion failure indicates a duplicate.
        }
    }
    return false;
}


std::vector<size_t> geneindex2clusterindex(const std::vector<std::vector<int>>& geneindex,
    const std::vector<Cluster> &clusterVector) {
    
    std::unordered_map<std::string, size_t> toClusterindex;
    for(size_t i = 0; i < clusterVector.size(); i++){
        for(const std::vector<int>& line : clusterVector[i].content){
            toClusterindex[vec2Key(line)] = i;
        }
    }

    std::vector<size_t> clusterindex;
    for (size_t i = 0; i < geneindex.size(); i++) {
        std::string key = vec2Key(geneindex[i]);
        auto it = toClusterindex.find(key);
        if (it != toClusterindex.end()) {
            // If the current value is not equal to the last value of clusterindex, then add it.
            if (clusterindex.empty() || clusterindex.back() != it->second) {
                clusterindex.push_back(it->second);
            }
            
        }
        else {
            //Error exit
            std::cout << "Error: gene " << i << " not found in any cluster" << std::endl;
            return std::vector<size_t>();
        }
    }
    if (hasDuplicate(clusterindex)) {
        std::cout << "cluster index error: has duplicate." << std::endl;
        return std::vector<size_t>();
    }
    
    return clusterindex;
}

void writeClusterVector(const std::vector<Cluster>& clusterVector,
    const std::string& outFile) {
    std::ofstream out(outFile);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file " << outFile << " for writing." << std::endl;
        return;
    }

    // Traversing a 3D vector
    for (size_t i = 0; i < clusterVector.size(); i++) {
        for (size_t j = 0; j < clusterVector[i].content.size(); j++) {
            // Traverse the innermost vector (numerical value).
            for (size_t k = 0; k < clusterVector[i].content[j].size(); k++) {
                const int& value = clusterVector[i].content[j][k];
                out << (value == -9 ? "-" : std::to_string(value));
                if (k != clusterVector[i].content[j].size() - 1) {
                    out << "\t";
                }
            }
            out << "\n";  // End the current subcluster
        }
    }
    out.close();
}

std::vector<std::vector<int>> clusters2geneindex(const std::vector<Cluster> &clusterVector) {
    std::vector<std::vector<int>> geneindex;
    for (const Cluster& cluster : clusterVector) {
        for (const std::vector<int>& line : cluster.content) {
            geneindex.push_back(line);
        }
    }
    return geneindex;
}


std::vector<std::vector<int>> readGeneIndex(const std::string& file) {
    std::vector<std::vector<int>> result;
    std::ifstream fin(file);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open file: " + file);
    }

    std::string line;
    while (std::getline(fin, line)) {
        // Remove the first and last blank spaces
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;
        std::vector<int> row;

        while (iss >> token) {
            if (token == "-") {
                row.push_back(-9);
            }
            else {
                row.push_back(std::stoi(token));
            }
        }

        if (!row.empty()) {
            result.push_back(row);
        }
    }

    return result;
}

// Example
int main(int argc, char* argv[]) {
    std::string clusterFile;
    std::string geneindexFile;
	std::string outFile;
    size_t geneThreshold=10;

    bool hasCluster = false,hasGeneindex = false, hasOutput = false;
    for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-c" && i + 1 < argc) {
			clusterFile = argv[++i];
			hasCluster = true;
		}
        else if (arg == "-g" && i + 1 < argc) {
			geneindexFile = argv[++i];
            hasGeneindex = true;
		}
		else if (arg == "-o" && i + 1 < argc) {
			outFile = argv[++i];
			hasOutput = true;
		}
        else if (arg == "-n" && i + 1 < argc) {
			geneThreshold = atoi(argv[++i]);
		}
        else if (arg == "-h"){
            std::cout << "Usage: " << argv[0] << " -c input.cluster -g input.geneindex -n 10-o output.geneindex" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "-c input.cluster: input cluster file(cluster program output)" << std::endl;
            std::cout << "-g input.geneindex: input geneindex file" << std::endl;
            std::cout << "-n 10: gene threshold(Minimum cluster block size)" << std::endl;
            std::cout << "-o output.geneindex: output geneindex file" << std::endl;
            return 0;
        }
        
		else {
            std::cout << "error: unknown argument: " << arg << std::endl;
            std::cout << "Usage: " << argv[0] << " -c input.cluster -g input.geneindex -n 10-o output.geneindex" << std::endl;
			return 1;
		}
	}

    // Check required parameters
	if (!hasCluster || !hasGeneindex || !hasOutput) {
        std::cout << "error: must provide (-i and -g and -o)" << std::endl;
        std::cout << "Usage: " << argv[0] << " -c input.cluster -g input.geneindex -n 10 -o output.geneindex" << std::endl;
		return 1;
	}

    // Program parameters...
    std::cout << "cluster file: " << clusterFile << std::endl;
    std::cout << "geneindex file: " << geneindexFile << std::endl;
    std::cout << "output file: " << outFile << std::endl;
    

    // Read cluster data
    std::vector<std::vector<std::vector<int>>> clusterMatrix= readCluster(clusterFile);
    if (!isNormMatrix(clusterMatrix)){
        std::cout<<"error: input file is not norm matrix"<<std::endl;
        return 1;
    }
    size_t clusterNum=clusterMatrix.size();
    std::vector<Cluster> clusterVector;
    for (size_t i=0;i<clusterNum;i++){
		clusterVector.emplace_back(Cluster::fromMultiline(clusterMatrix[i]));
	}

    // Read the geneindex file
    std::vector<std::vector<int>> geneindex = readGeneIndex(geneindexFile);
    if (geneindex.empty()) {
        std::cout << "error: geneindex file is empty" << std::endl;
        return 1;
    }

    //Get clusterindex
    std::vector<size_t> clusterIndex = geneindex2clusterindex(geneindex, clusterVector);
    if (clusterIndex.empty()) {
        std::cout << "error: cluster index is empty" << std::endl;
        return 1;
    }
    
    // Clustering
    std::vector<Cluster> newClusterVector;
    Cluster lastCluster = clusterVector[clusterIndex[0]];
    for (size_t i = 1; i < clusterIndex.size(); i++){
        size_t index = clusterIndex[i];
        const Cluster & thisCluster = clusterVector[index];
        float score = lastCluster.link_rate(thisCluster);
        if (score > 0){lastCluster = lastCluster + thisCluster;}
        else{
            newClusterVector.push_back(lastCluster);
            lastCluster = thisCluster;
        }
    }
    newClusterVector.push_back(lastCluster);

    //Separate cluster
    std::vector<Cluster> retainClusters;
    std::vector<Cluster> splitClusters;
    for (size_t i = 0; i < newClusterVector.size(); i++){
        if (newClusterVector[i].content.size()<=geneThreshold){splitClusters.push_back(newClusterVector[i]);}
        else{retainClusters.push_back(newClusterVector[i]);}
    }
    std::vector<std::vector<int>> retainGeneindex=clusters2geneindex(retainClusters);
    std::vector<size_t> retainIndex=geneindex2clusterindex(retainGeneindex,clusterVector);
    //Check data
    for (size_t i = 0; i < retainIndex.size(); i++){
        if(retainIndex[i]>=clusterNum){
            std::cout<<"error: retainIndex is out of range"<<std::endl;
            return 1;
        }
    }
    std::vector<std::vector<int>> splitGeneindex=clusters2geneindex(splitClusters);
    std::vector<size_t> splitIndex=geneindex2clusterindex(splitGeneindex,clusterVector);
    //Check data
    for (size_t i = 0; i < splitIndex.size(); i++){
        if(splitIndex[i]>=clusterNum){
            std::cout<<"error: splitIndex is out of range"<<std::endl;
            return 1;
        }
    }

    //Start Insertion Sort
    for (size_t i = 0; i < splitIndex.size(); i++){
        size_t insertIndex=splitIndex[i];
        const Cluster & insertCluster=clusterVector[insertIndex];
        Cluster lastCluster=clusterVector[retainIndex[0]];
        float bestScore=insertCluster.linkScore(lastCluster);
        size_t bestIndex=0;
        bool flag=false;

        for (size_t j = 0; j < retainIndex.size(); j++){
            const Cluster & thisCluster=clusterVector[retainIndex[j]];
            float score=thisCluster.linkScore(insertCluster);
            if (score>bestScore){bestScore=score;bestIndex=j;flag=true;}
        }
        if (flag){bestIndex++;}
        retainIndex.insert(retainIndex.begin()+bestIndex,insertIndex);
    }
    

    //Output
    std::vector<Cluster> result ;
    for (size_t i = 0; i < retainIndex.size(); i++){
        result.push_back(clusterVector[retainIndex[i]]);
    }
    writeClusterVector(result, outFile);

    return 0;
}
