#include<vector>
#include<iostream>
#include<string>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
using namespace std::chrono;


class Graph {
private:
    // Adjacency list
    std::unordered_map<int, std::vector<int>> outEdges; // Outward edge 
    std::unordered_map<int, std::vector<int>> inEdges;  // Inward edge

public:
    void addNode(int u) {
        outEdges[u];  // Ensure the key exists
        inEdges[u];
    }

    // Add edge u -> v
    void addEdge(int u, int v) {
        outEdges[u].push_back(v);
        inEdges[v].push_back(u);
    }

    // Return the list of outgoing edge nodes
    const std::vector<int>& nextNodes(int u) const {
        static std::vector<int> empty;
        auto it = outEdges.find(u);
        if (it != outEdges.end()) return it->second;
        return empty;
    }

    // Return to the list of inbound edge nodes
    const std::vector<int>& PrevNodes(int u) const {
        static std::vector<int> empty;
        auto it = inEdges.find(u);
        if (it != inEdges.end()) return it->second;
        return empty;
    }

    // Return out degree
    int outDegree(int u) const {
        auto it = outEdges.find(u);
        return it != outEdges.end() ? it->second.size() : 0;
    }

    // Return in-degree
    int inDegree(int u) const {
        auto it = inEdges.find(u);
        return it != inEdges.end() ? it->second.size() : 0;
    }
};


class Cluster {

public:
    std::vector<int> start;
    std::vector<int> end;
    int speciesNum;
    std::vector<std::vector<int>> content;


    // Constructor: Single-line initialization
    Cluster(const std::vector<int>& line) {
        start = line;
        end = line;
        speciesNum = line.size();
        content.push_back(line);
    }

    // Initialize from multiple lines
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

    // Merging two clusters
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
};

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


std::vector<int> imitateStart(const std::vector<std::vector<int>>& geneindex) {
    if (geneindex.empty()) return {};

    int speciesNum = geneindex[0].size();
    std::vector<int> start(speciesNum, -9);

    for (int speciesIndex = 0; speciesIndex < speciesNum; ++speciesIndex) {
        std::vector<int> a;
        for (size_t geneIndex = 0; geneIndex < geneindex.size(); ++geneIndex) {
            int val = geneindex[geneIndex][speciesIndex];
            if (val != -9) {
                a.push_back(val);
            }
        }
        if (!a.empty()) {
            start[speciesIndex] = *std::min_element(a.begin(), a.end()) - 1;
        }
    }
    return start;
}

bool removeImitateStartFromClusters(std::vector<Cluster>& clusters, const std::vector<int>& start) {
    bool startIsIsolated = false;
    int isolatedIndex = -1;
    bool removeImitateStart = false;

    for (size_t index = 0; index < clusters.size(); ++index) {
        Cluster& cluster = clusters[index];
        if (cluster.start == start) {
            if (cluster.content.size() > 1) {
                cluster.content.erase(cluster.content.begin());  // Remove the first line
                cluster.start = cluster.updateStart();
                cluster.end = cluster.updateEnd();
                cluster.speciesNum = cluster.start.size();
                removeImitateStart = true;
            }
            else {
                // cluster.content contains only one line, and becomes empty after deletion.
                startIsIsolated = true;
                isolatedIndex = index;
            }
        }
    }

    if (startIsIsolated && isolatedIndex >= 0 && isolatedIndex < clusters.size()) {
        clusters.erase(clusters.begin() + isolatedIndex);
        removeImitateStart = true;
    }

    return removeImitateStart;
}

std::vector<Cluster> geneindex2clusters(const std::vector<std::vector<int>>& geneindex) {
    std::vector<Cluster> clusters;
    clusters.reserve(geneindex.size());  // Pre-allocate memory to avoid multiple expansions.

    for (const auto& line : geneindex) {
        clusters.emplace_back(line);  // Calling the Cluster(vector<int>) constructor
    }

    return clusters;
}

std::vector<Cluster> mergeClusters(const std::vector<Cluster>& clusters) {
    int clustersNum = clusters.size();
    Graph G;

    // build graph
    for (int node = 0; node < clustersNum; ++node) {
        G.addNode(node);
    }
    for (int i = 0; i < clustersNum; ++i) {
        for (int j = 0; j < clustersNum; ++j) {
            if (i != j && clusters[i].isLinked(clusters[j])) {
                G.addEdge(i, j);
            }
        }
    }

    std::vector<Cluster> mergedClusters;
    std::unordered_set<int> visited;

    // First process the nodes where in_degree == 1 or out_degree == 1.
    for (int node = 0; node < clustersNum; ++node) {
        if (visited.count(node)) continue;

        int in_degree = G.inDegree(node);
        int out_degree = G.outDegree(node);

        if (in_degree == 1) {
            int pred_node = G.PrevNodes(node)[0];
            if (visited.count(pred_node)) continue;

            Cluster mergedCluster = clusters[pred_node] + clusters[node];
            mergedClusters.push_back(mergedCluster);

            visited.insert(node);
            visited.insert(pred_node);
        }
        else if (out_degree == 1) {
            int next_node = G.nextNodes(node)[0];
            if (visited.count(next_node)) continue;

            Cluster mergedCluster = clusters[node] + clusters[next_node];
            mergedClusters.push_back(mergedCluster);

            visited.insert(node);
            visited.insert(next_node);
        }
    }

    // Then add the unmerged nodes separately.
    for (int node = 0; node < clustersNum; ++node) {
        if (!visited.count(node)) {
            mergedClusters.push_back(clusters[node]);
            visited.insert(node);
        }
    }

    return mergedClusters;
}


void writeClusterList(const std::vector<Cluster>& clusters, const std::string& outFile) {
    std::ofstream f(outFile);
    if (!f.is_open()) {
        throw std::runtime_error("Failed to open file: " + outFile);
    }

    for (size_t index = 0; index < clusters.size(); ++index) {
        f << "#Cluster-" << index << ":\n";
        for (const auto& line : clusters[index].content) {
            for (size_t i = 0; i < line.size(); ++i) {
                if (line[i] == -9) {
                    f << "-";
                }
                else {
                    f << line[i];
                }
                if (i + 1 < line.size()) {
                    f << "\t";
                }
            }
            f << "\n";
        }
    }

    f.close();
}


int main(int argc, char* argv[]) {
	std::string infile;
	std::string outfile;
	
    // // Output parameters All parameters
    // for (size_t i = 0; i < argc-1; i++){
    //     std::cout << argv[i] << " ";
    // }
    // std::cout << argv[argc-1] << std::endl;


	// Parsing command line arguments
	bool hasInput = false, hasOutput = false;
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-i" && i + 1 < argc) {
			infile = argv[++i];
			hasInput = true;
		}
		else if (arg == "-o" && i + 1 < argc) {
			outfile = argv[++i];
			hasOutput = true;
		}
		else {
			std::cout << "error: unknown argument: " << arg << std::endl;
			std::cout<<"usage: cluster -i file.num -o file.cluster" << std::endl;
			return 1;
		}
	}

	// Check required parameters
	if (!hasInput || !hasOutput) {
		std::cout << "error: must provide (-i and -o argument)" << std::endl;
		std::cout << "usage: cluster -i file.num -o file.cluster" << std::endl;
		return 1;
	}

    // Record start time
    auto start_Time = high_resolution_clock::now();

    std::cout<<"read gene index file..."<<std::endl;
    std::vector<std::vector<int>> geneIndex = readGeneIndex(infile);

    std::cout<<"Get geneindex start..."<<std::endl;
	std::vector<int> start = imitateStart(geneIndex);

    std::cout << "Adding imitate start to geneindex..." << std::endl;
    geneIndex.insert(geneIndex.begin(), start);

    std::cout << "Clustering started..."<< std::endl;
    std::vector<Cluster> lastClusterList= geneindex2clusters(geneIndex);
    size_t lastClusterNum = lastClusterList.size();

    //Start iteration
    std::cout<<"initial clustering... "<< lastClusterNum<<" clusters (containing imitate start)."<<std::endl;
    auto iterStartTime = high_resolution_clock::now();
    std::vector<Cluster> newClusterList = mergeClusters(lastClusterList);
    auto iterEndTime = high_resolution_clock::now();
    auto iterDuration = duration_cast<seconds>(iterEndTime - iterStartTime).count();
    int m_hours = iterDuration / 3600;
    int m_minutes = (iterDuration % 3600) / 60;
    int m_seconds = iterDuration % 60;
    std::cout << "initial clustering "<< "Cluster number " << lastClusterNum<< " -> "<< newClusterList.size()<< ". " << "Runtime: " << m_hours << "h "<< m_minutes << "m "<< m_seconds << "s\n";

    size_t newClusterNum = newClusterList.size();

    size_t iteration = 1;
    while (lastClusterNum != newClusterNum){
        std::cout << "Iteration " << iteration << ": Cluster number " << lastClusterNum << " -> " << newClusterNum << ". Merging clusters ...";
        iteration++;
        lastClusterNum = newClusterNum;

        auto iterStartTime = high_resolution_clock::now();
        newClusterList = mergeClusters(newClusterList);
        iterEndTime = high_resolution_clock::now();
        iterDuration = duration_cast<seconds>(iterEndTime - iterStartTime).count();
        m_hours = iterDuration / 3600;
        m_minutes = (iterDuration % 3600) / 60;
        m_seconds = iterDuration % 60;
        std::cout << "Runtime: "<< m_hours << "h "<< m_minutes << "m "<< m_seconds << "s\n";
        
        newClusterNum = newClusterList.size();
    }
    std::cout << "Clustering completed. " << newClusterNum << " clusters." << std::endl;

    std::cout << "Removing imitate start from clusters..." << std::endl;
    bool isRemoved = removeImitateStartFromClusters(newClusterList, start);
    if (isRemoved){
        std::cout << "Imitate start removed from clusters." << std::endl;
    }
    else {
        std::cout<<"error:No imitate start found in clusters!!!"<<std::endl;
    }

    writeClusterList(newClusterList, outfile);
    std::cout << "Clusters written to " << outfile << "." << std::endl;

    // Record end time
    auto end_Time = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end_Time - start_Time).count();
    m_hours = duration / 3600;
    m_minutes = (duration % 3600) / 60;
    m_seconds = duration % 60;
    std::cout << "Runtime: "<< m_hours << "h "<< m_minutes << "m "<< m_seconds << "s\n";

    std::cout << "Done." << std::endl;

    return 0;
}