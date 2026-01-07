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


class Graph {
private:
    // 邻接表
    std::unordered_map<int, std::vector<int>> outEdges; // 出边
    std::unordered_map<int, std::vector<int>> inEdges;  // 入边

public:
    void addNode(int u) {
        outEdges[u];  // 确保 key 存在
        inEdges[u];
    }

    // 添加边 u -> v
    void addEdge(int u, int v) {
        outEdges[u].push_back(v);
        inEdges[v].push_back(u);
    }

    // 返回出边节点列表
    const std::vector<int>& nextNodes(int u) const {
        static std::vector<int> empty;
        auto it = outEdges.find(u);
        if (it != outEdges.end()) return it->second;
        return empty;
    }

    // 返回入边节点列表
    const std::vector<int>& PrevNodes(int u) const {
        static std::vector<int> empty;
        auto it = inEdges.find(u);
        if (it != inEdges.end()) return it->second;
        return empty;
    }

    // 返回出度
    int outDegree(int u) const {
        auto it = outEdges.find(u);
        return it != outEdges.end() ? it->second.size() : 0;
    }

    // 返回入度
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


    // 构造函数：单行初始化
    Cluster(const std::vector<int>& line) {
        start = line;
        end = line;
        speciesNum = line.size();
        content.push_back(line);
    }

    // 从多行初始化
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

    // 更新 start
    std::vector<int> updateStart() const {
        std::vector<int> newLine = content[0];
        for (const auto& line : content) {
            for (size_t i = 0; i < newLine.size(); ++i) {
                if (newLine[i] == -9) newLine[i] = line[i];
            }
        }
        return newLine;
    }

    // 更新 end
    std::vector<int> updateEnd() const {
        std::vector<int> newLine = content.back();
        for (auto it = content.rbegin(); it != content.rend(); ++it) {
            for (size_t i = 0; i < newLine.size(); ++i) {
                if (newLine[i] == -9) newLine[i] = (*it)[i];
            }
        }
        return newLine;
    }

    // 合并两个 Cluster
    Cluster operator+(const Cluster& nextCluster) const {
        std::vector<std::vector<int>> combined = content;
        combined.insert(combined.end(), nextCluster.content.begin(), nextCluster.content.end());
        return Cluster::fromMultiline(combined);
    }

    // 判断是否相邻
    bool isLinked(const Cluster& nextCluster) const {
        for (size_t i = 0; i < end.size(); ++i) {
            int lastEnd = end[i];
            int thisStart = nextCluster.start[i];
            if (lastEnd == -9 || thisStart == -9) continue;
            if (lastEnd + 1 == thisStart) return true;
        }
        return false;
    }

    // 支持的物种数量
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
        // 去掉首尾空白
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
                cluster.content.erase(cluster.content.begin());  // 移除第一行
                cluster.start = cluster.updateStart();
                cluster.end = cluster.updateEnd();
                cluster.speciesNum = cluster.start.size();
                removeImitateStart = true;
            }
            else {
                // cluster.content 只有一行，被删除后空
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
    clusters.reserve(geneindex.size());  // 预分配内存，避免多次扩容

    for (const auto& line : geneindex) {
        clusters.emplace_back(line);  // 调用 Cluster(vector<int>) 构造函数
    }

    return clusters;
}

std::vector<Cluster> mergeClusters(const std::vector<Cluster>& clusters) {
    int clustersNum = clusters.size();
    Graph G;

    // 建图
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

    // 先处理 in_degree == 1 或 out_degree == 1 的节点
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

    // 再把没合并的节点单独加入
    for (int node = 0; node < clustersNum; ++node) {
        if (!visited.count(node)) {
            mergedClusters.push_back(clusters[node]);
            visited.insert(node);
        }
    }

    return mergedClusters;
}

std::vector<Cluster> initMergeClusters(const std::vector<Cluster>& clusters) {
    int clustersNum = clusters.size();
    Graph G;

    // 建图
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

    // 对当前节点出度为1，下一节点入度为1 的节点进行合并
    // 对当前节点入度为1，前驱节点出度为1 的节点进行合并
    for (int node = 0; node < clustersNum; ++node) {
        if (visited.count(node)) continue;

        std::vector<int> chain;
        chain.push_back(node);

        // 向前收集
        int cur = node;
        while (G.inDegree(cur) == 1) {
            int pred = G.PrevNodes(cur)[0];
            if (visited.count(pred)||G.outDegree(pred) != 1) break;
            visited.insert(pred);
            visited.insert(cur);
            chain.insert(chain.begin(), pred);
            cur = pred;
        }

        // 向后收集
        cur = node;
        while (G.outDegree(cur) == 1) {
            int nxt = G.nextNodes(cur)[0];
            if (visited.count(nxt)||G.inDegree(nxt)!=1) break;
            visited.insert(nxt);
            visited.insert(cur);
            chain.push_back(nxt);
            cur = nxt;
        }

        // 合并链条
        if (chain.size() > 1) {
            Cluster mergedCluster = clusters[chain[0]];
            for (size_t i = 1; i < chain.size(); ++i) {
                mergedCluster = mergedCluster + clusters[chain[i]];
            }
            mergedClusters.push_back(mergedCluster);
        }
    }

    // 再把没合并的节点单独加入
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
	
    // 输出参数所有参数
    for (size_t i = 0; i < argc-1; i++){
        std::cout << argv[i] << " ";
    }
    std::cout << argv[argc-1] << std::endl;


	// 解析命令行参数
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

	// 检查必须参数
	if (!hasInput || !hasOutput) {
		std::cout << "error: must provide (-i and -o argument)" << std::endl;
		std::cout << "usage: cluster -i file.num -o file.cluster" << std::endl;
		return 1;
	}

    // 记录开始时间
    using namespace std::chrono;
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

    //开始迭代
    std::cout<<"initial clustering... "<< lastClusterNum<<" clusters (containing imitate start)."<<std::endl;
    std::vector<Cluster> newClusterList = initMergeClusters(lastClusterList);
    size_t newClusterNum = newClusterList.size();

    size_t iteration = 0;
    while (lastClusterNum != newClusterNum){
        std::cout << "Iteration " << iteration << ": Cluster number " << lastClusterNum << " -> " << newClusterNum << ". Merging clusters ..." << std::endl;
        iteration++;
        lastClusterNum = newClusterNum;
        newClusterList = mergeClusters(newClusterList);
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

    // 记录结束时间
    auto end_Time = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(end_Time - start_Time).count();
    int hours = duration / 3600;
    int minutes = (duration % 3600) / 60;
    int seconds = duration % 60;
    std::cout << "Total time taken: "<< hours << "h "<< minutes << "m "<< seconds << "s\n";

    std::cout << "Done." << std::endl;

    return 0;
}