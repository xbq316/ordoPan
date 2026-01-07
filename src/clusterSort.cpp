#include <vector>
#include <sstream>
#include <climits>
#include <stdexcept> 
#include <utility>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <random>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <queue>
#include <cfloat>



using namespace std::chrono;
class Cluster{
    public:
        size_t geneNum;    // The number of genes (rows) is not the actual number of genes; it's just for ease of calculation.
        size_t speciesNum; // Number of species (number of columns)
        std::vector<int> start; // The first non--9 value in each column
        std::vector<int> end;   // The last non--9 value in each column

        Cluster(const std::vector<std::vector<int>>& cluster)
        : geneNum(cluster.size()),
        speciesNum(cluster.empty() ? 0 : cluster[0].size()),
        start(std::vector<int>(speciesNum, -9)),  // Directly initialize start and end
        end(std::vector<int>(speciesNum, -9)) {
            // Traverse each column
            for (size_t i = 0; i < speciesNum; i++) {
                // Find the first non--9 value in each column (from front to back).
                for (size_t j = 0; j < geneNum; j++) {
                    if (cluster[j][i] != -9) {
                        start[i] = cluster[j][i];
                        break;
                    }
                }

                // Find the last non--9 value in each column (from right to left).
                for (size_t j = geneNum; j-- > 0; ) {
                    if (cluster[j][i] != -9) {
                        end[i] = cluster[j][i];
                        break;
                    }
                }
            }
        }

        Cluster(const std::vector<int>& line):geneNum(1),speciesNum(line.size()),start(line),end(line){}
        
        Cluster() : geneNum(0), speciesNum(0), start(), end() {}// Default constructor

        Cluster operator+(const Cluster& other) const {
            if (speciesNum != other.speciesNum) {
                throw std::invalid_argument("SpeciesNum mismatch in Cluster operator+");
            }
            std::vector<int> newStart(speciesNum,-9);
            std::vector<int> newEnd(speciesNum,-9);

            for (size_t i = 0; i < speciesNum; ++i) {
                if (this->start[i] == -9) newStart[i] = other.start[i];
                else {newStart[i] = this->start[i];}

                if (other.end[i] == -9) newEnd[i] = this->end[i];
                else {newEnd[i] = other.end[i];}
            }
            Cluster newCluster;
            newCluster.start = newStart;
            newCluster.end = newEnd;
            newCluster.geneNum = this->geneNum + other.geneNum;
            newCluster.speciesNum = this->speciesNum; 
            return newCluster;
        }

        float averageRank() const {
            float result = 0;
            for (size_t i = 0; i < this->speciesNum; i++){
                if (this->start[i] == -9){continue;}
                result += this->start[i];
            }
            return result / this->supportSpecies();
        }

        float linkNum(const Cluster& ortherCluster) const {
            float result = 0;
            for (size_t i = 0; i < this->speciesNum; i++){
                if (this->end[i]==-9 || ortherCluster.start[i]==-9){continue;}
                else if (this->end[i] +1 == ortherCluster.start[i]){ result++;}
            }
            return result;
        }

        float linkScore(const Cluster& ortherCluster) const {
            float result = 0;
            for (size_t i = 0; i < this->speciesNum; i++){
                if (this->end[i] == -9 || ortherCluster.start[i] == -9) { continue; }
                int diff = ortherCluster.start[i] - this->end[i];
                if (diff < 0) { continue; }
                else{result += exp(1-diff);}
            }
            return result;
        }

        size_t supportSpecies() const {
            size_t count = 0;
            for (int val : this->start) {
                count += (val != -9); // Implicitly converts to 0 or 1
            }
            return count;
        }

        bool islink(const Cluster& other) const {
            for (size_t i = 0; i < speciesNum; i++) {
                if (end[i] != -9 && other.start[i] != -9 && end[i] + 1 == other.start[i]) {
                    return true; // There is a connection
                }
            }
            return false; // No connection
        }
};

//Generate initial solution
std::vector<int> linkcluster2solution(const std::vector<Cluster> clustervector){
    size_t clustervectorSize=clustervector.size();
    std::vector<bool> visited(clustervectorSize,false);
    std::vector<int> solution={0};
    visited[0]=true;

    Cluster currentCluster;
    currentCluster=clustervector[0];
    bool flag=true;

    while (flag){
        float bestScore=0;
        size_t bestindex=-1;

        for (size_t i = 1; i < clustervectorSize; i++){
            if (visited[i]){continue;}
            float thisScore=currentCluster.linkScore(clustervector[i]);
            if (thisScore>bestScore){bestScore=thisScore;bestindex=i;}
        }
        if(bestScore>=1 && bestindex!=-1){
            visited[bestindex]=true;
            solution.push_back(bestindex);
            currentCluster=currentCluster+clustervector[bestindex];
        }
        else{flag=false;}
    }
    
    std::vector<int> remainVector ;
    for (size_t i = 1; i < visited.size(); i++){
        if (visited[i]==false){remainVector.push_back(i);}
    }
    
    std::stable_sort(remainVector.begin(), remainVector.end(), [&](int a, int b) {
        if (clustervector[a].supportSpecies() != clustervector[b].supportSpecies())
            return clustervector[a].supportSpecies() > clustervector[b].supportSpecies();
        return clustervector[a].geneNum > clustervector[b].geneNum;
    });

    if(solution.size()/(remainVector.size()+solution.size())<0.7){
        std::vector<int> result_solution(clustervectorSize);
        std::iota(result_solution.begin(), result_solution.end(), 0);

        std::stable_sort(result_solution.begin(), result_solution.end(), [&](int a, int b) {
            return clustervector[a].averageRank() < clustervector[b].averageRank();
        });

        return result_solution;
    }

    std::vector<int> iterVector=std::move(remainVector);
    remainVector.clear();

    while (true){
        // Start Insertion Sort
        for (const int &i : iterVector){
            float bestScore = 0;
            size_t bestPos = -1; 
            currentCluster=clustervector[0];
            for (size_t j = 0; j < solution.size(); j++){
                if (j!=0){currentCluster=currentCluster + clustervector[solution[j]];}
                float thisScore = currentCluster.linkScore(clustervector[i]);
                if (thisScore > bestScore) {
                    bestScore = thisScore;
                    bestPos = j + 1; // Insert after j
                }
            }
            if (bestScore<1 && bestPos==-1){remainVector.push_back(i);}
            else{solution.insert(solution.begin() + bestPos, i);}
            
        }
        if (remainVector.empty()){break;}
        else{
            iterVector=std::move(remainVector);
            remainVector.clear();
        }
    }
    
    return solution;
}

std::vector<int> rankSort(const std::vector<Cluster>& clusterVector) {
    std::vector<int> result(clusterVector.size());
    std::iota(result.begin(), result.end(), 0); // Initialize to 0, 1, ..., n-1
    // Sort by averageRank
    std::stable_sort(result.begin(), result.end(), [&](int a, int b) {
        return clusterVector[a].averageRank() < clusterVector[b].averageRank();
    });
    return result;
}

// Path length
float solution_score(const std::vector<int>& solution, const std::vector<Cluster>& clusterVector) {
    float total = 0.0;
    Cluster currentCluster = clusterVector[solution[0]];
    for (size_t i = 1; i < solution.size(); ++i) {
        total+= currentCluster.linkScore(clusterVector[solution[i]]);
        currentCluster = currentCluster + clusterVector[solution[i]];
    }
    return total;
}

bool isNormMatrix(const std::vector<std::vector<std::vector<int>>> &clusterVector) {
    size_t speciesNum;
    if (clusterVector.empty()){return false ;}
    if (clusterVector[0].empty()){return false ;}
    if (clusterVector[0][0].empty()){return false ;}
    speciesNum = clusterVector[0][0].size();
    
    for (int i = 0; i < clusterVector.size(); i++) {
        for (int j = 0; j < clusterVector[i].size(); j++) {
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

Cluster imitateStart(const std::vector<Cluster>& clusterVector) {
    if (clusterVector.empty()) { std::cerr << "clusterVector is empty !"; }
    const size_t& speciesNum = clusterVector[0].speciesNum;
    std::vector<int> line(speciesNum, -9);
    size_t clusterNum = clusterVector.size();
    for (size_t thisSpecies = 0; thisSpecies < speciesNum; thisSpecies++){
        int min_start = INT_MAX;
        for (size_t thisCluster = 0; thisCluster < clusterNum; thisCluster++){
            const int& this_start = clusterVector[thisCluster].start[thisSpecies];
            if (this_start == -9) { continue; }
            if (this_start < min_start) { min_start = this_start; }
        }
        if (min_start == INT_MAX) { continue; }
        else { line[thisSpecies] = min_start -1; }
    }
    return Cluster(line);
}

Cluster imitateEnd(const std::vector<Cluster>& clusterVector) {
    if (clusterVector.empty()) { std::cerr << "clusterVector is empty !"; }
    const size_t& speciesNum = clusterVector[0].speciesNum;
    std::vector<int> line(speciesNum, -9);
    size_t clusterNum = clusterVector.size();
    for (size_t thisSpecies = 0; thisSpecies < speciesNum; thisSpecies++) {
        int max_end = 0;  
        for (size_t thisCluster = 0; thisCluster < clusterNum; thisCluster++) {
            const int& this_end = clusterVector[thisCluster].end[thisSpecies];  // Use end instead of start
            if (this_end == -9) { continue; }
            if (this_end > max_end) { max_end = this_end; }  // The comparison directions are opposite
        }
        if (max_end == 0) { continue; }
        else { line[thisSpecies] = max_end + 1; }  // The ending position is usually +1
    }
    return Cluster(line);
}


void writeClusterVector(const std::vector<std::vector<std::vector<int>>>& clusterVector, const std::string& outFile) {
    std::ofstream out(outFile);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file " << outFile << " for writing." << std::endl;
        return;
    }

    // Traversing a 3D vector
    for (size_t i = 0; i < clusterVector.size(); i++) {
        for (size_t j = 0; j < clusterVector[i].size(); j++) {
            // Traverse the innermost vector (numerical value).
            for (size_t k = 0; k < clusterVector[i][j].size(); k++) {
                const int& value = clusterVector[i][j][k];
                out << (value == -9 ? "-" : std::to_string(value));
                if (k != clusterVector[i][j].size() - 1) {
                    out << "\t";
                }
            }
            out << "\n";  // End the current subcluster
        }
    }
    out.close();
}



std::vector<std::vector<float>> create_Matrix(const std::vector<Cluster> &clusterVector) {
    std::cout << "Calculate the distance matrix..." << std::endl;
    int n = clusterVector.size();
    std::vector<std::vector<float>> matrix(n, std::vector<float>(n, 0));

    // Calculate the distance matrix directly in sequence
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                matrix[i][j] = clusterVector[i].linkScore(clusterVector[j]);
            }
        }
    }
    return matrix;
}


std::vector<std::vector<float>> readMatrixe(const std::string& filename) {
    std::vector<std::vector<float>> matrix;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return matrix;
    }

    std::string line;
    while (getline(file, line)) {
        std::vector<float> row;
        std::istringstream iss(line);
        float num;

        while (iss >> num) {
            row.push_back(num);
        }

        if (!row.empty()) {
            matrix.push_back(row);
        }
    }

    file.close();
    return matrix;
}


void write_Matrix(const std::vector<std::vector<float>>& distance_matrix, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    for (const auto& row : distance_matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i];
            if (i < row.size() - 1) {
                outfile << "\t";
            }
        }
        outfile << "\n";
    }
    outfile.close();
}


// Delete node
std::pair<std::vector<int>, std::vector<int>> destroy_operator(const std::vector<int>& solution,const float & removal_rate, std::mt19937& rng) {
    std::vector<int> removable(solution.begin() + 1, solution.end());
    int num_remove = std::max(1, int(removal_rate * removable.size()));
    std::shuffle(removable.begin(), removable.end(), rng );
    std::vector<int> removed(removable.begin(), removable.begin() + num_remove);

    std::set<int> removed_set(removed.begin(), removed.end());
    std::vector<int> partial;
    for (int node : solution) {
        if (!removed_set.count(node)) {
            partial.push_back(node);
        }
    }
    return { partial, removed };
}

// Insert node
std::vector<int> repair_operator(std::vector<int> partial,
    const std::vector<int>& removed,
    const std::vector<std::vector<float>>& scoreMatrix,
    std::mt19937& rng
) {
    for (int node : removed) {
        float max_cost = 0;
        int best_pos = 1;
        for (size_t i = 0; i < partial.size() - 1; ++i) {
            float cost = scoreMatrix[partial[i]][node];// +scoreMatrix[node][partial[i + 1]] - scoreMatrix[partial[i]][partial[i + 1]];
            if (cost > max_cost) {
                max_cost = cost;
                best_pos = i + 1;
            }
        }
        partial.insert(partial.begin() + best_pos, node);
    }
    return partial;
}


// Main Algorithm
std::vector<int> LNS_TSP(const std::vector<std::vector<float>>& scoreMatrix,
    int iterations,
    std::mt19937& rng,
    const int & start,
    const std::vector<int> & initSolution,
    const std::vector<Cluster> & clusterVector,
    const float& rate
) {
    auto start_time = steady_clock::now();

    std::vector<int> current_solution = initSolution;
    float current_cost = solution_score(current_solution, clusterVector);
    std::vector<int> best_solution = current_solution;
    float best_cost = current_cost;

    std::uniform_int_distribution<int> dist(0, 9);

    std::cout<<"Start Iteration..." << std::endl;
    for (int iter = 1; iter <= iterations; ++iter) {
        auto result = destroy_operator(current_solution, rate, rng);
        std::vector<int> partial = result.first;
        std::vector<int> removed = result.second;

        std::vector<int> new_solution = repair_operator(partial, removed, scoreMatrix, rng);
        float new_cost = solution_score(new_solution, clusterVector);

        if (new_cost > best_cost) {
            best_cost = new_cost;
            best_solution = new_solution;
            current_solution = new_solution;
        }
        else if (dist(rng) < 1) {// There is a 10% probability of accepting the difference solution.
            current_solution = new_solution;
        }

        if (iter % 10 == 0 || iter == iterations) {
            auto now = steady_clock::now();
            double elapsed_sec = duration_cast<seconds>(now - start_time).count();
            double rate = double(iter) / elapsed_sec;
            int remaining = int((iterations - iter) / rate);
            std::cout << std::fixed << std::setprecision(2)
                << "\r[Iter: " << iter << "/" << iterations << "] Best: " << best_cost
                << " | Time: " << elapsed_sec << "s"
                << " | ETA: " << remaining << "s" <<
                "\r" << std::flush;
        }

        std::cout << "working-route ";
        if (iter % 10000 == 0 || iter == iterations) {
            for (size_t idx = 0; idx < new_solution.size(); ++idx) {
                size_t clusterIdx = new_solution[idx];
                if (clusterIdx == 0) { continue; }
                else{ clusterIdx -= 1; }
                std::cout << clusterIdx << " ";
            }
            std::cout << std::endl;
        }
    }

    return best_solution;
}

void reportBug(const std::string& programName,const std::string& message) {
	std::cerr << message << std::endl;
	std::cerr << "Usage: " << programName << " [-i input_file] [-o output_file] [-s seed(default: 42)] [-n iterations(default: 50, unit: 10 thousand)] [-r rate(default:0.1)]" << std::endl;
}

int main(int argc, char* argv[]) {
	std::string inFile;
	std::string outFile;
	size_t iteration=1;
	float rate=0.1;
	unsigned int seedValue=42;
	

	// Parsing command line arguments
	bool hasInput = false, hasOutput = false;
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-i" && i + 1 < argc) {
			inFile = argv[++i];
			hasInput = true;
		}
		else if (arg == "-o" && i + 1 < argc) {
			outFile = argv[++i];
			hasOutput = true;
		}
		else if (arg == "-s" && i + 1 < argc) {
			seedValue = std::atoi(argv[++i]);
		}
		else if (arg == "-n" && i + 1 < argc) {
			iteration = std::atoi(argv[++i]);
		}
		else if (arg == "-r" && i + 1 < argc) {
			rate = std::atof(argv[++i]);
		}
		else {
			reportBug(argv[0],"Unknown option or missing argument: " + arg);
			return 1;
		}
	}

	// Check required parameters
	if (!hasInput || !hasOutput) {
		reportBug(argv[0], "error: must provide (-i and -o argument)");
		return 1;
	}

    iteration = iteration * 10000;

	// Program Logic...
	std::cout << "input file: " << inFile << std::endl;
	std::cout << "output file: " << outFile << std::endl;
	std::cout << "random seed: " << seedValue << std::endl;
	std::cout << "iterations num : " << iteration << std::endl;
	std::cout << "rate: " << rate << std::endl;


	std::mt19937 seed(seedValue);

	if (rate < 0 || rate>1) {
		reportBug(argv[0], "error: rate must be between 0 and 1");
		return 1;
	}

	std::vector<std::vector<std::vector<int>>> clusterMatrix= readCluster(inFile);
    if (!isNormMatrix(clusterMatrix)){
        std::cout<<"error: input file is not norm matrix"<<std::endl;
        return 1;
    }
	size_t clusterNum=clusterMatrix.size();
	std::vector<Cluster> clusterVector;

	for (size_t i=0;i<clusterNum;i++){
		clusterVector.emplace_back(Cluster(clusterMatrix[i]));
	}

	Cluster imitateStartCluster=imitateStart(clusterVector);
	clusterVector.insert(clusterVector.begin(),imitateStartCluster);//Introducing a virtual starting point
	
	std::vector<std::vector<float>> matrix= create_Matrix(clusterVector);
	
    std::vector<int> initSolution = linkcluster2solution(clusterVector);//rankSort(clusterVector);

	std::vector<int> bestRoute = LNS_TSP(matrix,iteration,seed,0, initSolution,clusterVector,rate);

    std::cout << "best-route ";
	std::vector<std::vector<std::vector<int>>> bestClusterVector;
	for (size_t idx = 0; idx < bestRoute.size(); ++idx) {
		size_t clusterIdx = bestRoute[idx];
		if (clusterIdx == 0) { continue; }
		else{ clusterIdx -= 1; }
        std::cout << clusterIdx << " ";
		bestClusterVector.emplace_back(clusterMatrix[clusterIdx]);
	}
    std::cout<<std::endl;

	writeClusterVector(bestClusterVector,outFile);

	return 0;
}
