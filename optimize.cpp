#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>
#include <numeric>

using namespace std;


// 工件结构体，存储某工件在不同机器中的加工时间
struct Job { 
   vector<int> processingTimes;
};

// 对指定序号的测试用例执行优化
void startOpt(int testCaseNumber);

//通用工具类函数
void   generateRandomSchedule(vector<int>& solution);// 生成随机调度顺序
int    calculateMakespan(const vector<Job>& jobs, const vector<int>& schedule, int m);// 计算特定调度顺序花费的时间

//模拟退火相关函数
vector<int> simulatedAnnealing(vector<Job>& jobs, int m, double temp, double coolingRate, int maxIterations);// 模拟退火算法主程序

//遗传算法相关函数
vector<vector<int>>            initializePopulation(int populationSize, int numJobs);// 初始化种群
vector<pair<int, vector<int>>> evaluatePopulation(const vector<Job>& jobs, const vector<vector<int>>& population, int m);// 评估种群
vector<vector<int>>            rouletteWheelSelection(const vector<pair<int, vector<int>>>& evaluatedPopulation, int numToSelect);// 使用轮盘赌选择策略进行选择
vector<int>                    crossover(const vector<int>& parent1, const vector<int>& parent2);// 交叉操作
void                           mutate(vector<int>& individual, double mutationRate);// 突变操作
vector<int>                    geneticAlgorithm(const vector<Job>& jobs, int m, int populationSize, int generations, double mutationRate);// 遗传算法主程序

//输入输出函数
vector<Job> readJobsFromFile(const string& filename, int instanceNumber, int& m);// 从测试用例文件中读取工件信息
void printJobsInfo(const vector<Job>& jobs);// 打印工件信息
void printResult(const vector<int>& schedule, const vector<Job>& jobs, int m, int makespan);// 打印优化结果
//void printSchedule(const vector<int>& schedule);// 打印调度顺序

//初始化随机数，使用dis(gen)可得到0到1之间的随机数
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0,1);

int main() {
    srand(time(NULL));  // 初始化随机种子
    int testCaseNumber = 0;
    // cout << "请输入测试用例的序号(从0到10之间选一个): ";
    // cin >> testCaseNumber;
    startOpt(testCaseNumber);
    // for(int i = 0; i < 11; i++){
    //     startOpt(i);
    // }
    return 0;
}

void startOpt(int testCaseNumber){
    string filePath = "testcase.txt";// 指定测试用例的文件路径
    int m;// 机器数
    vector<Job> jobs = readJobsFromFile(filePath, testCaseNumber, m);
    //printJobsInfo(jobs);// 打印工件信息（用来测试一下输入是否正确）

    //模拟退火算法
    double temp = 20000;// 初始温度
    double coolingRate = 0.8;// 冷却率
    int maxIterations = 3000;// 最大迭代次数
    vector<int> bestSolutionSA = simulatedAnnealing(jobs, m, temp, coolingRate, maxIterations);
    int bestMakespanSA = calculateMakespan(jobs, bestSolutionSA, m);
    cout << "模拟退火：" << endl;
    printResult(bestSolutionSA, jobs, m, bestMakespanSA);

    //遗传算法
    int populationSize = 200;// 种群大小
    int generations = 300;// 迭代次数
    double mutationRate = 0.05;// 突变率
    vector<int> bestSolutionGA = geneticAlgorithm(jobs, m, populationSize, generations, mutationRate);
    int bestMakespanGA = calculateMakespan(jobs, bestSolutionGA, m);
    cout << "遗传算法：" << endl;
    printResult(bestSolutionGA, jobs, m, bestMakespanGA);
}

void generateRandomSchedule(vector<int>& schedule) {
    // 先将schedule初始化为 {0, 1, 2, ..., n-1}
    for (int i = 0; i < schedule.size(); i++) {
        schedule[i] = i;
    }
    // 将schedule{0, 1, 2, ..., n-1}内的顺序进行随机打乱
    shuffle(schedule.begin(), schedule.end(), mt19937(random_device()()));
}

int calculateMakespan(const vector<Job>& jobs, const vector<int>& schedule, int m) {
    int n = schedule.size();
    vector<vector<int>> completion_time(n+1, vector<int>(m+1, 0));
    // 为什么是m+1和n+1呢？这样是为了方便对边界进行计算，第一行和第一列都是0，而计算是从第二行第二列开始进行的，这样就不用特殊处理了
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            completion_time[i][j] = jobs[schedule[i-1]].processingTimes[j-1] +
                                    max(completion_time[i-1][j], completion_time[i][j-1]);
        }
    }
    return completion_time[n][m];
}

vector<int> simulatedAnnealing(vector<Job>& jobs, int m, double temp, double coolingRate, int maxIterations) {
    // temp是初始温度，coolingRate是冷却率，maxIterations是最大迭代次数
    vector<int> schedule(jobs.size());  // 用来存放工件的调度顺序
    generateRandomSchedule(schedule);// 生成初始解

    vector<int> currentSchedule = schedule;
    int minMakespan = calculateMakespan(jobs, currentSchedule, m);

    for (int i = 0; i < maxIterations; i++) {
        // 随机选择两个位置交换
        int pos1 = rand() % schedule.size();
        int pos2 = rand() % schedule.size();
        swap(currentSchedule[pos1], currentSchedule[pos2]);

        int currentMakespan = calculateMakespan(jobs, currentSchedule, m);
        int energyChange = currentMakespan - minMakespan;

        // 如果新解更好，或者根据概率接受更差的解
        if (energyChange < 0 || exp(-energyChange / temp) > (double)rand() / RAND_MAX) {
            minMakespan = currentMakespan;
            schedule = currentSchedule;  // 更新最优解
        } else {
            swap(currentSchedule[pos1], currentSchedule[pos2]);  // 恢复原状
        }

        temp *= coolingRate;  // 降低温度
    }
    return schedule;
}

vector<vector<int>> initializePopulation(int populationSize, int numJobs) {
    vector<vector<int>> population;// 初始化种群
    for (int i = 0; i < populationSize; ++i) {// populationSize对应种群内的个体数量
        vector<int> individual(numJobs);// numJobs对应每个个体内工件的数量
        generateRandomSchedule(individual);
        population.push_back(individual);
    }
    return population;
}

vector<pair<int, vector<int>>> evaluatePopulation(const vector<Job>& jobs, const vector<vector<int>>& population, int m) {
    vector<pair<int, vector<int>>> evaluated;
    for (const auto& individual : population) {
        int makespan = calculateMakespan(jobs, individual, m);// 每个个体的完工时间
        evaluated.push_back({makespan, individual});
    }
    sort(evaluated.begin(), evaluated.end());// 按照完工时间从小到大排序
    return evaluated;
}

vector<vector<int>> rouletteWheelSelection(const vector<pair<int, vector<int>>>& evaluatedPopulation, int numToSelect) {
    vector<vector<int>> selected;
    double totalFitness = 0.0;
    vector<double> probabilities;

    // 计算总适应度
    for (const auto& p : evaluatedPopulation) {
        // 设适应度是完工时间的倒数，完工时间越小，适应度越高
        totalFitness += 1.0 / p.first; 
    }

    // 计算选择概率
    for (const auto& p : evaluatedPopulation) {
        probabilities.push_back((1.0 / p.first) / totalFitness);
    }

    // 累积概率
    partial_sum(probabilities.begin(), probabilities.end(), probabilities.begin());// p={p1, p1+p2, p1+p2+p3, ...}

    // 随机选择个体

    for (int i = 0; i < numToSelect; ++i) {
        double randProb = dis(gen);
        int idx = 0; // 初始化索引
        // 循环遍历累积概率数组
        for (; idx < probabilities.size(); ++idx) {
            if (probabilities[idx] >= randProb) {
                break; // 如果找到第一个大于等于randProb的概率，停止搜索
            }
        }
        // 添加选中的个体到selected数组
        selected.push_back(evaluatedPopulation[idx].second);
    }

    return selected;
}

vector<int> crossover(const vector<int>& parent1, const vector<int>& parent2) {
    vector<int> child(parent1.size());
    int crossPoint = rand() % (parent1.size() - 2) + 1; // 保证交叉点不在末端
    copy_n(parent1.begin(), crossPoint, child.begin());
    copy_if(parent2.begin(), parent2.end(), child.begin() + crossPoint, [&child](int gene) {
        return find(child.begin(), child.end(), gene) == child.end();
    });
    return child;
}


void mutate(vector<int>& individual, double mutationRate = 0.05) {
    double chance = dis(gen);
    if (chance < mutationRate) {
        int index1 = rand() % individual.size();
        int index2 = rand() % individual.size();
        swap(individual[index1], individual[index2]);
    }
}

// 遗传算法主程序
vector<int> geneticAlgorithm(const vector<Job>& jobs, int m, int populationSize, int generations, double mutationRate) {
    vector<vector<int>> population = initializePopulation(populationSize, jobs.size());// 初始化种群
    for (int gen = 0; gen < generations; ++gen) {
        auto evaluated = evaluatePopulation(jobs, population, m);// 评估种群
        vector<vector<int>> parents = rouletteWheelSelection(evaluated, populationSize / 2);// 选择父代
        vector<vector<int>> newGeneration;// 新一代种群

        while (newGeneration.size() < populationSize) {
            int idx1 = rand() % parents.size();// 随机选择两个父代
            int idx2 = rand() % parents.size();
            vector<int> child = crossover(parents[idx1], parents[idx2]);// 交叉操作
            mutate(child, mutationRate);// 突变操作
            newGeneration.push_back(child);// 添加子代到新一代种群
        }

        population = newGeneration;// 更新种群
    }

    auto finalEvaluated = evaluatePopulation(jobs, population, m);
    return finalEvaluated.front().second;
}

vector<Job> readJobsFromFile(const string& filename, int instanceNumber, int& m) {
    ifstream file(filename);
    string line;
    vector<Job> jobs;
    bool startReading = false;  // 控制是否开始读取数据的标志
    int n;  // 工件数
    if (file.is_open()) {
        while (getline(file, line)) {
            // 检查是否到达指定的实例
            if (!startReading && line.find("instance " + to_string(instanceNumber)) != string::npos) {
                startReading = true;  // 开始读取数据
                getline(file, line);  // 读取下一行，即包含工件数和机器数的行
                istringstream iss(line);
                iss >> n >> m;
                continue;  // 跳过当前循环，进入数据读取部分
            }

            if (startReading) {
                if (n == 0) {
                    break;  // 读取完n个工件的数据后退出
                }
                Job job;
                job.processingTimes.resize(m);
                istringstream lineStream(line);
                int machineIndex, time;
                
                for (int i = 0; i < m; ++i) {
                    lineStream >> machineIndex >> time;
                    job.processingTimes[machineIndex] = time;
                }
                jobs.push_back(job);
                n--;
            }
        }
        file.close();
    } else {
        cerr << "测试用例文件异常" << endl;
    }
    
    return jobs;
}

void printJobsInfo(const vector<Job>& jobs) {
    int jobIndex = 0;
    for (const auto& job : jobs) {
        cout << "Job " << jobIndex << ": ";
        for (size_t i = 0; i < job.processingTimes.size(); ++i) {
            cout << "Machine " << i << " Time " << job.processingTimes[i] << "; ";
        }
        cout << endl;
        ++jobIndex;
    }
}

void printResult(const vector<int>& schedule, const vector<Job>& jobs, int m, int makespan) {
    std::cout << "调度顺序：";
    for (size_t i = 0; i < schedule.size()-1; i++) {
        std::cout << schedule[i]+1 << ", ";
    }
    std::cout << schedule[schedule.size()-1]+1 << "";
    std::cout << std::endl;
    cout << "所需加工时间: " << makespan << endl;
}

// void printSchedule(const vector<int>& schedule) {
//     std::cout << "目前的调度顺序是：";
//     for (size_t i = 0; i < schedule.size()-1; i++) {
//         std::cout << schedule[i]+1 << ", ";
//     }
//     std::cout << schedule[schedule.size()-1]+1 << "";
//     std::cout << std::endl;
// }