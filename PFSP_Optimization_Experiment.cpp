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
#include <chrono>
#include <iomanip>

using namespace std;

//前言
/*
代码和报告的实现思路略微有差别, 
区别在于各机器加工各工件所需时间的存储方式: 
在使用数学语言表达时, 我用二维数组的行代表机器, 列代表工件, 这样在分析总加工时长的时候很直观
后来写代码的时候, 发现测试用例文件正好反过来了, 每行代表的是各工件在各机器中对应的时间
为了方便输入, 代码就按照测试用例的输入格式来了
其他地方代码和报告的实现细节是完全一致的
*/

//代码使用模拟退火算法和遗传算法来解决置换流水车间调度问题

//工件结构体, 存储某工件在不同机器中的加工时间
struct Job { 
   vector<int> processingTimes;
};

//对指定序号的测试用例执行优化
void startOpt(int testCaseNumber);

//通用工具类函数
void   generateRandomSchedule(vector<int>& solution);// 生成随机调度顺序
int    calculateMakespan(const vector<Job>& jobs, const vector<int>& schedule, int m);// 计算特定调度顺序花费的时间
random_device rd;                      // 随机数种子
mt19937 gen(rd());                    //  生成随机数引擎
uniform_real_distribution<> dis(0,1);//   0到1间实数均匀分布, 使用dis(gen)可得到0到1之间的随机数

//模拟退火相关函数
vector<int> simulatedAnnealing(int testCaseNumber, vector<Job>& jobs, int m, double temp,  double endTemp, double coolingRate, int maxIterations);// 原始模拟退火算法主程序
vector<int> improvedSimulatedAnnealing(int testCaseNumber, vector<Job>& jobs, int m, double temp, double endTemp, double coolingRate, int maxIterations,int numSchedules);// 改进的模拟退火算法主程序

//遗传算法相关函数
vector<vector<int>>            initializePopulation(int populationSize, int numJobs);// 初始化种群
vector<pair<int, vector<int>>> evaluatePopulation(const vector<Job>& jobs, const vector<vector<int>>& population, int m);// 评估种群
vector<vector<int>>            rouletteWheelSelection(const vector<pair<int, vector<int>>>& evaluatedPopulation, int numToSelect);// 使用轮盘赌选择策略进行选择
vector<int>                    crossover(const vector<int>& parent1, const vector<int>& parent2);// 交叉操作
void                           mutate(vector<int>& individual, double mutationRate);// 突变操作
vector<int>                    geneticAlgorithm(int testCaseNumber, const vector<Job>& jobs, int m, int populationSize, int generations, double mutationRate, double eliteRate);// 遗传算法主程序

//输入输出函数
vector<Job> readJobsFromFile(const string& filename, int instanceNumber, int& m);// 从测试用例文件中读取工件信息
void        printJobsInfo(const vector<Job>& jobs);// 打印工件信息
void        printResult(const vector<int>& schedule, const vector<Job>& jobs, int m, int makespan);// 打印优化结果
void        printSchedule(const vector<int>& schedule);// 打印调度顺序

//调试函数
string generateDate();// 生成当前时间的字符串
void   writeDataToCSV(const string& filename, const vector<pair<int, int>>& data);// 将迭代数据写入CSV文件

//决定程序运行方式的全局变量
bool logOpen;          // 是否记录迭代数据（输出到csv文件）
bool repeatExperiment;//  是否重复实验
bool useSA;          //   是否使用模拟退火算法
bool useGA;         //    是否使用遗传算法
bool useiSA;       //     是否使用改进的模拟退火算法

int main() {
    srand(time(NULL));            // 初始化随机种子
    logOpen = false;             //  是否记录迭代数据到csv文件，默认不记录
    repeatExperiment = true;    //   是否重复实验，默认重复
    useSA = true;              //    是否使用模拟退火算法，默认使用
    useGA = true;             //     是否使用遗传算法，默认使用
    useiSA = false;          //      是否使用改进的模拟退火算法，默认不使用
    int testCaseNumber = -1;//       测试用例序号
    //输入测试用例的序号(输入0到10之间的整数测试单个用例,或者输入-1对全部用例进行测试
    // cout << "Please enter the test instance number (input an integer between 0 and 10 to test a single case, or input -1 to test all cases)";
    // cin >> testCaseNumber;
    if(testCaseNumber == -1){
        for(int i = 0; i < 11; i++){// 测试所有用例
            startOpt(i);
        }
    }else{
        startOpt(testCaseNumber);   // 测试指定用例
    }
    return 0;
}

// 对指定序号的测试用例执行优化
void startOpt(int testCaseNumber){
    cout << "----------------------------------" << endl;
    cout << "instance " << testCaseNumber << ":" << endl;
    string filePath = "testcase.txt";// 指定测试用例的文件路径
    int m;// 机器数
    vector<Job> jobs = readJobsFromFile(filePath, testCaseNumber, m);
    //printJobsInfo(jobs);// 打印文件的输入（用来测试一下输入是否正确）
    if(useSA){
        auto start = chrono::high_resolution_clock::now();// 计时开始(用于测试程序运行时间)
        //模拟退火算法
        cout << "  SA: " << endl;
        double temp = 10000;          // 初始温度
        double endTemp = 1e-30;      //  最低温度
        double coolingRate = 0.99;  //   降温百分比
        int maxIterations = 100000;//    最大迭代次数
        cout << "\tinitial temp: "<< temp << ", endTemp: " << endTemp << ", coolingRate: " << coolingRate << ", maxIterations: " << maxIterations << endl;
        if(repeatExperiment){    // 重复实验
            const int numExperiments = 10;// 重复实验的次数
            vector<int> bestSchedule;    //  用于存放最优排列
            vector<int> makespans;      //   用于存放每次实验的makespan
            int minMakespan = numeric_limits<int>::max(); // 初始化最小makespan为最大值
            int maxMakespan = numeric_limits<int>::min();//  初始化最大makespan为最小值
            int totalMakespan = 0;                      //   初始化总makespan为0
            for (int i = 0; i < numExperiments; ++i) {
                vector<int> schedule = simulatedAnnealing(testCaseNumber, jobs, m, temp, endTemp, coolingRate, maxIterations);// 调用模拟退火算法
                int currentMakespan = calculateMakespan(jobs, schedule, m);// 计算当前排列的makespan
                makespans.push_back(currentMakespan);
                if (currentMakespan < minMakespan) {
                    minMakespan = currentMakespan;// 更新最小makespan
                    bestSchedule = schedule;     //  更新最优排列
                }
                if (currentMakespan > maxMakespan) {
                    maxMakespan = currentMakespan;// 更新最大makespan
                }
                totalMakespan += currentMakespan;// 更新总makespan
            }
            double avgMakespan = static_cast<double>(totalMakespan) / numExperiments;// 计算平均makespan
            // 输出makespan结果
            cout << "\trepetition times of the experiment: " << numExperiments << ", makespan(min/avg/max):  " << minMakespan << "/" << avgMakespan  << "/" << maxMakespan << endl;
            // 输出最短makespan
            cout << "\tthe min makespan: " << minMakespan << endl;
            // 输出最优makespan对应的排列顺序
            cout << "\tthe optimal ";
            printSchedule(bestSchedule);
        }else{// 不重复实验
            vector<int> bestSchedule = simulatedAnnealing(testCaseNumber, jobs, m, temp, endTemp, coolingRate, maxIterations);// 调用模拟退火算法
            printResult(bestSchedule, jobs, m, calculateMakespan(jobs, bestSchedule, m));
        }
        auto end = chrono::high_resolution_clock::now();// 计时结束
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "\talgorithm runtime: " << duration << " milliseconds" << endl;// 输出算法运行时间
    }
    if(useGA){
        auto start = chrono::high_resolution_clock::now();// 计时开始(用于测试程序运行时间)
        //遗传算法
        cout << "  GA: " << endl;
        int populationSize = 500;   // 种群大小
        int generations = 500;     //  迭代次数
        double mutationRate = 0.5;//   突变率
        double eliteRate = 0.4;  //    精英留存率
        cout << "\tpopulationSize: "<< populationSize << ", generations: " << generations << ", mutation rate: " << mutationRate << ", eliteRate: " << eliteRate << endl;
        if(repeatExperiment){  // 重复实验
            const int numExperiments = 10;// 重复实验的次数
            vector<int> bestSchedule;    //  用于存放最优排列
            vector<int> makespans;      //   用于存放每次实验的makespan
            int minMakespan = numeric_limits<int>::max(); // 初始化最小makespan为最大值
            int maxMakespan = numeric_limits<int>::min();//  初始化最大makespan为最小值
            int totalMakespan = 0;
            for (int i = 0; i < numExperiments; ++i) {
                vector<int> schedule = geneticAlgorithm(testCaseNumber, jobs, m, populationSize, generations, mutationRate, eliteRate);// 调用遗传算法
                int currentMakespan = calculateMakespan(jobs, schedule, m);// 计算当前排列的makespan
                makespans.push_back(currentMakespan);
                if (currentMakespan < minMakespan) {
                    minMakespan = currentMakespan;// 更新最小makespan
                    bestSchedule = schedule;     //  更新最优排列
                }
                if (currentMakespan > maxMakespan) {
                    maxMakespan = currentMakespan;// 更新最大makespan
                }
                totalMakespan += currentMakespan;
            }
            double avgMakespan = static_cast<double>(totalMakespan) / numExperiments;// 计算平均makespan
            // 输出makespan结果
            cout << "\trepetition times of the experiment: " << numExperiments << ", makespan(min/avg/max):  " << minMakespan << "/" << avgMakespan  << "/" << maxMakespan << endl;
            // 输出最短makespan
            cout << "\tthe min makespan: " << minMakespan << endl;
            // 输出最优makespan对应的排列顺序
            cout << "\tthe optimal ";
            printSchedule(bestSchedule);
        }else{
            vector<int> bestSchedule = geneticAlgorithm(testCaseNumber, jobs, m, populationSize, generations, mutationRate, eliteRate);// 调用遗传算法
            printResult(bestSchedule, jobs, m, calculateMakespan(jobs, bestSchedule, m));
        }
        auto end = chrono::high_resolution_clock::now();// 计时结束
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "\talgorithm runtime: " << duration << " milliseconds" << endl;// 输出算法运行时间
    }
    if(useiSA){//改进的模拟退火算法
        auto start = chrono::high_resolution_clock::now();// 计时开始(用于测试程序运行时间)
        cout << "  iSA: " << endl;
        double temp = 10000;          // 初始温度
        double endTemp = 1e-30;      //  最低温度
        double coolingRate = 0.99;  //   降温百分比
        int maxIterations = 100000;//    最大迭代次数
        int numSchedules = 10;    //     改进的模拟退火算法中的同时迭代的解的个数
        cout << "\tinitial temp: "<< temp << ", endTemp: " << endTemp << ", coolingRate: " << coolingRate << ", maxIterations: " << maxIterations << ", numSchedules: " << numSchedules << endl;
        if(repeatExperiment){    // 重复实验
            const int numExperiments = 10;// 重复实验的次数
            vector<int> bestSchedule;    //  用于存放最优排列
            vector<int> makespans;      //   用于存放每次实验的makespan
            int minMakespan = numeric_limits<int>::max(); // 初始化最小makespan为最大值
            int maxMakespan = numeric_limits<int>::min();//  初始化最大makespan为最小值
            int totalMakespan = 0;                      //   初始化总makespan为0
            for (int i = 0; i < numExperiments; ++i) {
                vector<int> schedule = improvedSimulatedAnnealing(testCaseNumber, jobs, m, temp, endTemp, coolingRate, maxIterations, numSchedules);// 调用改进的模拟退火算法
                int currentMakespan = calculateMakespan(jobs, schedule, m);// 计算当前排列的makespan
                makespans.push_back(currentMakespan);
                if (currentMakespan < minMakespan) {
                    minMakespan = currentMakespan;// 更新最小makespan
                    bestSchedule = schedule;     //  更新最优排列
                }
                if (currentMakespan > maxMakespan) {
                    maxMakespan = currentMakespan;// 更新最大makespan
                }
                totalMakespan += currentMakespan;// 更新总makespan
            }
            double avgMakespan = static_cast<double>(totalMakespan) / numExperiments;// 计算平均makespan
            // 输出makespan结果
            cout << "\trepetition times of the experiment: " << numExperiments << ", makespan(min/avg/max):  " << minMakespan << "/" << avgMakespan  << "/" << maxMakespan << endl;
            // 输出最短makespan
            cout << "\tthe min makespan: " << minMakespan << endl;
            // 输出最优makespan对应的排列顺序
            cout << "\tthe optimal ";
            printSchedule(bestSchedule);
        }else{// 不重复实验
            vector<int> bestSchedule = improvedSimulatedAnnealing(testCaseNumber, jobs, m, temp, endTemp, coolingRate, maxIterations, numSchedules);// 调用改进的模拟退火算法
            printResult(bestSchedule, jobs, m, calculateMakespan(jobs, bestSchedule, m));
        }
        auto end = chrono::high_resolution_clock::now();// 计时结束
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "\talgorithm runtime: " << duration << " milliseconds" << endl;// 输出算法运行时间
    }
    cout << endl;
}

// 生成随机调度顺序
void generateRandomSchedule(vector<int>& schedule) {
    // 先将schedule初始化为 {0, 1, 2, ..., n-1}
    for (int i = 0; i < schedule.size(); i++) {
        schedule[i] = i;
    }
    // 将schedule{0, 1, 2, ..., n-1}内的顺序进行随机打乱
    shuffle(schedule.begin(), schedule.end(), mt19937(random_device()()));
}

// 计算特定调度顺序花费的时间
int calculateMakespan(const vector<Job>& jobs, const vector<int>& schedule, int m) {
    int n = schedule.size();
    vector<vector<int>> completion_time(n+1, vector<int>(m+1, 0));
    // 为什么是m+1和n+1呢？这样是为了方便对边界进行计算, 第一行和第一列都是0, 而计算是从第二行第二列开始进行的, 这样就不用特殊处理了
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            completion_time[i][j] = jobs[schedule[i-1]].processingTimes[j-1] +
                                    max(completion_time[i-1][j], completion_time[i][j-1]);
        }
    }
    return completion_time[n][m];
}

// 模拟退火算法主程序
vector<int> simulatedAnnealing(int testCaseNumber, vector<Job>& jobs, int m, double temp, double endTemp, double coolingRate, int maxIterations) {
    // temp是初始温度, coolingRate是降温百分比, maxIterations是最大迭代次数, endTemp是结束温度
    vector<pair<int, int>> iterationData;// 用于调试, 记录每次迭代的makespan和迭代次数
    string filename = "log_"+to_string(testCaseNumber) + "_SA_("+to_string(temp)+","+to_string(endTemp)+","+to_string(coolingRate)+","+to_string(maxIterations)+")_" + generateDate()+".csv";
    vector<int> schedule(jobs.size());  // 用来存放工件的调度顺序
    generateRandomSchedule(schedule);// 生成初始解
    vector<int> currentSchedule = schedule;
    int minMakespan = calculateMakespan(jobs, currentSchedule, m);
    for (int i = 0; i < maxIterations; i++) {
        // 随机选择两个位置交换
        int pos1 = rand() % schedule.size();
        int pos2 = rand() % schedule.size();
        swap(currentSchedule[pos1], currentSchedule[pos2]);
        // 计算新解的makespan
        int currentMakespan = calculateMakespan(jobs, currentSchedule, m);
        // 计算两次解的差值（模拟退火的能量变化）
        int energyChange = currentMakespan - minMakespan;
        // 如果新解更好, 或者根据概率接受更差的解
        if (energyChange < 0 || exp(-energyChange / temp) > (double)rand() / RAND_MAX) {
            minMakespan = currentMakespan;
            schedule = currentSchedule;  // 更新最优解
            if(logOpen){
                iterationData.push_back(make_pair(i + 1, minMakespan));
            }
        } else {
            swap(currentSchedule[pos1], currentSchedule[pos2]);  // 恢复原状
        }
        temp *= coolingRate;  // 降低温度
        if(temp < endTemp){  //  温度过低时提前结束
            break;
        }
    }
    if(logOpen){
        // 写入数据到csv文件
        writeDataToCSV(filename, iterationData);
    }
    return schedule;
}

// 初始化种群
vector<vector<int>> initializePopulation(int populationSize, int numJobs) {
    vector<vector<int>> population;            // 初始化种群
    for (int i = 0; i < populationSize; ++i) {//  populationSize对应种群内的个体数量
        vector<int> individual(numJobs);     //   numJobs对应每个个体内工件的数量
        generateRandomSchedule(individual); //    生成随机调度顺序
        population.push_back(individual);  //     将个体添加到种群中
    }
    return population;
}

// 评估种群
vector<pair<int, vector<int>>> evaluatePopulation(const vector<Job>& jobs, const vector<vector<int>>& population, int m) {
    vector<pair<int, vector<int>>> evaluated;
    for (const auto& individual : population) {
        int makespan = calculateMakespan(jobs, individual, m);// 每个个体的完工时间
        evaluated.push_back({makespan, individual});
    }
    // 按照完工时间从小到大排序
    sort(evaluated.begin(), evaluated.end(), [](const pair<int, vector<int>>& a, const pair<int, vector<int>>& b) {
        return a.first < b.first;
    });
    return evaluated;// 返回排序后的种群
}

// 轮盘赌选择策略
vector<vector<int>> rouletteWheelSelection(const vector<pair<int, vector<int>>>& evaluatedPopulation, int numToSelect) {
    vector<vector<int>> selected;
    double totalFitness = 0.0;
    vector<double> probabilities;
    // 计算总适应度
    for (const auto& p : evaluatedPopulation) {
        // 设适应度是完工时间的倒数, 完工时间越小, 适应度越高
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
                break; // 如果找到第一个大于等于randProb的概率, 停止搜索
            }
        }
        // 添加选中的个体到selected数组
        selected.push_back(evaluatedPopulation[idx].second);
    }
    return selected;
}

// 交叉操作
vector<int> crossover(const vector<int>& parent1, const vector<int>& parent2) {
    vector<int> child(parent1.size());
    int crossPoint = rand() % (parent1.size() - 2) + 1; // 保证交叉点不在末端
    copy_n(parent1.begin(), crossPoint, child.begin());// 从parent1复制前crossPoint个元素到child
    copy_if(parent2.begin(), parent2.end(), child.begin() + crossPoint, [&child](int gene) {// 从parent2复制剩余元素到child
        return find(child.begin(), child.end(), gene) == child.end();                      // 保证child中没有重复元素
    });
    return child;
}

// 突变操作
void mutate(vector<int>& individual, double mutationRate) {
    double chance = dis(gen);
    if (chance < mutationRate) {// 根据mutationRate的概率进行突变
        int index1 = rand() % individual.size();
        int index2 = rand() % individual.size();
        swap(individual[index1], individual[index2]);// 交换两个位置的元素
    }
}

// 遗传算法主程序
vector<int> geneticAlgorithm(int testCaseNumber, const vector<Job>& jobs, int m, int populationSize, int generations, double mutationRate, double eliteRate) {
    // populationSize是种群大小, generations是迭代次数, mutationRate是突变率, eliteRate是精英留存率
    vector<pair<int, int>> iterationData; // 用于调试, 记录每代的最佳makespan和迭代次数
    string filename = "log_"+to_string(testCaseNumber) + "_GA_("+to_string(populationSize)+","+to_string(generations)+","+to_string(mutationRate)+","+to_string(eliteRate)+")_" + generateDate()+".csv";
    // 初始化种群
    vector<vector<int>> population = initializePopulation(populationSize, jobs.size());
    int eliteSize = static_cast<int>(populationSize * eliteRate); // 保留40%的精英
    for (int gen = 0; gen < generations; ++gen) {
        auto evaluated = evaluatePopulation(jobs, population, m); // 评估种群,返回结果按照makespan从小到大排序
        vector<vector<int>> newGeneration; // 新一代种群
        // 保留精英个体
        for (int i = 0; i < eliteSize; ++i) {
            newGeneration.push_back(evaluated[i].second);
        }
        // 选择父代
        vector<vector<int>> parents = rouletteWheelSelection(evaluated, populationSize / 2); // 选择父代
        // 填充其余的新一代种群
        while (newGeneration.size() < populationSize) {
            int idx1 = rand() % parents.size();
            int idx2 = rand() % parents.size();                           // 随机选择两个父代
            vector<int> child = crossover(parents[idx1], parents[idx2]); // 交叉操作
            mutate(child, mutationRate);                                // 突变操作
            newGeneration.push_back(child);                            // 添加子代到新一代种群
        }
        population = newGeneration; // 更新种群
        if(logOpen){
            // 记录当前代中最佳makespan
            int bestMakespan = evaluatePopulation(jobs, population, m).front().first;
            iterationData.push_back(make_pair(gen + 1, bestMakespan));
        }
    }
    if(logOpen){
        // 写入数据到csv文件
        writeDataToCSV(filename, iterationData);
    }
    // 最终评估当前种群, 返回适应度最好的个体
    auto finalEvaluated = evaluatePopulation(jobs, population, m);
    return finalEvaluated.front().second;
}

// 从文件中读取工件信息
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
                startReading = true;   // 开始读取数据
                getline(file, line);  // 读取下一行, 即包含工件数和机器数的行
                istringstream iss(line);
                iss >> n >> m;
                continue;  // 跳过当前循环, 进入数据读取部分
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
        cerr << "file error" << endl;
    }
    return jobs;
}

// 打印工件信息
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

// 打印优化结果
void printResult(const vector<int>& schedule, const vector<Job>& jobs, int m, int makespan) {
    cout << "\tscheduling sequence is: ";
    for (size_t i = 0; i < schedule.size()-1; i++) {
        cout << schedule[i]+1 << ", ";
    }
    cout << schedule[schedule.size()-1]+1 << "";
    cout << endl;
    cout << "\tmakespan: " << makespan << endl;
}

// 打印调度顺序
void printSchedule(const vector<int>& schedule) {
    cout << "scheduling sequence is: ";
    for (size_t i = 0; i < schedule.size()-1; i++) {
        cout << schedule[i]+1 << ", ";
    }
    cout << schedule[schedule.size()-1]+1 << "";
    cout << endl;
}

// 生成当前时间的字符串
string generateDate() {
    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);
    stringstream ss;
    ss << put_time(localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}

// 用于将迭代数据写入CSV文件
void writeDataToCSV(const string& filename, const vector<pair<int, int>>& data) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    // 写入头部
    file << "Iteration,Makespan\n";
    // 写入数据
    for (const auto& entry : data) {
        file << entry.first << "," << entry.second << "\n";
    }
    file.close();
}

// 改进的模拟退火算法主程序
vector<int> improvedSimulatedAnnealing(int testCaseNumber, vector<Job>& jobs, int m, double temp, double endTemp, double coolingRate, int maxIterations,int numSchedules) {
    // temp是初始温度, coolingRate是降温百分比, maxIterations是最大迭代次数, endTemp是结束温度
    vector<vector<int>> schedules(numSchedules);
    for(int i = 0; i < numSchedules; i++){
        schedules[i].resize(jobs.size());
        generateRandomSchedule(schedules[i]);// 生成初始解
    }
    vector<vector<int>> currentSchedules = schedules;
    vector<int> minMakespans(numSchedules);
    for(int i = 0; i < numSchedules; i++){
        minMakespans[i] = calculateMakespan(jobs, currentSchedules[i], m);
    }
    for (int i = 0; i < maxIterations; i++) {
        for(int j = 0; j < numSchedules; j++){
            // 随机选择两个位置交换
            int pos1 = rand() % schedules[0].size();
            int pos2 = rand() % schedules[0].size();
            swap(currentSchedules[j][pos1], currentSchedules[j][pos2]);
            // 计算新解的makespan
            int currentMakespan = calculateMakespan(jobs, currentSchedules[j], m);
            // 计算两次解的差值（模拟退火的能量变化）
            int energyChange = currentMakespan - minMakespans[j];
            // 如果新解更好, 或者根据概率接受更差的解
            if (energyChange < 0 || exp(-energyChange / temp) > (double)rand() / RAND_MAX) {
                minMakespans[j] = currentMakespan;
            } else {
                swap(currentSchedules[j][pos1], currentSchedules[j][pos2]);  // 恢复原状
            }
        }
        temp *= coolingRate;  // 降低温度
        if(temp < endTemp){  //  温度过低时提前结束
            break;
        }
    }

    int minMakespan = minMakespans[0];
    int minIndex = 0;
    for(int i = 1; i < numSchedules; i++){
        if(minMakespans[i] < minMakespan){
            minMakespan = minMakespans[i];
            minIndex = i;
        }
    }
    return currentSchedules[minIndex];
}