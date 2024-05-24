#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

struct Job {
    vector<int> processingTimes;
};

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

vector<Job> readJobsFromFile(const string& filename, int& m) {
    ifstream file(filename);
    string line;
    vector<Job> jobs;
    
    if (file.is_open()) {
        getline(file, line);
        istringstream iss(line);
        int n;
        iss >> n >> m;
        
        while (getline(file, line)) {
            Job job;
            job.processingTimes.resize(m);
            istringstream lineStream(line);
            int machineIndex, time;
            
            for (int i = 0; i < m; ++i) {
                lineStream >> machineIndex >> time;
                job.processingTimes[machineIndex] = time;
            }
            jobs.push_back(job);
        }
        file.close();
    } else {
        cerr << "Could not open the file." << endl;
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

void simulatedAnnealing(vector<Job>& jobs, vector<int>& schedule, int m, int& minMakespan) {
    srand(time(NULL));  // 初始化随机种子
    double temp = 10000;  // 初始温度
    double coolingRate = 0.99;  // 冷却率
    int maxIterations = 1000;  // 最大迭代次数

    vector<int> currentSchedule = schedule;
    minMakespan = calculateMakespan(jobs, currentSchedule, m);

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
}
void printResult(const vector<int>& schedule, const vector<Job>& jobs, int m, int makespan) {
    std::cout << "调度顺序是：";
    for (size_t i = 0; i < schedule.size(); i++) {
        std::cout << schedule[i]+1 << " ";
    }
    std::cout << std::endl;
    cout << "Minimum makespan is: " << makespan << endl;

}

int main() {
    int testCaseNumber = 0;
    cout << "请输入测试用例的序号(从0到10之间选一个): ";
    cin >> testCaseNumber;
    string filePath = "testcase/" + to_string(testCaseNumber) + ".txt"; // 指定测试用例的文件路径
    int m;
    vector<Job> jobs = readJobsFromFile(filePath, m);
    printJobsInfo(jobs);// 打印工件信息（用来测试一下输入是否正确）

    vector<int> schedule(jobs.size());  // 用来存放工件的调度顺序
    for (int i = 0; i < schedule.size(); i++) {// 初始化调度顺序为顺序索引
        schedule[i] = i;
    }
    //schedule = {7, 2, 4, 3, 10, 1, 6, 9, 5, 0, 8};
    schedule = {7,4,8,2,3,0,1,0,1,9,6,5};
    //下面开始使用模拟退火算法进行优化
    int minMakespan;
    //simulatedAnnealing(jobs, schedule, m, minMakespan);
    //模拟退火算法优化结束
    minMakespan = calculateMakespan(jobs, schedule, m);
    printResult(schedule, jobs, m, minMakespan);
    

    return 0;
}

