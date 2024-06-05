import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# 存储每个工件的被每台机器加工的时间
class Job:
    def __init__(self, processing_times):
        self.processing_times = processing_times

# 计算调度顺序的完成时间
def calculate_makespan(jobs, schedule, m):
    n = len(schedule)
    completion_time = np.zeros((n+1, m+1))
    for i in range(1, n+1):
        for j in range(1, m+1):
            completion_time[i][j] = jobs[schedule[i-1]].processing_times[j-1] + \
                                    max(completion_time[i-1][j], completion_time[i][j-1])
    return int(completion_time[n][m]), completion_time

# 从文件中读取工件信息
def read_jobs_from_file(filename):
    jobs = []
    m = 0
    try:
        with open(filename, 'r') as file:
            first_line = next(file).strip()
            n, m = map(int, first_line.split())
            for line in file:
                parts = list(map(int, line.strip().split()))
                processing_times = [0]*m
                for i in range(m):
                    machine_index = parts[2*i]
                    time = parts[2*i+1]
                    processing_times[machine_index] = time
                jobs.append(Job(processing_times))
    except Exception as e:
        print(f"Could not open or read the file: {e}")
    return jobs, m

# 绘制甘特图
def plot_gantt_chart(jobs, schedule, m, completion_time ,test_case_number, makespan):
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = plt.cm.viridis(np.linspace(0, 0.8, m))
    # colors = plt.cm.cividis(np.linspace(0, 0.8, m))
    # 将高度设为1并且确保Y坐标使矩形紧密排列
    for i in range(1, len(schedule) + 1):
        for j in range(1, m + 1):
            start_time = completion_time[i][j] - jobs[schedule[i-1]].processing_times[j-1]
            duration = jobs[schedule[i-1]].processing_times[j-1]
            ax.add_patch(patches.Rectangle((start_time, m - j), duration, 1, edgecolor='black', facecolor=colors[j-1]))
            ax.text(start_time + duration / 2, m - j + 0.5, f'{schedule[i-1]+1}', color='white', weight='bold',
                    va='center', ha='center', fontsize=12)

    ax.set_xlim(0, completion_time[len(schedule)][m])
    ax.set_ylim(0, m)
    ax.set_yticks([m - 0.5 - i for i in range(m)])
    ax.set_yticklabels([f'M{i}' for i in range(m)])
    ax.set_xlabel('Time')
    ax.set_ylabel('Machine')
    ax.set_title('Gantt Chart for instance ' + str(test_case_number) +' (makespan:' + str(makespan)+')')
    plt.gca().invert_yaxis()
    plt.show()

def main():
    test_case_number = int(input("请输入测试用例的序号(从0到10之间选一个): "))
    file_path = f"testcase/{test_case_number}.txt"
    jobs, m = read_jobs_from_file(file_path)
    schedule_input = input("请输入工件的调度顺序，使用逗号分隔（例如: 8,1,3,11,2,5,7,9,4,6,10）: ")
    schedule = [int(num) - 1 for num in schedule_input.split(',') if num.strip().isdigit()]  # 将输入字符串分割，转换为整数，并减1
    makespan, completion_time = calculate_makespan(jobs, schedule, m)
    plot_gantt_chart(jobs, schedule, m, completion_time, test_case_number, makespan)
    
if __name__ == "__main__":
    main()
