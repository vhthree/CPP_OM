import random
import math
from typing import List, Tuple

class Job:
    def __init__(self, processing_times: List[int]):
        self.processing_times = processing_times

def read_jobs_from_file(filename: str, instance_number: int) -> Tuple[List[Job], int]:
    jobs = []
    start_reading = False
    m = 0  # Number of machines
    with open(filename, 'r') as file:
        for line in file:
            if not start_reading and f"instance {instance_number}" in line:
                start_reading = True
                line = next(file)  # Read the line with job count and machine count
                m = int(line.split()[1])  # Extract the number of machines
                continue
            if start_reading:
                if line.strip() == "":
                    break
                processing_times = list(map(int, line.split()[1::2]))
                jobs.append(Job(processing_times))
    return jobs, m

def calculate_makespan(jobs: List[Job], schedule: List[int], m: int) -> int:
    n = len(schedule)
    completion_time = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            completion_time[i][j] = (
                jobs[schedule[i - 1]].processing_times[j - 1] +
                max(completion_time[i - 1][j], completion_time[i][j - 1])
            )
    return completion_time[n][m]

def generate_initial_solution(n: int) -> List[int]:
    solution = list(range(n))
    random.shuffle(solution)
    return solution

def simulated_annealing(jobs: List[Job], m: int, temp: float, cooling_rate: float, max_iterations: int) -> Tuple[List[int], int]:
    schedule = generate_initial_solution(len(jobs))
    min_makespan = calculate_makespan(jobs, schedule, m)
    for _ in range(max_iterations):
        pos1, pos2 = random.sample(range(len(schedule)), 2)
        # Swap two elements
        schedule[pos1], schedule[pos2] = schedule[pos2], schedule[pos1]
        new_makespan = calculate_makespan(jobs, schedule, m)
        if new_makespan < min_makespan or math.exp((min_makespan - new_makespan) / temp) > random.random():
            min_makespan = new_makespan
        else:
            # Swap back
            schedule[pos1], schedule[pos2] = schedule[pos2], schedule[pos1]
        temp *= cooling_rate
    return schedule, min_makespan

def start_optimize(filename: str, instance_number: int):
    jobs, m = read_jobs_from_file(filename, instance_number)
    final_schedule, min_makespan = simulated_annealing(jobs, m, 20000, 0.8, 3000)
    print("测试用例", instance_number, "的优化结果：")
    print("调度顺序是：", " ".join(str(s + 1) for s in final_schedule))
    print("Minimum makespan is:", min_makespan)

if __name__ == "__main__":
    filename = "testcase2023.txt"  # Change to your actual file name
    for i in range(11):
        start_optimize(filename, i)
