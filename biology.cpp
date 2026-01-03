#define _USE_MATH_DEFINES
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <map>
#include <iomanip>

using namespace std;

// 生物个体相关
struct Organism {
    int id;
    double x, y;
    bool marked;
    int age;
    double moveSpeed;
    
    // 初始化个体函数实现
    Organism(int _id = 0, double _x = 0, double _y = 0, bool _marked = false, 
             int _age = 0, double _moveSpeed = 0.0) 
        : id(_id), x(_x), y(_y), marked(_marked), age(_age), moveSpeed(_moveSpeed) {}
};

// 种群相关
class Population {
private:
    int nextId;     // 下一个个体ID
    double area;    // 区域面积
    int Time;       // 时间
    
public:
    vector<Organism> organisms;    // 所有个体
    int initialCount = 0;
    
    // 初始化种群
    Population(int initial = -1) {
        nextId = 1;
        area = 1.0;
        Time = 0;
        
        // 初始化个体数量
        if (initial > 0) {
            initialCount = initial;
        } else {
            initialCount = rand() % 202 + 900; // 900-1100
        }
        
        // 初始化具体个体
        for (int i = 0; i < initialCount; ++i) {
            Organism org;
            org.id = nextId++;
            org.x = (double)rand() / RAND_MAX * area;
            org.y = (double)rand() / RAND_MAX * area;
            org.marked = false;
            org.age = rand() % 731;
            org.moveSpeed = 0.01 + ((double)rand() / RAND_MAX) * 0.02; // 0.01-0.03
            
            organisms.push_back(org);
        }
    }
    
    // 获取当前真实种群数量
    int getCount() const {
        return organisms.size();
    }
    
    // 获取面积
    double getArea() const {
        return area;
    }
    
    // 获取真实种群密度
    double getTrueDensity() const {
        if (area <= 0) return 0.0;
        return getCount() / area;
    }
    
    // 获取种群年龄结构
    void getAgeStructure(int& young, int& adult, int& old) const {
        young = 0;
        adult = 0;
        old = 0;
        
        for (int i = 0; i < organisms.size(); ++i) {
            int age = organisms[i].age;
            if (age <= 180) young++;
            else if (age <= 540) adult++;
            else old++;
        }
    }
    
    // 检查点是否在样方内
    bool isPointInRectangle(double x, double y, double recX, double recY, 
                            double width, double height) const {
        return (x >= recX && x <= recX + width && y >= recY && y <= recY + height);
    }
    
    // 移动个体
    void moveOrganisms() {
        for (int i = 0; i < organisms.size(); ++i) {
            // 根据年龄调整移动速度
            double speedFactor = 1.0;
            if (organisms[i].age <= 180) speedFactor = 1.2;
            else if (organisms[i].age > 540) speedFactor = 0.8;
            
            // 随机移动方向
            double angle = ((double)rand() / RAND_MAX) * 2.0 * M_PI;
            double distance = organisms[i].moveSpeed * speedFactor;
            
            // 计算新位置
            double newX = organisms[i].x + distance * cos(angle);
            double newY = organisms[i].y + distance * sin(angle);
            
            // 边界检查
            if (newX < 0) newX = -newX;
            if (newX > area) newX = 2 * area - newX;
            if (newY < 0) newY = -newY;
            if (newY > area) newY = 2 * area - newY;
            
            // 更新位置
            organisms[i].x = newX;
            organisms[i].y = newY;
        }
    }
    
    // 重置所有标记
    void resetAllMarks() {
        for (int i = 0; i < organisms.size(); ++i) {
            organisms[i].marked = false;
        }
    }
    
    // 获取移动统计信息
    void getMovementStats(double& avgSpeed, double& minSpeed, double& maxSpeed) const {
        if (organisms.empty()) {
            avgSpeed = minSpeed = maxSpeed = 0;
            return;
        }
        
        double totalSpeed = 0;
        minSpeed = organisms[0].moveSpeed;
        maxSpeed = organisms[0].moveSpeed;
        
        for (int i = 0; i < organisms.size(); ++i) {
            double speed = organisms[i].moveSpeed;
            totalSpeed += speed;
            if (speed < minSpeed) minSpeed = speed;
            if (speed > maxSpeed) maxSpeed = speed;
        }
        
        avgSpeed = totalSpeed / organisms.size();
    }
    
    // 种群更新程序
    void updatePopulation() {
        Time++;
        int currentCount = getCount();
        if (currentCount == 0) return; // 种群灭绝
        
        // 1. 个体年龄增加
        for (int i = 0; i < organisms.size(); ++i) {
            organisms[i].age++;
        }
        
        // 2. 个体自然死亡
        vector<int> toRemove;
        
        // 死亡概率与判断
        for (int i = 0; i < organisms.size(); ++i) {
            // 基础死亡概率
            double deathRate = 0.002;
            
            // 年龄影响死亡率
            int age = organisms[i].age;
            if (age <= 180) deathRate = 0.003;
            else if (age > 540) deathRate = 0.005;
            
            // 随机扰动
            double randomFactor = ((double)rand() / RAND_MAX) * 0.001 - 0.0005;
            deathRate += randomFactor;
            if (deathRate < 0) deathRate = 0.0;
            
            // 个体死亡判定
            if (((double)rand() / RAND_MAX) < deathRate) {
                toRemove.push_back(i);
            }
        }
        
        // 移除死亡个体（从后往前，避免索引问题）
        for (int i = toRemove.size() - 1; i >= 0; --i) {
            int idx = toRemove[i];
            organisms.erase(organisms.begin() + idx);
        }
        
        // 3. 个体繁殖
        // 统计成年个体数
        int adultCount = 0;
        for (int i = 0; i < organisms.size(); ++i) {
            int age = organisms[i].age;
            if (age > 180 && age <= 540) adultCount++;
        }
        
        // 密度影响出生率
        double density = getTrueDensity();
        double birthRate = 0.005; // 基础出生率
        
        double densityFactor = 1.0;
        if (density > 1000.0) densityFactor = 1.0 - (density - 1000.0) / 5000.0;
        if (densityFactor < 0.1) densityFactor = 0.1;
        double actualBirthRate = birthRate * densityFactor;
        
        int expectedBirths = (int)(adultCount * actualBirthRate);
        
        // 随机波动
        int minBirths = (int)(expectedBirths * 0.7);
        if (minBirths < 0) minBirths = 0;
        int maxBirths = (int)(expectedBirths * 1.3);
        
        int births = 0;
        if (maxBirths > 0) {
            births = minBirths + rand() % (maxBirths - minBirths + 1);
        }
        
        // 4. 生成后代
        for (int i = 0; i < births; ++i) {
            if (organisms.empty()) break; // 无成年个体无法繁殖
            
            // 随机选择亲代
            int parentIdx = rand() % organisms.size();
            Organism parent = organisms[parentIdx];
            double parentX = parent.x;
            double parentY = parent.y;
            
            // 后代位置在亲代附近
            double offsetX = ((double)rand() / RAND_MAX - 0.5) * 0.1; // -0.05到0.05
            double offsetY = ((double)rand() / RAND_MAX - 0.5) * 0.1; // -0.05到0.05
            
            double childX = parentX + offsetX;
            double childY = parentY + offsetY;
            
            // 确保后代位置在区域内
            if (childX < 0) childX = 0;
            if (childX > area) childX = area;
            if (childY < 0) childY = 0;
            if (childY > area) childY = area;
            
            // 创建后代个体
            Organism child;
            child.id = nextId++;
            child.x = childX;
            child.y = childY;
            child.marked = false;
            child.age = 0;
            
            organisms.push_back(child);
        }
        
        // 5. 迁移过程
        // 迁出率与密度相关
        double emigrationRate = 0.003 * (density / 1000.0);
        int emigrants = (int)(currentCount * emigrationRate);
        if (emigrants > organisms.size()) emigrants = getCount();
        
        // 随机选择迁出个体
        for (int i = 0; i < emigrants; ++i) {
            if (organisms.empty()) break;
            int idx = rand() % organisms.size();
            organisms.erase(organisms.begin() + idx);
        }
        
        // 迁入个体
        double immigrationRate = 0.005;
        int immigrants = (int)(currentCount * immigrationRate);
        
        for (int i = 0; i < immigrants; ++i) {
            Organism immigrant;
            immigrant.id = nextId++;
            immigrant.x = (double)rand() / RAND_MAX;
            immigrant.y = (double)rand() / RAND_MAX;
            immigrant.marked = false;
            immigrant.age = 180 + rand() % 361; // 180-540, 成年个体
            immigrant.moveSpeed = 0.01 + ((double)rand() / RAND_MAX) * 0.02;
            
            organisms.push_back(immigrant);
        }
        
        // 6. 个体移动
        moveOrganisms();
        
        // 7. 环境承受限制
        int maxCapacity = (int)(area * 5000);
        if (getCount() > maxCapacity) {
            int toRemoveCount = getCount() - maxCapacity;
            for (int i = 0; i < toRemoveCount; ++i) {
                if (organisms.empty()) break;
                int idx = rand() % organisms.size();
                organisms.erase(organisms.begin() + idx);
            }
        }
    }
};

// 计数方法类
class CountMethod {
public:
    // 1. 逐个计数法
    static double totalCount(const Population& population) {
        if (population.getArea() <= 0) return 0.0;
        int count = population.getCount();
        return count / population.getArea();
    }
    
    // 2. 五点样方法
    static double fivePointQuadratMethod(const Population& population) {
        // 样方面积
        double sampleArea = 0.01;
        
        // 样方中心点坐标
        double centers[5][2] = {
            {0.5, 0.5},   // 中心
            {0.25, 0.25}, // 左下
            {0.75, 0.25}, // 右下
            {0.25, 0.75}, // 左上
            {0.75, 0.75}  // 右上
        };
        
        double side = sqrt(sampleArea);
        int counts[5] = {0, 0, 0, 0, 0};
        
        // 统计每个样方内的个体数
        for (int i = 0; i < population.organisms.size(); ++i) {
            double x = population.organisms[i].x;
            double y = population.organisms[i].y;
            
            for (int j = 0; j < 5; ++j) {
                double centerX = centers[j][0];
                double centerY = centers[j][1];
                
                // 计算样方边界
                double recX = centerX - side / 2.0;
                double recY = centerY - side / 2.0;
                if (recX < 0) recX = 0;
                if (recY < 0) recY = 0;
                if (recX + side > population.getArea()) recX = population.getArea() - side;
                if (recY + side > population.getArea()) recY = population.getArea() - side;
                
                // 检查个体是否在样方内
                if (population.isPointInRectangle(x, y, recX, recY, side, side)) {
                    counts[j]++;
                    break; // 一个个体只计入一个样方
                }
            }
        }
        
        // 计算平均密度
        double total = 0;
        for (int i = 0; i < 5; ++i) {
            total += counts[i];
        }
        double averageCount = total / 5.0;
        return averageCount / sampleArea;
    }
    
    // 3. 等距样方法
    static double equidistanceQuadratMethod(const Population& population) {
        // 基本参数
        double sampleArea = 0.01;
        const int numSamples = 10;
        double side = sqrt(sampleArea);
        
        // 等距布置样方
        double spacing = (population.getArea() - side) / (numSamples - 1); // 避免超出边界
        int counts[numSamples] = {0};
        
        // 统计样方内个体数
        for (int i = 0; i < population.organisms.size(); i++) {
            double x = population.organisms[i].x;
            double y = population.organisms[i].y;
            
            for (int j = 0; j < numSamples; ++j) {
                // 样方中心
                double centerX = spacing * (j + 1);
                double centerY = 0.5; // y方向居中
                
                // 计算样方边界
                double recX = centerX - side / 2;
                double recY = centerY - side / 2;
                
                if (recX < 0) recX = 0;
                if (recY < 0) recY = 0;
                if (recX > 1.0) recX = 1.0;
                if (recY > 1.0) recY = 1.0;
                
                // 检查点是否在样方内
                if (population.isPointInRectangle(x, y, recX, recY, side, side)) {
                    counts[j]++;
                    break;
                }
            }
        }
        
        // 计算平均密度
        double total = 0;
        for (int i = 0; i < numSamples; i++) {
            total += counts[i];
        }
        double average = total / numSamples;
        return average / sampleArea;
    }
    
    // 4. 标记重捕法
    static double markRecaptureMethod(Population& population) {
        int captureNum = 75;
        int populationSize = population.getCount();
        
        // 检测种群大小
        if (populationSize < captureNum) {
            return -1; // 返回无效值
        }
        
        // 重置标记
        population.resetAllMarks();
        
        // 第一次捕捉
        int markedCount = 0;
        for (int i = 0; i < captureNum; i++) {
            if (population.organisms.size() == 0) break;
            int idx = rand() % population.organisms.size();
            population.organisms[idx].marked = true;
            markedCount++;
        }
        
        // 第二次捕捉
        int recaptureMarked = 0;
        int totalRecaptureCount = 0;
        
        for (int i = 0; i < captureNum; i++) {
            if (population.organisms.size() == 0) break;
            int idx = rand() % population.organisms.size();
            totalRecaptureCount++;
            
            if (population.organisms[idx].marked) {
                recaptureMarked++;
            }
        }
        
        // 如果重捕无标记
        if (recaptureMarked == 0) {
            population.resetAllMarks();
            return -1;
        }
        
        // 计算种群数
        double estimatedCount = (markedCount * captureNum / recaptureMarked);
        
        population.resetAllMarks();
        // 计算密度
        return estimatedCount / population.getArea();
    }
    
    // 执行所有调查方法并返回结果
    static void runAllMethods(const Population& originalPopulation, 
                            double& totalCountResult, double& fivePointResult, 
                            double& equidistanResult, double& markRecaptureResult) {
        // 逐个计数法
        totalCountResult = totalCount(originalPopulation);
        // 五点样方法
        fivePointResult = fivePointQuadratMethod(originalPopulation);
        // 等距样方法
        equidistanResult = equidistanceQuadratMethod(originalPopulation);
        // 标记重捕法
        Population populationCopy = originalPopulation; // 创建副本
        markRecaptureResult = markRecaptureMethod(populationCopy);
    }
    
};

// 实验模拟类
class PopulationExperiment {
private:
    Population population;
    int Time;
    
    // 打印单天结果
    void printDayResult(int day) {
        double trueDensity = population.getTrueDensity();
        int trueCount = population.getCount();
        
        cout << "│" << setw(5) << day << " │";
        cout << setw(11) << trueCount << " │";
        cout << fixed << setprecision(2) << setw(11) << trueDensity << " │";
        
        // 估计密度
        double totalCountResult, fivePointResult, equidistanResult, markRecaptureResult;
        CountMethod::runAllMethods(population, totalCountResult, fivePointResult, 
                                  equidistanResult, markRecaptureResult);
        
        cout << fixed << setprecision(2) << setw(11) << totalCountResult << " │";
        cout << setw(11) << fivePointResult << " │";
        cout << setw(11) << equidistanResult << " │";
        
        if (markRecaptureResult >= 0) {
            cout << setw(11) << markRecaptureResult << " │" << endl;
        } else {
            cout << setw(11) << "N/A" << " │" << endl;
        }
        
        // 每5天添加一个分隔线，增加可读性
        if ((day + 1) % 5 == 0 && day < Time) {
            cout << "├──────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤" << endl;
        }
    }
    
    // 打印最终结果比较
    void printFinalComparison() {
        double trueDensity = population.getTrueDensity();
        double totalCountResult, fivePointResult, equidistanResult, markRecaptureResult;
        
        CountMethod::runAllMethods(population, totalCountResult, fivePointResult, 
                                  equidistanResult, markRecaptureResult);
        
        cout << "\n最终结果比较:" << endl;
        cout << "  ┌─────────────────┬────────────┬──────────────┐" << endl;
        cout << "  │     调查方法    │  估计密度  │   误差(%)    │" << endl;
        cout << "  ├─────────────────┼────────────┼──────────────┤" << endl;
        
        cout << "  │     真实密度    │" << fixed << setprecision(2) << setw(11) << trueDensity << " │" << setw(13) << " - " << " │" << endl;
        
        cout << "  ├─────────────────┼────────────┼──────────────┤" << endl;
        cout << "  │    逐个计数法   │" << setw(11) << totalCountResult << " │" << setw(13) << "0.00" << " │" << endl;
        
        if (fivePointResult >= 0) {
            double error = fabs((fivePointResult - trueDensity) / trueDensity * 100);
            cout << "  │    五点样方法   │" << setw(11) << fivePointResult << " │" << fixed << setprecision(2) << setw(13) << error << " │" << endl;
        }
        
        if (equidistanResult >= 0) {
            double error = fabs((equidistanResult - trueDensity) / trueDensity * 100);
            cout << "  │    等距样方法   │" << setw(11) << equidistanResult << " │" << setw(13) << error << " │" << endl;
        }
        
        if (markRecaptureResult >= 0) {
            double error = fabs((markRecaptureResult - trueDensity) / trueDensity * 100);
            cout << "  │    标记重捕法   │" << setw(11) << markRecaptureResult << " │" << setw(13) << error << " │" << endl;
        }
        
        cout << "  └─────────────────┴────────────┴──────────────┘" << endl;
    }
    
public:
    PopulationExperiment(int initialCount = -1, int simulationDays = 10) 
        : population(initialCount) {
        Time = simulationDays;
    }
    
    void run() {
        // 打印实验标题
        cout << "\n";
        cout << "╔══════════════════════════════════════════════════════════════════════╗" << endl;
        cout << "║                  生物种群密度估计模拟实验                            ║" << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════╣" << endl;
        cout << "║  初始种群数量: " << setw(10) << population.initialCount << "                                            ║" << endl;
        cout << "║  模拟天数: " << setw(12) << Time << "                                              ║" << endl;
        cout << "║  区域面积: " << setw(12) << "\t1.0 平方米" << "                                             ║" << endl;
        cout << "╚══════════════════════════════════════════════════════════════════════╝" << endl;
        cout << "\n";
        
        // 打印种群信息标题
        cout << "当前种群信息:" << endl;
        int young = 0, adult = 0, old = 0;
        population.getAgeStructure(young, adult, old);
        double avgSpeed, minSpeed, maxSpeed;
        population.getMovementStats(avgSpeed, minSpeed, maxSpeed);
        
        cout << "  种群数量: " << population.getCount() << " [幼年: " << young << ", 成年: " << adult << ", 老年: " << old << "]" << endl;
        cout << "  真实密度: " << fixed << setprecision(2) << population.getTrueDensity() << " 个/平方米" << endl;
        cout << "  移动速度: 平均=" << fixed << setprecision(4) << avgSpeed 
             << ", 最小=" << minSpeed << ", 最大=" << maxSpeed << endl;
        cout << "\n";
        
        // 打印表格标题
        cout << "┌──────┬────────────┬────────────┬────────────┬────────────┬────────────┬────────────┐" << endl;
        cout << "│ 天数 │  真实数量  │  真实密度  │ 逐个计数法 │ 五点样方法 │ 等距样方法 │ 标记重捕法 │" << endl;
        cout << "├──────┼────────────┼────────────┼────────────┼────────────┼────────────┼────────────┤" << endl;
        
        // 逐天实验
        for (int day = 0; day <= Time; day++) {
            if (day == 0) {
                printDayResult(day);
            } else {
                population.updatePopulation();
                printDayResult(day);
            }
        }
        
        cout << "└──────┴────────────┴────────────┴────────────┴────────────┴────────────┴────────────┘" << endl;
        
        // 打印最终结果比较
        printFinalComparison();
        
        // 实验总结
        cout << "\n";
        cout << "╔══════════════════════════════════════════════════════════════════════╗" << endl;
        cout << "║                             实验总结                                 ║" << endl;
        cout << "╠══════════════════════════════════════════════════════════════════════╣" << endl;
        
        
        cout << "\n生态学参数设置:" << endl;
        cout << "  ┌──────────────────────────────────┬──────────────────────────────────┐" << endl;
        cout << "  │          参数类别                │            参数值                │" << endl;
        cout << "  ├──────────────────────────────────┼──────────────────────────────────┤" << endl;
        cout << "  │  基础出生率                      │  0.005 (密度制约)                │" << endl;
        cout << "  │  基础死亡率                      │  0.002                           │" << endl;
        cout << "  │  幼年死亡率                      │  0.003                           │" << endl;
        cout << "  │  老年死亡率                      │  0.005                           │" << endl;
        cout << "  │  基础迁出率                      │  0.0002 (密度相关)               │" << endl;
        cout << "  │  基础迁入率                      │  0.0005                          │" << endl;
        cout << "  │  环境承载力                      │  5000 个/平方米                  │" << endl;
        cout << "  │  年龄分界点                      │  幼年<180, 成年180-540, 老年>540 │" << endl;
        cout << "  │  个体移动速度                    │  0.01-0.03                       │" << endl;
        cout << "  │  初始数量范围                    │  900-1100                        │" << endl;
        cout << "  │  寿命范围                        │  0-730 天                        │" << endl;
        cout << "  │  捕捉数量                        │  75                              │" << endl;
        cout << "  │  模拟天数                        │  " << Time << "                              │" << endl;
        cout << "  └──────────────────────────────────┴──────────────────────────────────┘" << endl;
        cout << "╚══════════════════════════════════════════════════════════════════════╝" << endl;
    }
};

// 主函数
int main() {
    // 初始化随机种子
    srand(static_cast<unsigned int>(time(nullptr)));
    
    cout << fixed << setprecision(2);
    
    PopulationExperiment experiment(-1, 100);
    experiment.run();
    
    cout << "\n模拟实验完成，按任意键退出..." << endl;
    cin.get();
    
    return 0;
}