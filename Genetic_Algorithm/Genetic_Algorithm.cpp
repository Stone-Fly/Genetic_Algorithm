#include<stdio.h>
#include<string.h>
#include<string>
#include<algorithm>
#include<queue>
#include<stack>
#include<map>
#include<set>
#include<functional>
#include<vector>
#include<iostream>
#include<time.h>
#include<math.h>
using namespace std;
//-------------------------------------------------------
//#define LL long long
#define LL __int64
template<typename T>T Max(T a, T b){ return a>b ? a : b; }
template<typename T>T Min(T a, T b){ return a<b ? a : b; }
template<typename T>void Swap(T &a, T &b){ T c = a; a = b; b = c; }
template<typename T>T Abs(T a){ return a >= 0 ? a : -a; }
#define clr_all(a,b) memset(a,b,sizeof(a))
#define clr(a,b,n) memset(a,b,sizeof(a[0])*n)
const int MAXN = 2000 + 5;
const int MAXM = 2000 + 5;
const int mod = 1000000007;
const int INF = 0x7FFFFFFF;
//-------------------------------------------------------

typedef struct node{
	double x, y;
}Node;
typedef struct population{
	int Order[MAXN];
	double Length;
	double Adapt;
	double Probability;
}Pop;
Node City[MAXN];					//城市坐标
int City_top;						//城市总数
double dis[MAXN][MAXN];				//城市距离
Pop Population[MAXM];				//种群，供遗传算法使用
int Best_City_Order[MAXN];			//最优城市排列
double Min_Length = (double)INF;	//最短距离

void Input(){
	//输入城市坐标
	FILE *fp = fopen("../ch318.tsp", "r");
	if (fp == NULL){
		printf("请将数据集文件放入相应位置\n");
	}
	else{
		int id;
		double x, y;
		City_top = 0;
		while (fscanf(fp, "%d%lf%lf", &id, &x, &y) != EOF){
			City[id].x = x;
			City[id].y = y;
			City_top = id;
		}
		printf("City_top = %d\n", City_top);
		fclose(fp);
	}
}

double Cal_Dis(Node a, Node b){
	//计算两城市间的距离
	double x = a.x - b.x;
	double y = a.y - b.y;
	return sqrt(x*x + y*y);
}

double Cal_Length(int *Order){
	//计算城市排列总长度
	double Length = 0;
	for (int i = 1; i <= City_top; i++){
		Length += dis[Order[i]][Order[i + 1]];
	}
	return Length;
}

void Checking(){
	//利用最优答案测试距离计算是否正确
	//FILE *fp = fopen("../ch130.opt.tour.txt", "r");
	FILE *fp = fopen("../test.txt", "r");
	if (fp == NULL){
		printf("请将数据集文件放入相应位置\n");
	}
	else{
		bool flag[MAXN] = { 0 };
		int Order[MAXN];
		int Order_top = 0;
		while (fscanf(fp, "%d", &Order[++Order_top]) != EOF);
		for (int i = 1; i <= City_top; i++){
			if (!flag[Order[i]])
				flag[Order[i]] = true;
			else{
				printf("城市有重复，重复城市为：%d\n", Order[i]);
				break;
			}
		}
		Order[City_top + 1] = Order[1];
		double Length = Cal_Length(Order);
		printf("官方最优解的长度为：%.10f\n", Length);
		fclose(fp);
	}
}

void Init_City_Order(int *Order){
	//生成一个城市的排列
	//排列为1..ci..cj..cn..1
	//排列总数City_top + 1
	bool Flag[MAXN] = { 1, 1 };
	int Order_top = 0;
	Order[++Order_top] = 1;
	for (int i = 2; i <= City_top; i++){
		int tmp;
		while (tmp = rand() % City_top + 1){
			if (!Flag[tmp])
				break;
		}
		Order[++Order_top] = tmp;
		Flag[tmp] = true;
	}
	Order[Order_top + 1] = 1;
}

void Improved_circle(int *Order){
	//改良圈算法
	bool flag = true;
	while (flag){
		flag = false;
		for (int i = 1; i <= City_top - 2; i++){
			for (int j = i + 2; j <= City_top; j++){
				if (dis[Order[i]][Order[j]] + dis[Order[i + 1]][Order[j + 1]] < dis[Order[i]][Order[i + 1]] + dis[Order[j]][Order[j + 1]]){
					for (int l = i + 1, r = j; l < r; l++, r--)
						Swap(Order[l], Order[r]);
					flag = true;
				}
			}
		}
	}
}

void Print(int M){
	for (int i = 1; i <= M; i++){
		for (int j = 1; j <= City_top; j++){
			printf("%d%c", Population[i].Order[j], j % 5 == 0 ? '\n' : ' ');
		}
		printf("\n");
	}
}

void Init_Population(int M){
	//生成初始种群，合计 M 个
	for (int i = 1; i <= M; i++){
		Init_City_Order(Population[i].Order);
		Improved_circle(Population[i].Order);
	}
}

void Init(int M){
	//计算每个城市之间的距离
	for (int i = 1; i <= City_top; i++){
		for (int j = 1; j <= City_top; j++){
			dis[i][j] = Cal_Dis(City[i], City[j]);
		}
	}
	//Checking();
	Init_Population(M);
}

double Cal_Adapt(int M){
	double mMax = 0;
	for (int i = 1; i <= M; i++){
		//计算每组城市列表的适应度，值为 1 / Length
		Population[i].Length = Cal_Length(Population[i].Order);
		Population[i].Adapt = 1 / Population[i].Length;
		mMax = Max(mMax, Population[i].Adapt);
	}
	return mMax;
}

void Compare_Min_Length(int M){
	for (int i = 1; i <= M; i++){
		if (Population[i].Length < Min_Length){
			Min_Length = Population[i].Length;
			for (int j = 1; j <= City_top; j++)
				Best_City_Order[j] = Population[i].Order[j];
		}
	}
}

Pop tmp[MAXM];
void Select(int M){
	double sum = 0;
	for (int i = 1; i <= M; i++)
		sum += Population[i].Adapt;
	//计算概率
	for (int i = 1; i <= M; i++)
		Population[i].Probability = Population[i].Adapt / sum;
	//开启轮盘赌
	for (int i = 1; i <= M; i++){
		double pick = ((double)rand()) / RAND_MAX;
		for (int j = 1; j <= M; j++){
			pick -= Population[j].Probability;
			if (pick <= 0){
				tmp[i] = Population[j];
			}
		}
	}
	for (int i = 1; i <= M; i++)
		Population[i] = tmp[i];
}

double Pcross = 0.6;
void Cross(int M){
	int choice1, choice2;
	int pos1, pos2;
	int conflict1[MAXN];
	int conflict2[MAXN];
	for (int move = 1; move < M; move += 2){
		double pick = ((double)rand()) / RAND_MAX;
		if (pick>Pcross)
			continue;
		choice1 = move;
		choice2 = move + 1;
		pos1 = rand() % (City_top - 1) + 2;
		pos2 = rand() % (City_top - 1) + 2;
		if (pos1 > pos2)
			Swap(pos1, pos2);
		//交换pos1到pos2间的基因序列
		for (int i = pos1; i <= pos2; i++)
			Swap(Population[choice1].Order[i], Population[choice2].Order[i]);
		int num1 = 0, num2 = 0;
		for (int i = 2; i <= City_top; i++){
			if (i >= pos1&&i <= pos2)
				continue;
			for (int j = pos1; j <= pos2; j++){
				if (Population[choice1].Order[i] == Population[choice1].Order[j])
					conflict1[num1++] = i;
				if (Population[choice2].Order[i] == Population[choice2].Order[j])
					conflict2[num2++] = i;
			}
		}
		if ((num1 == num2) && num1 > 0){
			for (int i = 0; i < num1; i++)
				Swap(Population[choice1].Order[conflict1[i]], Population[choice2].Order[conflict2[i]]);
		}
	}
}

double Pvariation = 0.6;
void Variation(int M){
	for (int i = 1; i <= M; i++){
		double pick = ((double)rand()) / RAND_MAX;
		if (pick>Pvariation)
			continue;
		int pos1, pos2;
		pos1 = rand() % (City_top - 1) + 2;
		pos2 = rand() % (City_top - 1) + 2;
		Swap(Population[i].Order[pos1], Population[i].Order[pos2]);
	}
}

Pop tmp2[MAXM];
void Genetic(int T, int M, int N, int K){
	//T为进化的最高代数，M为种群总数
	FILE *fp = fopen("../City_Order_318.txt", "w");
	double Max_Adapt = 0;
	double MAX_ADAPT = 1 / 42030.0;
	Max_Adapt = Cal_Adapt(M);
	int n = 0, k = 0;
	for (int t = 0; t <= T&&Max_Adapt < MAX_ADAPT; t++){
		Cal_Adapt(M);
		Compare_Min_Length(M);
		for (int i = 1; i <= M; i++)
			tmp2[i] = Population[i];
		Select(M);
		n = 0;
		while (1){
			Cross(M);
			Variation(M);
			for (int i = 1; i <= M; i++)
				Improved_circle(Population[i].Order);
			double max_Adapt = Cal_Adapt(M);
			if (max_Adapt >= Max_Adapt){
				//录优
				Max_Adapt = max_Adapt;
				k = 0;
				break;
			}
			else{
				n++;
				if (n == N){
					k++;
					if (k == K){
						//重构
						Init_Population(M);
						break;
					}
					else{
						//重选
						/*for (int i = 1; i <= M; i++)
						Population[i] = tmp2[i];*/
						break;
					}
				}
				else{
					//复原
					for (int i = 1; i <= M; i++)
						Population[i] = tmp2[i];
				}
			}
		}
		printf("迭代%d的最优解长度为：%.10f\n", t, Min_Length);
		fprintf(fp, "--------------------------------------------------\n");
		fprintf(fp, "迭代%d的最优解长度为：%.10f\n", t, Min_Length);
		for (int i = 1; i <= City_top; i++)
			fprintf(fp, "%d%c", Best_City_Order[i], i % 10 == 0 ? '\n' : ' ');
	}
	fclose(fp);
}

void Output(int M){
	FILE *fp = fopen("../City_Order_1002.txt", "a");
	Cal_Adapt(M);
	Compare_Min_Length(M);
	printf("--------------------------------------------------\n");
	printf("最终最优解长度为：%.10f\n", Min_Length);
	printf("--------------------------------------------------\n对应城市排列为：\n");
	for (int i = 1; i <= City_top; i++)
		printf("%d%c", Best_City_Order[i], i % 10 == 0 ? '\n' : ' ');
	fprintf(fp, "--------------------------------------------------\n");
	fprintf(fp, "最终最优解长度为：%.10f\n", Min_Length);
	for (int i = 1; i <= City_top; i++)
		fprintf(fp, "%d%c", Best_City_Order[i], i % 10 == 0 ? '\n' : ' ');
	fclose(fp);
}

int main(){
	int T = 100000, M = 1000;
	srand((unsigned)time(NULL));
	clock_t begin, end;
	begin = clock();
	Input();
	Init(M);
	double k1, k2;
	if (T > 5000)
		k1 = 0.001, k2 = 0.001;
	else
		k1 = 0.01, k2 = 0.01;
	Genetic(T, M, (int)(k1*T), (int)(k2*T));
	Output(M);
	end = clock();
	printf("%.0f秒\n", (double)(end - begin) / CLOCKS_PER_SEC);
	return 0;
}

/*
                   _ooOoo_
                  o8888888o
                  88" . "88
                  (| -_- |)
                  O\  =  /O
               ____/`---'\____
             .'  \\|     |//  `.
            /  \\|||  :  |||//  \
           /  _||||| -:- |||||-  \
           |   | \\\  -  /// |   |
           | \_|  ''\---/''  |   |
           \  .-\__  `-`  ___/-. /
         ___`. .'  /--.--\  `. . ___
      ."" '<  `.___\_<|>_/___.'  >' "".
     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
     \  \ `-.   \_ __\ /__ _/   .-` /  /
======`-.____`-.___\_____/___.-`____.-'======
                   `=---='
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       I have a dream!A AC deram!!
orz orz orz orz orz orz orz orz orz orz orz
orz orz orz orz orz orz orz orz orz orz orz
orz orz orz orz orz orz orz orz orz orz orz

*/