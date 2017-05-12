#include "workflow.h"

int main() {
	struct ant ants[POP];    // 全局蚁群
	struct ant best_ant;     // 最优秀蚂蚁
	struct node *nodes = NULL;
	double time_constraint, tau0, max_cost, min_cost, min_makespan;
	int nodes_length, i, iter, best_ant_index;


	nodes = read_data("DNA.txt", &nodes_length, &time_constraint);

	tau0 = calculate_tau0(nodes, nodes_length, &max_cost, &min_cost);
	min_makespan = calculate_est_and_best(nodes, nodes_length, time_constraint);
	pheromone_initialize(nodes, nodes_length, tau0);  //tau0要怎么算？
	ants_apply_memory(ants, nodes_length);
	best_ant.solutions = (struct solution *)malloc(sizeof(struct solution) * nodes_length);

	for (iter = 0; iter < ITERATION; iter++) {
		for (i = 0; i < POP; i++) {
			solution_construction(ants + i, nodes, nodes_length, time_constraint, min_makespan);
			time_evaluation(ants + i, nodes, nodes_length);
			cost_evaluation(ants + i, nodes, nodes_length);
			score_evaluation(ants + i, time_constraint, min_cost, max_cost);
			local_update_pheromone(&ants[i], tau0, nodes_length);

		}

		best_ant_index = find_the_best_ant(ants);
		if (iter == 0 || best_ant.score < ants[best_ant_index].score) {
			memcpy(best_ant.solutions, ants[best_ant_index].solutions, sizeof(solution)* nodes_length);
			best_ant.cost = ants[best_ant_index].cost;
			best_ant.score = ants[best_ant_index].score;
			best_ant.makespan = ants[best_ant_index].makespan;

		}
		global_update_pheromone(&best_ant, nodes_length);

	}

	printf("makespan: %lf, cost: %lf\n", best_ant.makespan, best_ant.cost);
	for (i = 0; i < nodes_length; i++) {
		printf("%d:%d %lf %lf %lf\n", i, best_ant.solutions[i].inst->id, best_ant.solutions[i].inst->time, best_ant.solutions[i].st, best_ant.solutions[i].ct);
	}
	free(best_ant.solutions);
	free_memory(nodes, nodes_length, ants);
	return 0;
}



/************************************************************************
* Function Name : calculate_est_and_best
* Create Date : 2017/2/8
* Author/Corporation : bigzhao
**
Description : 计算有向图的est 及反过来的 best
*
* Param :*nodes：节点数组
* nodes_length: 节点的数量
**
************************************************************************/
double calculate_est_and_best(struct node *nodes, int nodes_length, double time_constrain)
{
	int i, j, min_index;
	struct ant a, backward_a;
	double min_makespan;

	// 为蚂蚁的 solution[nodes_length] 申请内存
	a.solutions = (struct solution *)malloc(sizeof(struct solution) * nodes_length);
	if (a.solutions == NULL) {
		fprintf(stderr, "ants initialize failed.");
		exit(1);
	}
	backward_a.solutions = (struct solution *)malloc(sizeof(struct solution) * nodes_length);
	if (backward_a.solutions == NULL) {
		fprintf(stderr, "ants initialize failed.");
		exit(1);
	}

	for (i = 0; i < nodes_length; i++) {
		min_index = 0;
		for (j = 1; j < nodes[i].num_instances; j++) {
			if (nodes[i].instances[j].time < nodes[i].instances[min_index].time)
				min_index = j;
		}
		a.solutions[i].inst = &nodes[i].instances[min_index];
		backward_a.solutions[i].inst = &nodes[i].instances[min_index];
	}
	time_evaluation(&a, nodes, nodes_length);
	for (i = 0; i < nodes_length; i++) {
		nodes[i].est = a.solutions[i].st;
	}
	min_makespan = a.makespan;
	backward_time_evaluation(&backward_a, nodes, nodes_length);
	for (i = 0; i < nodes_length; i++) {
		nodes[i].best = backward_a.solutions[i].st;
	}
	printf("\n EST:");
	for (i = 0; i < nodes_length; i++) {
		printf(" %lf ", nodes[i].est);
	}
	printf("\n");
	printf("\n BEST:");
	for (i = 0; i < nodes_length; i++) {
		printf(" %lf ", nodes[i].best);
	}
	printf("\n");
	free(a.solutions);
	free(backward_a.solutions);
	return min_makespan;
}

/************************************************************************
* Function Name : calculate_tau0
* Create Date : 2017/2/8（注释时间）
* Author/Corporation : bigzhao
**
Description : 计算信息素的初始值，ε0，这里使用的是贪心取得的值
*
* Param :*nodes：节点数组
* nodes_length: 节点的数量
**
************************************************************************/
double calculate_tau0(struct node *nodes, int nodes_length, double *max_cost, double *min_cost)
{
	int i, j;
	int max_index = 0, min_index = 0;

	*max_cost = 0.0;
	*min_cost = 0.0;

	for (i = 0; i < nodes_length; i++) {
		//max_index = 0, min_index = 0;
		//for (j = 1; j < nodes[i].num_instances; j++) {
		//	if (nodes[i].instances[j].cost < nodes[i].instances[min_index].cost)
		//		min_index = j;
		//	if (nodes[i].instances[max_index].cost < nodes[i].instances[j].cost)
		//		max_index = j;
		//}
		(*max_cost) += nodes[i].max_cost;
		(*min_cost) += nodes[i].min_cost;
	}
	return (*min_cost) / (*max_cost);
}

void ants_apply_memory(struct ant *a, int nodes_length)
{
	int i;

	for (i = 0; i < POP; i++) {
		a[i].solutions = (struct solution *)malloc(sizeof(struct solution) * nodes_length);
		if (a[i].solutions == NULL) {
			fprintf(stderr, "ants initialize failed.");
			exit(1);
		}
	}
}

/************************************************************************
* Function Name : pheromone_initialize
* Create Date : 2017/1/16
* Author/Corporation : bigzhao
**
Description : 信息素初始化
*
* Param :
* tau0: 
**
************************************************************************/
void pheromone_initialize(struct node *nodes, int nodes_length, double tau0)
{
	int i, j;

	for (i = 0; i < nodes_length; i++) {
		for (j = 0; j < nodes[i].num_instances; j++) {
			nodes[i].instances[j].pheromone = tau0;
		}
	}
}

/************************************************************************
* Function Name : solution_construction
* Create Date : 2017/1/16
* Author/Corporation : bigzhao
**
Description : 蚂蚁的解的构造，按概率来操作，小于Q则挑选最大的，大于Q则
以赌轮盘的形式进行随机选择
*
* Param : *a: 蚁群
* distance: 距离矩阵
* pheromone：信息素矩阵
**
************************************************************************/
void solution_construction(struct ant *a, struct node *nodes, int nodes_length, double time_constraint, double min_makespan)
{
	int i, j, max_index;
	double q, sum = 0.0, probability, base, avg_min_time, SD, MAT;
	double *eta, eta_c, eta_t;
	for (i = 0; i < nodes_length; i++) {
		q = u01();
		eta = (double *)malloc(nodes[i].num_instances * sizeof(double));
		avg_min_time = calculate_avg_min_time(nodes[i]);
		SD = avg_min_time * time_constraint / min_makespan;
		MAT = max(abs(nodes[i].max_time - SD), abs(SD - nodes[i].min_time));
		if (q < Q) {
			max_index = 0;
			eta_c = (nodes[i].max_cost - nodes[i].instances[0].cost + 1) /
				(nodes[i].max_cost - nodes[i].min_cost + 1);
			//if (nodes[i].instances[0].time <= SD)
			//	eta_t = 1;
			//else
				eta_t = (MAT - abs(nodes[i].instances[0].time - SD) + 1) / (MAT + 1);
			eta[0] = 0.5 * (eta_c + eta_t);
			for (j = 1; j < nodes[i].num_instances; j++) {
				// 计算 eta
				eta_c = (nodes[i].max_cost - nodes[i].instances[j].cost + 1) /
					(nodes[i].max_cost - nodes[i].min_cost + 1);
				//if (nodes[i].instances[j].time <= SD)
				//	eta_t = 1;
				//else
					eta_t = (MAT - abs(nodes[i].instances[j].time - SD) + 1) / (MAT + 1);
				eta[j] = 0.5 * (eta_c + eta_t);

				if (nodes[i].instances[max_index].pheromone * pow(BETA, eta[max_index])
					< nodes[i].instances[j].pheromone * pow(BETA, eta[j]))
					max_index = j;
			}

			a->solutions[i].inst = &nodes[i].instances[max_index];
		}
		else {
			//否则就赌轮盘
			sum = 0.0;
			for (j = 0; j < nodes[i].num_instances; j++) {
				// 计算 eta
				eta_c = (nodes[i].max_cost - nodes[i].instances[j].cost + 1) /
					(nodes[i].max_cost - nodes[i].min_cost + 1);
				//if (nodes[i].instances[j].time <= SD)
				//	eta_t = 1;
				//else
					eta_t = (MAT - abs(nodes[i].instances[j].time - SD) + 1) / (MAT + 1);
				eta[j] = 0.5 * (eta_c + eta_t);
				sum += nodes[i].instances[j].pheromone * pow(BETA, eta[j]);
			}

			probability = u01();

			base = 0.0;
			for (j = 0; j < nodes[i].num_instances; j++) {
				base += nodes[i].instances[j].pheromone * pow(BETA, eta[j]) / sum;
				if (probability < base) {
					a->solutions[i].inst = &nodes[i].instances[j];
					break;
				}
			}
		}
		free(eta);
	}
	//for (i = 0; i < nodes_length; i++)
	//	printf(" %d ", a->solutions[i].inst->id);
	//printf("\n");
}

/************************************************************************
* Function Name : local_update_pheromone
* Create Date : 2017/2/7
* Author/Corporation : bigzhao
**
Description : 局部信息素更新
*
* Param : a: 需要更新的蚂蚁
* pheromone: 信息素矩阵
* tau0：初始值
**
************************************************************************/
void local_update_pheromone(struct ant *a,  double tau0, int nodes_length)
{
	int i;
	for (i = 0; i < nodes_length; i++) {
		a->solutions[i].inst->pheromone *= (1 - ROU);
		a->solutions[i].inst->pheromone += ROU * tau0;
	}
}

/************************************************************************
* Function Name : global_update_pheromone
* Create Date : 2017/1/16
* Author/Corporation : bigzhao
**
Description : 全局信息素更新，全局更新只针对于最好的蚂蚁
*
* Param : best_ant: 需要更新的蚂蚁
* pheromone: 信息素矩阵
**
************************************************************************/
void global_update_pheromone(struct ant *best_ant, int nodes_length)
{
	int i;
	int city_index = 0;
	for (i = 0; i < nodes_length; i++) {
		best_ant->solutions[i].inst->pheromone *= (1 - ALPHA);
		best_ant->solutions[i].inst->pheromone += ALPHA * best_ant->score;
	}
}


void free_memory(struct node *nodes, int length, struct ant *ants)
{
	int i;

	for (i = 0; i < length; i++) {
		free(nodes[i].instances);
	}
	free(nodes);
	for (i = 0; i < POP; i++) {
		free(ants[i].solutions);
	}
}

struct node * read_data(char *file_name, int *num_nodes, double *time_constraint)
{
	FILE *fp;
	int num_node, num_edges, i, j, temp, current, next;
	int max_cost_index, min_cost_index, max_time_index, min_time_index;
	double temp1;
	void *new_ptr = NULL;
	struct node *nodes;
	fp = fopen(file_name, "r");
	if (fp == NULL) {
		fprintf(stderr, "make sure the file exist");
		exit(1);
	}
	fscanf(fp, "%d", &num_node);
	fscanf(fp, "%d", &num_edges);
	*num_nodes = num_node;
	nodes = (struct node *) malloc(sizeof(struct node) * num_node);

	for (i = 0; i < num_node; i++) {
		nodes[i].id = i;
		//nodes[i].succ = (struct node**)malloc((num_node-1) * sizeof(struct node*));
		//nodes[i].pred = (struct node**)malloc((num_node - 1) * sizeof(struct node*));
		nodes[i].succ = (struct node**)malloc(1 * sizeof(struct node*));
		nodes[i].pred = (struct node**)malloc(1 * sizeof(struct node*));
		nodes[i].num_pred = 0;
		nodes[i].num_succ = 0;
		nodes[i].pred_len = 1;
		nodes[i].succ_len = 1;
	}
	for (i = 0; i < num_edges; i++) {
		fscanf(fp, "%d %d %d", &temp, &current, &next);
		// 检查后继是否越界，如果是则realloc内存,为了方便每次扩大两倍 应该能够满足需求
		while (nodes[current - 1].succ_len <= nodes[current - 1].num_succ) {
			nodes[current - 1].succ_len *= 2;
			new_ptr = realloc(nodes[current - 1].succ, sizeof(struct node *) * (nodes[current - 1].succ_len));
			if (!new_ptr) {
				fprintf(stderr, "succ realloc error!");
				exit(-1);
			}
			nodes[current - 1].succ = (struct node **)new_ptr;
		}
		nodes[current - 1].succ[nodes[current - 1].num_succ] = &nodes[next - 1];
		nodes[current - 1].num_succ++;

		// 检查前驱是否越界，如果是则realloc内存,为了方便每次扩大两倍 应该能够满足需求
		while (nodes[next - 1].pred_len <= nodes[next - 1].num_pred) {
			nodes[next - 1].pred_len *= 2;
			new_ptr = realloc(nodes[next - 1].pred, sizeof(struct node *) * nodes[next - 1].pred_len);
			if (!new_ptr) {
				fprintf(stderr, "pred realloc error!");
				exit(-1);
			}
			nodes[next - 1].pred = (struct node **)new_ptr;
		}
		nodes[next - 1].pred[nodes[next - 1].num_pred] = &nodes[current - 1];
		nodes[next - 1].num_pred++;
	}

	for (i = 0; i < num_node; i++) {
		fscanf(fp, "%d %d", &temp, &nodes[i].num_instances);
		nodes[i].instances = (struct instance *) malloc(sizeof(struct instance) * nodes[i].num_instances);
		max_cost_index = 0; min_cost_index = 0; max_time_index = 0; min_time_index = 0;
		for (j = 0; j < nodes[i].num_instances; j++) {
			fscanf(fp, "%d %lf %lf %lf", &nodes[i].instances[j].id,
				&nodes[i].instances[j].reliability, &nodes[i].instances[j].time, &nodes[i].instances[j].cost);
			if (0 == j)
				continue;     //因为下面四个index默认值为0 所以没必要再去判断
			if (nodes[i].instances[max_cost_index].cost < nodes[i].instances[j].cost)
				max_cost_index = j;
			if (nodes[i].instances[j].cost < nodes[i].instances[min_cost_index].cost)
				min_cost_index = j;
			if (nodes[i].instances[max_time_index].time < nodes[i].instances[j].time)
				max_time_index = j;
			if (nodes[i].instances[j].time < nodes[i].instances[min_time_index].time)
				min_time_index = j;
		}
		// 赋值
		nodes[i].max_time = nodes[i].instances[max_time_index].time;
		nodes[i].min_time = nodes[i].instances[min_time_index].time;
		nodes[i].max_cost = nodes[i].instances[max_cost_index].cost;
		nodes[i].min_cost = nodes[i].instances[min_cost_index].cost;
		printf("/n %lf %lf %lf %lf \n", nodes[i].max_cost, nodes[i].min_cost, nodes[i].max_time, nodes[i].min_time);
	}
	fscanf(fp, "%lf %lf %lf", &temp1, time_constraint, &temp1);
	fclose(fp);
	return nodes;
}

double calculate_avg_min_time(struct node n)
{
	int min_est_succ_index = 0, min_best_pred_index = 0;
	int i;
	if (n.num_succ == 0) {
		for (i = 1; i < n.num_pred; i++) {
			if (n.pred[i]->best < n.pred[min_best_pred_index]->best)
				min_best_pred_index = i;
		}
		return (n.pred[min_best_pred_index]->best - n.best);
	}

	if (n.num_pred == 0) {
		for (i = 1; i < n.num_succ; i++) {
			if (n.succ[i]->est < n.succ[min_est_succ_index]->est)
				min_est_succ_index = i;
		}
		return (n.succ[min_est_succ_index]->est - n.est);
	}

	for (i = 1; i < n.num_succ; i++) {
		if (n.succ[i]->est < n.succ[min_est_succ_index]->est)
			min_est_succ_index = i;
	}
	for (i = 1; i < n.num_pred; i++) {
		if (n.pred[i]->best < n.pred[min_best_pred_index]->best)
			min_best_pred_index = i;
	}
	return (n.succ[min_est_succ_index]->est - n.est + n.pred[min_best_pred_index]->best - n.best) / 2.0;
}

/************************************************************************
* Function Name : backward_time_evaluation
* Create Date : 2017/2/7
* Author/Corporation : bigzhao
**
Description : 测量反方向的时间，也就是论文里的BEST，与上面函数不一样的地方
在于这里的pred和succ颠倒即可
*
* Param :
* :
**
************************************************************************/
void backward_time_evaluation(struct ant *a, struct node *nodes, int nodes_length)
{
	std::queue<struct node> ready_queue;
	int i, max_index;
	struct node operating_node;
	int *push_on = NULL;    //记录队列状态，节点是否已经进队列

	push_on = (int *)malloc(sizeof(int)* nodes_length);
	memset(push_on, 0, sizeof(int)* nodes_length);
	for (i = 0; i < nodes_length; i++) {
		if (0 == nodes[i].num_succ)
			ready_queue.push(nodes[i]);
	}

	while (!ready_queue.empty()) {
		operating_node = ready_queue.front();
		if (0 == operating_node.num_succ)
			a->solutions[operating_node.id].ct = 0;
		else {
			max_index = 0;
			for (i = 1; i < operating_node.num_succ; i++) {
				if (a->solutions[operating_node.succ[max_index]->id].st < a->solutions[operating_node.succ[i]->id].st) {
					max_index = i;
				}
			}
			a->solutions[operating_node.id].ct = a->solutions[operating_node.succ[max_index]->id].st;
		}

		a->solutions[operating_node.id].st = a->solutions[operating_node.id].ct + a->solutions[operating_node.id].inst->time;
		ready_queue.pop();
		for (i = 0; i < operating_node.num_pred; i++) {
			if (0 == push_on[operating_node.pred[i]->id]) {
				ready_queue.push(*(operating_node.pred[i]));
				push_on[operating_node.pred[i]->id] = 1;
			}
		}
		//operating_node.succ;
		//}

	}

	//max_index = 0;
	//for (i = 1; i < nodes_length; i++) {
	//	if (a->solutions[max_index].ct < a->solutions[i].ct)
	//		max_index = i;
	//	//printf("%lf  ", a->solutions[i].ct);
	//}
	//// 给蚂蚁a 赋 makespan
	//a->makespan = a->solutions[max_index].ct;
}

/************************************************************************
* Function Name : time_evaluation
* Create Date : 2017/2/7
* Author/Corporation : bigzhao
**
Description : 测量时间，正常方向
*
* Param :
* :
**
************************************************************************/
void time_evaluation(struct ant *a, struct node *nodes, int nodes_length)
{
	std::queue<struct node> ready_queue;
	int i, max_index;
	struct node operating_node;
	int *push_on;    //记录队列状态，节点是否已经进队列

	push_on = (int *)malloc(sizeof(int)* nodes_length);
	memset(push_on, 0, sizeof(int)* nodes_length);

	for (i = 0; i < nodes_length; i++) {
		if (0 == nodes[i].num_pred)
			ready_queue.push(nodes[i]);
	}

	while (!ready_queue.empty()) {
		operating_node = ready_queue.front();
		if (0 == operating_node.num_pred)
			a->solutions[operating_node.id].st = 0;
		else {
			max_index = 0;
			for (i = 1; i < operating_node.num_pred; i++) {
				if (a->solutions[operating_node.pred[max_index]->id].ct < a->solutions[operating_node.pred[i]->id].ct) {
					max_index = i;
				}
			}
			a->solutions[operating_node.id].st = a->solutions[operating_node.pred[max_index]->id].ct;
		}

		a->solutions[operating_node.id].ct = a->solutions[operating_node.id].st + a->solutions[operating_node.id].inst->time;
		ready_queue.pop();
		for (i = 0; i < operating_node.num_succ; i++) {
			if (0 == push_on[operating_node.succ[i]->id]) {
				ready_queue.push(*(operating_node.succ[i]));
				push_on[operating_node.succ[i]->id] = 1;
			}
		}
		//operating_node.succ;
		//}

	}

	max_index = 0;
	for (i = 1; i < nodes_length; i++) {
		if (a->solutions[max_index].ct < a->solutions[i].ct)
			max_index = i;
		//printf("%lf  ", a->solutions[i].ct);
	}
	// 给蚂蚁a 赋 makespan
	a->makespan = a->solutions[max_index].ct;
	//printf("%lf\n", a->makespan);
}


/************************************************************************
* Function Name : cost_evaluation
* Create Date : 2017/2/7
* Author/Corporation : bigzhao
**
Description : 测量cost
*
* Param :
* :
**
************************************************************************/
void cost_evaluation(struct ant *a, struct node *nodes, int nodes_length)
{
	int i;
	double cost = 0;
	for (i = 0; i < nodes_length; i++) {
		cost += a->solutions[i].inst->cost;
	}
	a->cost = cost;
}

void score_evaluation(struct ant *a, double time_constraint, double min_cost, double max_cost)
{
	if (a->makespan <= time_constraint)
		a->score = 1 + min_cost / a->cost;
	else
		a->score = time_constraint / a->makespan + min_cost / max_cost;
}

/************************************************************************
* Function Name : find_the_best_ant
* Create Date : 2017/1/16
* Author/Corporation : bigzhao
**
Description : 用遍历的方法找蚁群里面最优秀的蚂蚁
*
* Param : *a: 蚁群
* return: 最优的下标
**
************************************************************************/
int find_the_best_ant(struct ant *a)
{
	int i, best_ant_index = 0;

	for (i = 1; i < POP; i++) {
		if (a[best_ant_index].score < a[i].score)
			best_ant_index = i;
	}
	return best_ant_index;
}