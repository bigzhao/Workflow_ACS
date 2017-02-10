#ifndef __WORKFLOW__
#define __WORKFLOW__

#define ROU       0.1     // ρ
#define POP       30      // 蚁群的种群大小
#define BETA      1.2      // β
#define Q         0.9     // 概率
#define ALPHA     0.1     // α
#define ITERATION 400

#define max(a,b)         (((a) > (b)) ? (a) : (b))

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>

boost::mt19937 rng((unsigned)time(0));
boost::uniform_01<boost::mt19937&> u01(rng);

struct ant {
	struct solution *solutions = NULL;
	double cost;
	double makespan;
	double score;
	//double fitness;
	//int start;
};

struct node {
	struct node **succ;
	struct node **pred;
	struct instance *instances = NULL;
	int num_instances;
	int num_succ = 0;
	int num_pred = 0;
	int id;
	double est;
	double best;
	double max_time, min_time, max_cost, min_cost;
};

struct instance {
	double reliability;
	double time;
	double cost;
	int id;
	double pheromone;
};

struct solution {
	struct node *n;
	struct instance *inst;
	double st;
	double ct;
};


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


void score_evaluation(struct ant *a, double time_constraint, double min_cost, double max_cost)
{
	if (a->makespan <= time_constraint)
		a->score = 1 + min_cost / a->cost;
	else
		a->score = time_constraint / a->makespan + min_cost / max_cost;
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
		for (i = 0; i < operating_node.num_succ; i++)
			ready_queue.push(*(operating_node.succ[i]));
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

	for (i = 0; i < nodes_length; i++) {
		if (0 == nodes[i].num_succ)
			ready_queue.push(nodes[i]);
	}

	while (!ready_queue.empty()) {
		operating_node = ready_queue.front();
		if (0 == operating_node.num_succ)
			a->solutions[operating_node.id].st = 0;
		else {
			max_index = 0;
			for (i = 1; i < operating_node.num_succ; i++) {
				if (a->solutions[operating_node.succ[max_index]->id].ct < a->solutions[operating_node.succ[i]->id].ct) {
					max_index = i;
				}
			}
			a->solutions[operating_node.id].st = a->solutions[operating_node.succ[max_index]->id].ct;
		}

		a->solutions[operating_node.id].ct = a->solutions[operating_node.id].st + a->solutions[operating_node.id].inst->time;
		ready_queue.pop();
		for (i = 0; i < operating_node.num_pred; i++)
			ready_queue.push(*(operating_node.pred[i]));
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

struct node * read_data(char *file_name, int *num_nodes, double *time_constraint) {
	FILE *fp;
	int num_node, num_edges, i, j, temp, current, next;
	int max_cost_index, min_cost_index, max_time_index, min_time_index;
	double temp1;
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
		nodes[i].succ = (struct node**)malloc((num_node-1) * sizeof(struct node*));
		nodes[i].pred = (struct node**)malloc((num_node - 1) * sizeof(struct node*));
		nodes[i].num_pred = 0;
		nodes[i].num_succ = 0;
	}
	for (i = 0; i < num_edges; i++) {
		fscanf(fp, "%d %d %d", &temp, &current, &next);
		nodes[current - 1].succ[nodes[current - 1].num_succ] = &nodes[next - 1];
		nodes[current - 1].num_succ++;
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

#endif