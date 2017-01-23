#include "workflow.h"



double calculate_tau0(struct node *nodes, int nodes_length, double *max_cost, double *min_cost)
{
	int i, j;
	int max_index = 0, min_index = 0;

	*max_cost = 0.0;
	*min_cost = 0.0;

	for (i = 0; i < nodes_length; i++) {
		max_index = 0, min_index = 0;
		for (j = 1; j < nodes[i].num_instances; j++) {
			if (nodes[i].instances[j].cost < nodes[i].instances[min_index].cost)
				min_index = j;
			if (nodes[i].instances[max_index].cost < nodes[i].instances[j].cost)
				max_index = j;
		}
		(*max_cost) += nodes[i].instances[max_index].cost;
		(*min_cost) += nodes[i].instances[min_index].cost;
	}
	return (*min_cost) / (*max_cost);
}

void ants_initialize(struct ant *a, int nodes_length)
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
void solution_construction(struct ant *a, struct node *nodes, int nodes_length)
{
	int i, j, k, max_index, max_instance_index = 0, min_instance_index = 0;
	double q, sum = 0.0, probability, base;
	double eta, eta_c, eta_t;
	for (i = 0; i < nodes_length; i++) {
		q = u01();
		max_index = 0;
		for (j = 1; j < nodes[i].num_instances; j++) {
			if (nodes[i].instances[max_instance_index].cost < nodes[i].instances[j].cost)
				max_instance_index = j;
			if (nodes[i].instances[min_instance_index].cost < nodes[i].instances[j].cost)
				min_instance_index = j;
		}

		if (q < Q) {
			for (j = 1; j < nodes[i].num_instances; j++) {
				eta_c = (nodes[i].instances[max_instance_index].cost - nodes[i].instances[j].cost + 1) /
					(nodes[i].instances[max_instance_index].cost - nodes[i].instances[min_instance_index].cost + 1);
				if (nodes[i].instances[max_index].pheromone * pow(1 / eta_c, BETA)
					< nodes[i].instances[j].pheromone * pow(1 / eta_c, BETA))
					max_index = j;
			}
		}
		else {
			//否则就赌轮盘
			sum = 0.0;
			for (j = 0; j < CITY_NUM; j++) {
				if (1 == remain_city[j])
					continue;
				sum += pheromone[a->city[i - 1]][j] * pow(1.0 / distance[a->city[i - 1]][j], BETA);
				city_not_be_selected[length_city_not_be_selected++] = j;
			}

			probability = u01();

			base = 0.0;
			for (j = 0; j < length_city_not_be_selected; j++) {
				base += (pheromone[a->city[i - 1]][city_not_be_selected[j]]
					* pow(1.0 / distance[a->city[i - 1]][city_not_be_selected[j]], 2) / sum);
				if (probability < base) {
					a->city[i] = city_not_be_selected[j];
					remain_city[city_not_be_selected[j]] = 1;
					break;
				}
			}
		}
	}
	a->fitness = fitness_evaluation(a->city, distance);
}

int main() {
	struct ant ants[POP];    // 全局蚁群
	struct ant best_ant;     // 最优秀蚂蚁
	struct node *nodes = NULL;
	double time_constraint, tau0, max_cost, min_cost;
	int nodes_length, i, j, iter, ant_index;
	nodes = read_data("9.txt", &nodes_length, &time_constraint);
	//for (i = 0; i < nodes_length; i++) {
	//	printf("id: %d\n", nodes[i].id);
	//	for (j = 0; j < nodes[i].num_instances; j++) {
	//		printf("instance id:%d %lf %lf %lf\n", nodes[i].instances[j].id, nodes[i].instances[j].reliability, nodes[i].instances[j].time, nodes[i].instances[j].cost);
	//	}
	//	printf("\n");
	//}
	//for (i = 0; i < nodes_length; i++) {
	//if (nodes[0].pred == NULL) {
	//	j = 100;
	//}
	//else
	//	j = nodes[0].pred->id;
	//printf("%d", nodes[0].succ->id);
	//printf("%d %d %d", nodes[6].num_pred, nodes[6].pred[0]->id, nodes[6].pred[1]->id);
	
	tau0 = calculate_tau0(nodes, nodes_length, &max_cost, &min_cost);
	//printf("%lf", tau0);
	pheromone_initialize(nodes,nodes_length, tau0);  //tau0要怎么算？
	ants_initialize(ants, nodes_length);
	for (iter = 0; iter < ITERATION; iter++) {
		for (i = 0; i < nodes_length; i++) {
			ants[0].solutions[i].n = &nodes[i];
			ants[0].solutions[i].inst = &nodes[i].instances[0];
		}

		for (i = 0; i < POP; i++) {
			time_evaluation(ants + i, nodes, nodes_length);
			score_evaluation(ants + i, time_constraint, min_cost, max_cost);
		}
	}
	//for (i = 0; i < nodes_length; i++)
	//	printf("%lf\n", ants[0].solutions[i].ct);
	free_memory(nodes, nodes_length, ants);
	return 0;
}