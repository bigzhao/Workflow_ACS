#ifndef __WORKFLOW__
#define __WORKFLOW__

#define ROU       0.1     // ρ
#define POP       1      // 蚁群的种群大小
#define ITER      2500    // 迭代次数
#define BETA      2       // β
#define Q         0.9     // 概率
#define ALPHA     0.1     // α
#define ITERATION 100

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

void score_evaluation(struct ant *a, double time_constraint, double min_cost, double max_cost)
{
	if (a->makespan <= time_constraint)
		a->score = 1 + min_cost / a->cost;
	else
		a->score = time_constraint / a->makespan + min_cost / max_cost;
}

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
}

struct node * read_data(char *file_name, int *num_nodes, double *time_constraint) {
	FILE *fp;
	int num_node, num_edges, i, j, temp, current, next;
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
		//if (NULL = nodes[current].succ)
		//	nodes[current].succ++;
		//if (nodes[current - 1].num_succ == 0) {
		//	nodes[current - 1].succ = (struct node**)malloc(sizeof(struct node*));
		//	nodes[current - 1].succ[0] = &nodes[next - 1];
		//	nodes[current - 1].num_succ++;
		//}
		//else
		//{
		//	nodes[current - 1].num_succ++;
		//	nodes[current - 1].succ = (struct node**)realloc(nodes[current - 1].succ, sizeof(struct node*));
		//	nodes[current - 1].succ[nodes[current - 1].num_succ - 1] = &nodes[next - 1];
		//}
		//if (NULL != nodes[current].pred)
		//	nodes[current].pred++;
		//if (nodes[current - 1].num_pred == 0) {
		//	nodes[current - 1].pred = (struct node**)malloc(sizeof(struct node*));
		//	nodes[next - 1].pred[0] = &nodes[current - 1];
		//	nodes[next - 1].num_pred++;
		//}
		//else
		//{
		//	nodes[current - 1].num_pred++;
		//	nodes[current - 1].pred = (struct node**)realloc(nodes[current - 1].pred, sizeof(struct node*) * nodes[current - 1].num_pred);
		//	nodes[current - 1].pred[nodes[current - 1].num_pred - 1] = &nodes[next - 1];
		//}
	}

	for (i = 0; i < num_node; i++) {
		fscanf(fp, "%d %d", &temp, &nodes[i].num_instances);
		printf("%d %d\n", temp, nodes[i].num_instances);

		nodes[i].instances = (struct instance *) malloc(sizeof(struct instance) * nodes[i].num_instances);

		for (j = 0; j < nodes[i].num_instances; j++) {
			fscanf(fp, "%d %lf %lf %lf", &nodes[i].instances[j].id,
				&nodes[i].instances[j].reliability, &nodes[i].instances[j].time, &nodes[i].instances[j].cost);
		}
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