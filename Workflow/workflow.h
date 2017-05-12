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
	int pred_len, succ_len;
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

double calculate_est_and_best(struct node *nodes, int nodes_length, double time_constrain);

double calculate_tau0(struct node *nodes, int nodes_length, double *max_cost, double *min_cost);

void ants_apply_memory(struct ant *a, int nodes_length);

void solution_construction(struct ant *a, struct node *nodes, int nodes_length, double time_constraint, double min_makespan);

void local_update_pheromone(struct ant *a, double tau0, int nodes_length);

void global_update_pheromone(struct ant *best_ant, int nodes_length);

void pheromone_initialize(struct node *nodes, int nodes_length, double tau0);

int find_the_best_ant(struct ant *a);

void score_evaluation(struct ant *a, double time_constraint, double min_cost, double max_cost);

void cost_evaluation(struct ant *a, struct node *nodes, int nodes_length);

void time_evaluation(struct ant *a, struct node *nodes, int nodes_length);

void backward_time_evaluation(struct ant *a, struct node *nodes, int nodes_length);

double calculate_avg_min_time(struct node n);

struct node * read_data(char *file_name, int *num_nodes, double *time_constraint);

void free_memory(struct node *nodes, int length, struct ant *ants);

#endif