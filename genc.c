#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

//
#define PI               3.14159265f
//
#define BITMASK(b)       (1 << ((b) % CHAR_BIT))
#define BITSLOT(b)       ((b) / CHAR_BIT)
#define BITSET(a, b)     ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b)   ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b)    ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb)    ((nb + CHAR_BIT - 1) / CHAR_BIT)
//
#define N_GENES          2

typedef struct {
     double min;
     double max;
     unsigned int bits;
     unsigned int offset;
} gene_t;

typedef struct {
     double prob_couple;
     double prob_mutation;
     gene_t genes[N_GENES];
     unsigned int size;
     unsigned int spec_bits;
     char **spec_pool;
     char *best_spec;
     unsigned int generation_cnt;
     unsigned int mutations_cnt;
     unsigned int coupling_cnt;
} population_t;

FILE *dbg_cuts;
FILE *dbg_couples;

double get_rand_normal() {
     return ((double)rand() + 1.0f)/((double)RAND_MAX + 1.0f);
}
void start_gene_pool(population_t *p) {
     int i, j;

     p->spec_pool = (char **)malloc(sizeof(char *)*p->size);
     for(i=0; i<p->size; i++) {
          p->spec_pool[i] = (char *)malloc(BITNSLOTS(p->spec_bits));
          for(j=0; j<p->spec_bits; j++) {
               if(rand() % 2)
                    BITSET(p->spec_pool[i], j);
               else
                    BITCLEAR(p->spec_pool[i], j);
          }
     }

     p->best_spec = (char *)malloc(BITNSLOTS(p->spec_bits));
     memset(p->best_spec, 0, BITNSLOTS(p->spec_bits));
}
void end_gene_pool(population_t *p) {
     int i;

     for(i=0; i<p->size; i++)
          free(p->spec_pool[i]);
     free(p->spec_pool);
     free(p->best_spec);
}
double get_gene_value(char *spec, gene_t gene) {
     int i = gene.offset;
     int end = gene.offset + gene.bits;
     double tmp = 0.0f;

     for(; i<end; i++)
          if(BITTEST(spec, i))
               tmp += (1 << (i-gene.offset));
     return gene.min + tmp*(gene.max - gene.min)/(pow(2, gene.bits) - 1);
}
double get_fitness(population_t *p, char *spec) {
     double gene_value[N_GENES];
     int i;

     for(i=0; i<N_GENES; i++)
          gene_value[i] = get_gene_value(spec, p->genes[i]);
     return (21.5f + gene_value[0]*sin(4.0f*PI*gene_value[0]) 
               + gene_value[1]*sin(20.0f*PI*gene_value[1]));
}
void save_population_king(population_t *p) {
     double tmp = 0.0f;
     double max_fitness = get_fitness(p, p->spec_pool[0]);
     double curr_fitness = get_fitness(p, p->best_spec);
     int max_idx = -1; 
     int i;

     for(i=0; i<p->size; i++) {
          tmp = get_fitness(p, p->spec_pool[i]);
          if((tmp > max_idx) && (tmp > curr_fitness)) {
               max_fitness = tmp;
               max_idx = i;
          }
     }

     if(max_idx > 0)
          memcpy(p->best_spec, p->spec_pool[max_idx], BITNSLOTS(p->spec_bits));
}
void couple_spec(population_t *p, char *spec1, char *spec2) {
     int cut_pos = rand() % (p->spec_bits);
     char *tmp = (char *)malloc(BITNSLOTS(p->spec_bits));
     int i;

     // PARENTS:
     //   g1 = [ a | b ]
     //   g2 = [ c | d ]
     //   len(a) = len(c) = cut_pos (bits)
     //   len(b) = len(d) = p->spec_bits - cut_pos (bits)
     //
     // CHILDS:
     //   g1 = [ a | d ]
     //   g2 = [ c | b ]
     //
     memcpy(tmp, spec1, BITNSLOTS(p->spec_bits));
     for(i=cut_pos; i<p->spec_bits; i++) {
          if(BITTEST(spec2, i))
               BITSET(spec1, i);
          else
               BITCLEAR(spec1, i);

          if(BITTEST(tmp, i))
               BITSET(spec2, i);
          else
               BITCLEAR(spec2, i);
     }

     p->coupling_cnt++;
     free(tmp);
}
void mutate_population(population_t *p) {
     int i, j;

     for(i=0; i<p->size; i++) {
          for(j=0; j<p->spec_bits; j++) {
               if(get_rand_normal() < p->prob_mutation) {
                    if(BITTEST(p->spec_pool[i], j))
                         BITCLEAR(p->spec_pool[i], j);
                    else
                         BITSET(p->spec_pool[i], j);
                    p->mutations_cnt++;
               }
          }
     }
}
void evolve_population(population_t *p) {
     int spec_bytes = BITNSLOTS(p->spec_bits);
     double *spec_fitness = (double *)malloc(sizeof(double)*p->size);
     double *spec_cum = (double *)malloc(sizeof(double)*p->size);
     char *pool_tmp = (char *)malloc(spec_bytes*p->size);
     char *spec_couple = NULL;
     double scale_factor = 0.0f;
     double tmp = 0.0f;
     int i, j;

     for(i=0; i<p->size; i++) {
          spec_fitness[i] = get_fitness(p, p->spec_pool[i]);
          scale_factor += spec_fitness[i];
          memcpy((pool_tmp+(i*spec_bytes)), p->spec_pool[i], spec_bytes);
     }

     spec_cum[0] = spec_fitness[0]/scale_factor;
     for(i=1; i<p->size; i++)
          spec_cum[i] = (spec_fitness[i]/scale_factor) + spec_cum[i-1];

     for(i=0; i<p->size; i++) {
          tmp = get_rand_normal();
          for(j=0; j<p->size; j++) {
               if(spec_cum[j] > tmp) {
                    //fprintf(dbg_couples, "%d:\t%d\t: %lf\t : %lf\n", i, j, spec_fitness[j], spec_fitness[i]/scale_factor);
                    memcpy(p->spec_pool[i], (pool_tmp+(j*spec_bytes)), spec_bytes);
                    break;
               }
          }
     }

     for(i=0; i<p->size; i++) {
          tmp = get_rand_normal();
          if(tmp > p->prob_couple) {
               // Not selected to evolve
          } else {
               //
               if(spec_couple != NULL) {
                    // PARENTS: p->spec_pool[i], spec_couple
                    //
                    couple_spec(p, p->spec_pool[i], spec_couple);
                    // Reset couple and try to find a new pair
                    spec_couple = NULL;
               } else
                    spec_couple = p->spec_pool[i];
          }
     }

     p->generation_cnt++;
     free(spec_fitness);
     free(spec_cum);
     free(pool_tmp);
}
void print_spec(population_t *p, char *spec) {
     int i;

     printf("Specimen Information:\n");
     printf("\t genome=[ ");
     for(i=0; i<p->spec_bits; i++)
          printf("%d", (BITTEST(spec, i) ? 1 : 0));
     printf(" ]\n");
     printf("\t fitness=%lf\n", get_fitness(p, spec));
     printf("\t Genes Value:\n");
     for(i=0; i<N_GENES; i++)
          printf("\t\t gene[%d]=%lf\n", i, get_gene_value(spec, p->genes[i]));
}
void print_population(population_t *p) {
     int i;

     printf("Population Detailed Status:\n");
     for(i=0; i<p->size; i++)
          print_spec(p, p->spec_pool[i]);
}
void print_stats(population_t *p) {
     double total_fitness = 0.0f;
     double tmp = 0.0f;
     double max_fitness = get_fitness(p, p->spec_pool[0]);
     int max_idx = 0;
     double min_fitness = max_fitness;
     int i;

     printf("Population Status:\n");
     for(i=0; i<p->size; i++) {
          tmp = get_fitness(p, p->spec_pool[i]);
          total_fitness += tmp;
          if(tmp > max_fitness) {
               max_fitness = tmp;
               max_idx = i;
          }
          if(tmp < min_fitness)
               min_fitness = tmp;
     }

     printf("\t generation_cnt=%d \t coupling_cnt=%d \t mutations_cnt=%d\n",
               p->generation_cnt, p->coupling_cnt, p->mutations_cnt);
     printf("\t sum_fitness=%lf \t avg_fitness=%lf\n", total_fitness, total_fitness/p->size);
     printf("\t max_fitness=%lf \t [ x1=%lf\t \t x2=%lf ]\n", max_fitness,
               get_gene_value(p->spec_pool[max_idx], p->genes[0]),
               get_gene_value(p->spec_pool[max_idx], p->genes[1]));
     printf("\t min_fitness=%lf\n", min_fitness);
}
int main(int argc, char **argv) {
     population_t pop;
     int i;

     dbg_cuts = fopen("dbg_cuts.log", "w");
     dbg_couples = fopen("dbg_couples.log", "w");

     // GENE 0: [-3.0, 12.1]
     // Resolution: 0.0001 
     // Range: 12.1 - (-3.0) = 15.1
     // Bits: 2^17 < 15.1 x (1/0.0001) < 2^18 = 18 Bits
     pop.genes[0].min = -3.0f;
     pop.genes[0].max = 12.1f;
     pop.genes[0].bits = 18;
     pop.genes[0].offset = 0;

     // GENE 0: [4.1, 5.8]
     // Resolution: 0.0001 
     // Range: 5.8 - 4.1 = 1.7
     // Bits: 2^14 < 1.7 x (1/0.0001) < 2^15
     // 15 Bits
     pop.genes[1].min = 4.1f;
     pop.genes[1].max = 5.8f;
     pop.genes[1].bits = 15;
     pop.genes[1].offset = 18;

     pop.size = 20;
     pop.spec_bits = 18 + 15;
     pop.prob_couple = 0.4f;//0.25f;
     pop.prob_mutation = 0.01f;
     pop.generation_cnt = 0;
     pop.mutations_cnt = 0;
     pop.coupling_cnt = 0;

     srand(time(NULL));
     printf("--------------------\n");
     printf("\t population_size=%d\t spec_bits=%d\n", pop.size, pop.spec_bits);
     printf("\t coupling_probability=%lf\n", pop.prob_couple);
     printf("\t mutation_probability=%lf\n", pop.prob_mutation);
     printf("\t goal: \t x1=%lf \t x2=%lf\n \t fitness=%lf\n",
               9.6236f, 4.4278f,
               (21.5f + 9.6236f*sin(4.0f*PI*9.6236f) 
               + 4.4278f*sin(20.0f*PI*4.4278f)));
     printf("\t N_GENES=%d\n", N_GENES);
     printf("\t genes info:\n");
     for(i=0; i<N_GENES; i++) {
          printf("\t gene[%d]:\n", i);
          printf("\t\t min=%lf\n \t\t max=%lf\n \t\t bits=%d\n \t\t offset=%d\n",
               pop.genes[i].min, pop.genes[i].max, pop.genes[i].bits, pop.genes[i].offset);
     }

     printf("--------------------\n");
     start_gene_pool(&pop);
     print_population(&pop);
     print_stats(&pop);

     printf("--------------------\n");
     for(i=0; i<100000; i++) {
          evolve_population(&pop);
          mutate_population(&pop);
          save_population_king(&pop);
          //print_stats(&pop);
     }
     
     printf("--------------------\n");
     print_population(&pop);
     print_stats(&pop);
     printf("--------------------\n");
     printf("Population King:\n");
     print_spec(&pop, pop.best_spec);

     end_gene_pool(&pop);
     fclose(dbg_couples);
     fclose(dbg_cuts);
     return 0;
}
