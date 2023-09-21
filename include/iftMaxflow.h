//
// Created by jordao on 16/03/18.
//

// TODO DOCUMENTAR essa pagina e adicionar maxflow iterativo

#ifndef IFT_IFTMAXFLOW_H_H
#define IFT_IFTMAXFLOW_H_H

#include "ift/core/dtypes/GQueue.h"
#include "ift/core/dtypes/LabeledSet.h"
#include "iftMemory.h"
#include "iftFImage.h"
#include "iftMImage.h"
#include "iftRadiometric.h"

#define IFT_TERMINAL ( (iftMaxflowArc*) 1)
#define IFT_ORPHAN ( (iftMaxflowArc*) 2)


#define iftTryInsertGQueue(Q, index) if(Q->L.elem[index].color != IFT_GRAY) iftInsertGQueue(&Q, index);

typedef struct ift_maxflow_node iftMaxflowNode;
typedef struct ift_maxflow_arc iftMaxflowArc;
typedef struct ift_maxflow_nodeptr iftMaxflowNodePtr;


/**
 * @brief Arc between nodes of graph to compute maxflow algorithm
 * @author Jordao Bragantini
 * @date April, 2018
 * @struct head         Node the arc points to
 * @struct next         Next arc with the same originating node
 * @struct sister       Reverse arc
 * @struct r_cap        Residual capacity
 */

struct ift_maxflow_arc {
    iftMaxflowNode *head;
    iftMaxflowArc *next;
    iftMaxflowArc *sister;
    float r_cap;
};


/**
 * @brief Node of graph to compute maxflow algorithm
 * @author Jordao Bragantini
 * @date April, 2018
 * @struct first        First outcoming arc of arc linked list
 * @struct parent       Arc to node's parent
 * @struct next         Next active node (or to itself if it is the last node in the list)
 * @struct time_stamp   Time stamp show when *dist* was computed
 * @struct dist         Distance from terminal (source or sink)
 * @struct is_sink      Flag showing whether the node belongs to source or sink
 * @struct color        IFT_WHITE = node is new or was change, IFT_BLACK node already was use
 * in a maxflow iteration
 * @struct tr_cap       Tree residual capacity, if tr_cap > 0 then tr_cap is the residual capacity of the arc
 * from source, otherwise -tr_cap is the residual capacity of the arc from SINK
 */

struct ift_maxflow_node {
    iftMaxflowArc *first;
    iftMaxflowArc *parent;
    iftMaxflowNode *next;
    int	time_stamp;
    int	dist;
    int	is_sink;
    char color;
    float tr_cap;
    int index; /* <- */
};


/**
 * @brief Structure to produce a linked list of nodes
 * @author Jordao Bragantini
 * @date April, 2018
 * @struct node         Pointer to node
 * @struct next         Pointer to next node
 */

struct ift_maxflow_nodeptr {
    iftMaxflowNode *node;
    iftMaxflowNodePtr *next;
};


/**
 * @brief Graph used to compute and store information of maxflow algorithm
 * @author Jordao Bragantini
 * @date April, 2018
 * @struct nodes        Array of nodes (nodes are allocated contiguously)
 * @struct arcs         Array of arcs (arcs are allocated contiguously)
 * @struct orphan_first First of orphans nodes linked list
 * @struct orphan_last  Last of orphans nodes linked list
 * @struct arc_last     Index of last arc + 1
 * @struct arc_max      Size of *arcs* array
 * @struct node_last    Index of last node + 1
 * @struct node_max     Size of *nodes* array
 * @struct iteration    Count iterations of maxflow algorithms usage
 * @struct time         Records the time to save on the *time_stamp*
 * @struct queue_first  Firsts on the lists of active nodes
 * @struct queue_last   Lasts on the lists of active nodes
 * @struct flow         Total flow
 * @note                There are two active nodes queues. Active nodes are added to the end
 * of the second queue and read from the front of the first queue. If the first queue is empty,
 * it is replaced by the second queue (and the second queue becomes empty).
 */

typedef struct ift_maxflow_graph {
    iftMaxflowNode *nodes;
    iftMaxflowArc *arcs;
    iftMaxflowNodePtr *orphan_first;
    iftMaxflowNodePtr *orphan_last;
    int arc_last;
    int arc_max;
    int node_last;
    int node_max;
    int iteration;
    int time;
    iftMaxflowNode *queue_first[2];
    iftMaxflowNode *queue_last[2];
    float flow;
} iftMaxflowGraph;


iftMaxflowGraph* iftCreateMaxflowGraph(int node_max, int edge_max);

void iftDestroyMaxflowGraph(iftMaxflowGraph** graph);

void iftAddNodeMaxflow(iftMaxflowGraph *graph);

void iftAddEdgeMaxflow(iftMaxflowGraph* graph, int p, int q, float cap, float rev_cap);

void iftAddTWeightsMaxflow(iftMaxflowGraph *graph, int _p, float cap_source, float cap_sink);

void iftActivateMaxflowNode(iftMaxflowGraph* graph, iftMaxflowNode* p);

iftMaxflowNode* iftNextActivateMaxflowNode(iftMaxflowGraph* graph);

void iftSetOrphanFrontMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p);

void iftSetOrphanRearMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p);

void iftAugmentMaxflow(iftMaxflowGraph *graph, iftMaxflowArc *a);

void iftProcessSinkOrphanMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p);

void iftProcessSourceOrphanMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p);

void iftMaxflowInit(iftMaxflowGraph *graph);

float iftMaxflow(iftMaxflowGraph *graph);

int iftMaxflowNodeLabel(iftMaxflowGraph *graph, int p);

iftMaxflowGraph* iftMaxflowGraphFromGradient(const iftImage *gradient, iftLabeledSet* seeds, int beta);

iftMaxflowGraph* iftMaxflowGraphFromMImage(const iftMImage * mimg, iftLabeledSet* seeds, int beta);

iftMaxflowGraph* iftMaxflowGraphFromMImageProb(const iftMImage * mimg, const iftImage *regions, const iftFImage* obj_prob,
                                               const iftFImage* bkg_prob, double beta);

iftMaxflowGraph* iftMaxflowGraphWithObjmap(const iftMImage * mimg, const iftImage *objmap,
                                           iftLabeledSet* seeds, float alpha, int beta);

//! swig(newobject)
iftImage* iftGraphCut(const iftImage *gradient, iftLabeledSet *seeds, int beta);

//! swig(newobject)
iftImage* iftGraphCutFromMImage(const iftMImage * mimg, iftLabeledSet* seeds, int beta);

//! swig(newobject)
iftImage* iftGraphCutWithObjmap(const iftMImage * mimg, const iftImage *objmap, iftLabeledSet* seeds, float alpha, int beta);

iftImage *iftMMaxflowSegment(iftMaxflowGraph *graph, const iftMImage *mimg);

iftImage *iftMaxflowSegment(iftMaxflowGraph *graph, const iftImage *img);

// ----------------------------------------------------------------

void iftProcessSinkOrphanOptPathMaxflow(iftMaxflowGraph *graph, iftGQueue *Q, iftMaxflowNode *p);

void iftProcessSourceOrphanOptPathMaxflow(iftMaxflowGraph *graph, iftGQueue *Q, iftMaxflowNode *p);

iftMaxflowNode *iftGetActiveNodeMaxflow(iftMaxflowGraph *graph, iftGQueue *Q);

float iftOptimumPathMaxflow(iftMaxflowGraph *graph, const iftImage *gradient);

//! swig(newobject)
iftImage* iftOptimumPathGraphCut(const iftImage *gradient, iftLabeledSet* seeds, int beta);

//! swig(newobject)
double iftGraphCutBeta(iftMImage *mimg);

#endif //IFT_IFTMAXFLOW_H_H
