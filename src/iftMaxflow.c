//
// Created by jordao on 16/03/18.
//

#include "iftMaxflow.h"

#include "ift/core/io/Stream.h"

iftMaxflowGraph *iftCreateMaxflowGraph(int node_max, int edge_max){
    iftMaxflowGraph *graph = NULL;

    graph = iftAlloc(1, sizeof *graph);

    graph->nodes = iftAlloc(node_max, sizeof *graph->nodes);
    graph->arcs = iftAlloc(edge_max,  sizeof *graph->arcs);

    if (!graph->nodes || !graph->arcs)
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMaxflowGraph");

    graph->node_last = 0;
    graph->node_max = node_max;

    graph->arc_last = 0;
    graph->arc_max = edge_max;

    graph->iteration = 0;
    graph->flow = 0.0;

    return graph;
}


void iftDestroyMaxflowGraph(iftMaxflowGraph** graph){
    iftMaxflowGraph *aux = *graph;
    if(aux != NULL) {
        iftFree(aux->nodes);
        iftFree(aux->arcs);
        iftFree(aux);
        *graph = NULL;
    }
}


void iftAddNodeMaxflow(iftMaxflowGraph* graph){

    if(graph->node_last == graph->node_max) iftError("Maxflow graph is full!", "iftAddNodeMaxflow");

    iftMaxflowNode* last = &graph->nodes[graph->node_last];

    last->first = NULL;
    last->tr_cap = 0;
    last->color = IFT_WHITE;

    graph->node_last++;
}


void iftAddEdgeMaxflow(iftMaxflowGraph* graph, int _p, int _q, float cap, float rev_cap){

    if( (_p == _q) || (_p < 0) || (_p > graph->node_last) || (_q < 0) || (_q > graph->node_last) || (cap < 0) || (rev_cap < 0))
        iftError("Invalid Edge", "iftAddEdgeMaxflow");

    // TODO CREATE REALLOC ARCS

    if(graph->arc_last == graph->arc_max)
        iftError("Cannot add more arcs, need to alloc more memory", "iftAddEdgeMaxflow");

    iftMaxflowArc *a = &graph->arcs[graph->arc_last];
    graph->arc_last++;
    iftMaxflowArc *a_rev = &graph->arcs[graph->arc_last];
    graph->arc_last++;

    iftMaxflowNode *p = &graph->nodes[_p];
    iftMaxflowNode *q = &graph->nodes[_q];

    a->sister = a_rev;
    a_rev->sister = a;

    a->next = p->first;
    p->first = a;

    a_rev->next = q->first;
    q->first = a_rev;

    a->head = q;
    a_rev->head = p;

    a->r_cap = cap;
    a_rev->r_cap = rev_cap;
}


void iftAddTWeightsMaxflow(iftMaxflowGraph *graph, int p, float cap_source, float cap_sink){

    if( (p < 0) || (p > graph->node_last) )
        iftError("Invalid Edge", "iftAddTWeightsMaxflow");

    float delta = graph->nodes[p].tr_cap;
    if (delta > 0) cap_source += delta;
    else           cap_sink -= delta;
    graph->flow += (cap_source < cap_sink) ? cap_source : cap_sink;
    graph->nodes[p].tr_cap = cap_source - cap_sink;
}


void iftActivateMaxflowNode(iftMaxflowGraph* graph, iftMaxflowNode* node){

    if(!node->next){
        if(graph->queue_last[1]) graph->queue_last[1]->next = node;
        else                     graph->queue_first[1] = node;
        graph->queue_last[1] = node;
        node->next = node;
    }

}


iftMaxflowNode* iftNextActivateMaxflowNode(iftMaxflowGraph* graph){
    iftMaxflowNode *p;

    while(true){
        p = graph->queue_first[0];

        if(!p){
            p = graph->queue_first[1];
            graph->queue_first[0] = p;
            graph->queue_last[0] = graph->queue_last[1];
            graph->queue_first[1] = NULL;
            graph->queue_last[1] = NULL;
            if (!p) return NULL;
        }

        if(p->next == p){
            graph->queue_last[0] = NULL;
            graph->queue_first[0] = NULL;
        } else {
            graph->queue_first[0] = p->next;
        }
        p->next = NULL;

        if(p->parent != NULL) return p;
    }
}


void iftSetOrphanFrontMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p){
    iftMaxflowNodePtr* ptr;
    p->parent = IFT_ORPHAN;
    ptr = iftAlloc(1, sizeof *ptr);
    ptr->node = p;
    ptr->next = graph->orphan_first;
    graph->orphan_first = ptr;
}


void iftSetOrphanRearMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p){
    iftMaxflowNodePtr* ptr;
    p->parent = IFT_ORPHAN;
    ptr = iftAlloc(1, sizeof *ptr);
    ptr->node = p;
    if(graph->orphan_last != NULL) graph->orphan_last->next = ptr;
    else                           graph->orphan_first = ptr;
    graph->orphan_last = ptr;
    ptr->next = NULL;
}


void iftAugmentMaxflow(iftMaxflowGraph *graph, iftMaxflowArc *middle_arc){
    iftMaxflowNode* p;
    iftMaxflowArc * a;
    float bottleneck;

    /* 1. Finding bottleneck capacity
     * 1a - the source tree
     */
    bottleneck = middle_arc->r_cap;
    for(p = middle_arc->sister->head; ; p = a->head){
        a = p->parent;
        if(a == IFT_TERMINAL) break;
        if(bottleneck > a->sister->r_cap) bottleneck = a->sister->r_cap;
    }
    if(bottleneck > p->tr_cap) bottleneck = p->tr_cap;
    /* 1b - the sink tree */
    for(p = middle_arc->head; ; p = a->head){
        a = p->parent;
        if(a == IFT_TERMINAL) break;
        if(bottleneck > a->r_cap) bottleneck = a->r_cap;
    }
    if(bottleneck > -p->tr_cap) bottleneck = -p->tr_cap;

    /* 2. Augmenting
     * 2a - the source tree
     */
    middle_arc->sister->r_cap += bottleneck;
    middle_arc->r_cap -= bottleneck;
    for(p = middle_arc->sister->head; ; p = a->head){
        a = p->parent;
        if(a == IFT_TERMINAL) break;
        a->r_cap += bottleneck;
        a->sister->r_cap -= bottleneck;
        if(a->sister->r_cap == 0){
            iftSetOrphanFrontMaxflow(graph, p);
        }
    }
    p->tr_cap -= bottleneck;
    if(p->tr_cap == 0){
        iftSetOrphanFrontMaxflow(graph, p);
    }
    /* 2b - the sink tree */
    for(p = middle_arc->head; ; p = a->head){
        a = p->parent;
        if(a == IFT_TERMINAL) break;
        a->sister->r_cap += bottleneck;
        a->r_cap -= bottleneck;
        if(a->r_cap == 0){
            iftSetOrphanFrontMaxflow(graph, p);
        }
    }
    p->tr_cap += bottleneck;
    if(p->tr_cap == 0){
        iftSetOrphanFrontMaxflow(graph, p);
    }

    graph->flow += bottleneck;
}


void iftProcessSinkOrphanMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p){

    iftMaxflowNode *q;
    iftMaxflowArc *a0, *a0_min = NULL, *a;
    int d, d_min = IFT_INFINITY_INT;

    for(a0 = p->first; a0 != NULL; a0 = a0->next) {
        if (a0->r_cap != 0) {
            q = a0->head;
            a = q->parent;
            if (q->is_sink && a != NULL) {
                d = 0;
                while (true) {
                    if (q->time_stamp == graph->time) {
                        d += q->dist;
                        break;
                    }
                    a = q->parent;
                    d++;
                    if (a == IFT_TERMINAL) {
                        q->time_stamp = graph->time;
                        q->dist = 1;
                        break;
                    }
                    if (a == IFT_ORPHAN) {
                        d = IFT_INFINITY_INT;
                        break;
                    }
                    q = a->head;
                }

                if (d < IFT_INFINITY_INT) { /* q originates from the sink - done */

                    if (d < d_min) {
                        a0_min = a0;
                        d_min = d;
                    }
                    /* set marks along the path */
                    for (q = a0->head; q->time_stamp != graph->time; q = q->parent->head) {
                        q->time_stamp = graph->time;
                        q->dist = d;
                        d--;
                    }
                }
            }
        }
    }

    p->parent = a0_min;
    if(p->parent != NULL){
        p->time_stamp = graph->time;
        p->dist = d_min + 1;
    } else {
        /* TODO addd to changed list */
        for (a0 = p->first; a0 != NULL; a0 = a0->next) {
            q = a0->head;
            a = q->parent;
            if (q->is_sink && a != NULL) {
                if (a0->r_cap != 0) iftActivateMaxflowNode(graph, q);
                if (a != IFT_TERMINAL && a != IFT_ORPHAN && a->head == p) {
                    iftSetOrphanRearMaxflow(graph, q);
                }
            }
        }
    }
}


void iftProcessSourceOrphanMaxflow(iftMaxflowGraph *graph, iftMaxflowNode *p){

    iftMaxflowNode *q;
    iftMaxflowArc *a0, *a0_min = NULL, *a;
    int d, d_min = IFT_INFINITY_INT;

    for(a0 = p->first; a0 != NULL; a0 = a0->next) {
        if (a0->sister->r_cap != 0) {
            q = a0->head;
            a = q->parent;
            if (!q->is_sink && a != NULL) {
                d = 0;
                while (true) {
                    if (q->time_stamp == graph->time) {
                        d += q->dist;
                        break;
                    }
                    a = q->parent;
                    d++;
                    if (a == IFT_TERMINAL) {
                        q->time_stamp = graph->time;
                        q->dist = 1;
                        break;
                    }
                    if (a == IFT_ORPHAN) {
                        d = IFT_INFINITY_INT;
                        break;
                    }
                    q = a->head;
                }

                if (d < IFT_INFINITY_INT) { /* q originates from the sink - done */

                    if (d < d_min) {
                        a0_min = a0;
                        d_min = d;
                    }
                    /* set marks along the path */
                    for (q = a0->head; q->time_stamp != graph->time; q = q->parent->head) {
                        q->time_stamp = graph->time;
                        q->dist = d;
                        d--;
                    }
                }
            }
        }
    }

    p->parent = a0_min;
    if(p->parent != NULL){
        p->time_stamp = graph->time;
        p->dist = d_min + 1;
    } else {
        /* TODO add to changed list */
        for(a0 = p->first; a0 != NULL; a0 = a0->next){
            q = a0->head;
            a = q->parent;
            if(!q->is_sink && a != NULL){
                if(a0->sister->r_cap != 0) iftActivateMaxflowNode(graph, q);
                if(a != IFT_TERMINAL && a != IFT_ORPHAN && a->head == p){
                    iftSetOrphanRearMaxflow(graph, q);
                }
            }
        }
    }
}


void iftMaxflowInit(iftMaxflowGraph *graph){

    graph->queue_first[0] = graph->queue_last[0] = NULL;
    graph->queue_first[1] = graph->queue_last[1] = NULL;
    graph->orphan_first = NULL;

    graph->time = 0;

    for(int i = 0; i < graph->node_last; i++){

        iftMaxflowNode *p = &graph->nodes[i];

        p->next = NULL;
        p->color = IFT_WHITE;
        p->time_stamp = graph->time;
        if (p->tr_cap > 0) {
            /* p is connected to the source */
            p->is_sink = 0;
            p->parent = IFT_TERMINAL;
            iftActivateMaxflowNode(graph, p);
            p->dist = 1;
        } else if (p->tr_cap < 0) {
            /* p is connected to the sink */
            p->is_sink = 1;
            p->parent = IFT_TERMINAL;
            iftActivateMaxflowNode(graph, p);
            p->dist = 1;
        } else {
            p->parent = NULL;
        }
    }
}


float iftMaxflow(iftMaxflowGraph *graph){

    iftMaxflowNode *p, *q, *current_node  = NULL;
    iftMaxflowArc *a;
    iftMaxflowNodePtr *ptr, *ptr_next;

    iftMaxflowInit(graph);

    int count = 0;

    while(true){

        count++;

        p = current_node;
        if(p){
            p->next = NULL; /* remove active flag */
            if(!p->parent) p = NULL;
        }
        if(!p){
            p = iftNextActivateMaxflowNode(graph);
            if(!p) break;
        }

        /* growth */
        if(!p->is_sink){
            /* grow source tree */
            for(a = p->first; a != NULL; a = a->next){
                if(a->r_cap != 0){
                    q = a->head;
                    if(!q->parent){
                        q->is_sink = 0;
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                        iftActivateMaxflowNode(graph, q);
                        /* TODO add to change list */
                    } else if (q->is_sink) {
                         break;
                    } else if (q->time_stamp <= p->time_stamp && q->dist > p->dist){
                        /* heuristic - trying to make the distance from j to the source shorter */
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                    }
                }
            }
        } else {
            /* grow sink tree */
            for(a = p->first; a != NULL; a = a->next){
                if(a->sister->r_cap != 0){
                    q = a->head;
                    if(!q->parent){
                        q->is_sink = 1;
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                        iftActivateMaxflowNode(graph, q);
                        /* TODO add to changed list */
                    } else if (!q->is_sink) {
                        a = a->sister;
                        break;
                    } else if (q->time_stamp <= p->time_stamp && q->dist > p->dist){
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                    }
                }
            }
        }

        graph->time++;

        if(a != NULL) {
            p->next = p; /* set active flag */
            current_node = p;

            iftAugmentMaxflow(graph, a);

            /* adoption */
            while( (ptr = graph->orphan_first) != NULL){
                ptr_next = ptr->next;
                ptr->next = NULL;

                while( (ptr = graph->orphan_first) != NULL){
                    graph->orphan_first = ptr->next;
                    p = ptr->node;
                    iftFree(ptr);
                    if(graph->orphan_first == NULL) graph->orphan_last = NULL;
                    if(p->is_sink) {
                        iftProcessSinkOrphanMaxflow(graph, p);
                    } else {
                        iftProcessSourceOrphanMaxflow(graph, p);
                    }
                }

                graph->orphan_first = ptr_next;
            }
            /* adoption end */
        } else {
            current_node = NULL;
        }

    }

//    printf("count %d\n", count);

    graph->iteration++;
    return graph->flow;
}


int iftMaxflowNodeLabel(iftMaxflowGraph *graph, int p){
    if(graph->nodes[p].parent != NULL){
        return (graph->nodes[p].is_sink) ? 0 : 1;
    } else {
        return 1;
    }
}


iftMaxflowGraph* iftMaxflowGraphFromGradient(const iftImage* gradient, iftLabeledSet* seeds, int beta){

    int n_edges = 0;
    iftAdjRel *A = NULL;

    if (iftIs3DImage(gradient)) {
        n_edges = 26 * gradient->xsize * gradient->ysize * gradient->zsize
                  - 18 * (gradient->xsize * gradient->ysize + gradient->xsize * gradient->zsize + gradient->ysize * gradient->zsize)
                  + 12 * (gradient->xsize + gradient->ysize + gradient->zsize) - 8;
        A = iftSpheric(sqrtf(3.0f));
    } else {
        n_edges = 8 * gradient->xsize * gradient->ysize - 6 * (gradient->xsize + gradient->ysize) + 4;
        A = iftCircular(sqrtf(2.0f));
    }

    iftMaxflowGraph *graph = iftCreateMaxflowGraph(gradient->n, n_edges);

    iftFImage *normalized_grad = iftCreateFImage(gradient->xsize, gradient->ysize, gradient->zsize);

    int p;
    float maxval = iftMaximumValue(gradient);
    for(p = 0; p < gradient->n; p++){
        normalized_grad->val[p] = gradient->val[p]/maxval;
        iftAddNodeMaxflow(graph);
    }

    for(p = 0; p < gradient->n; p++){
        iftVoxel u = iftGetVoxelCoord(gradient,p);

        int i;
        for(i = 1; i < A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if(iftValidVoxel(gradient, v)){
                int q = iftGetVoxelIndex(gradient,v);

                if(p < q){
                    float edge_weight = powf(1.0f - (normalized_grad->val[p] + normalized_grad->val[q])/2.0f, beta);
                    iftAddEdgeMaxflow(graph, p, q, edge_weight, edge_weight);
                }
            }
        }
    }

    for(iftLabeledSet *S = seeds; S != NULL; S = S->next){
        p = S->elem;
        int label = S->label;
        if(label != 0)
            iftAddTWeightsMaxflow(graph, p, IFT_INFINITY_FLT, 0);
        else
            iftAddTWeightsMaxflow(graph, p, 0, IFT_INFINITY_FLT);
    }

    iftDestroyAdjRel(&A);
    iftDestroyFImage(&normalized_grad);

    return graph;
}


iftMaxflowGraph* iftMaxflowGraphFromMImage(const iftMImage* mimg, iftLabeledSet* seeds, int beta){

    int n_edges = 0;
    iftAdjRel *A = NULL;

    if (iftIs3DMImage(mimg)) {
        n_edges = 26 * mimg->xsize * mimg->ysize * mimg->zsize
                  - 18 * (mimg->xsize * mimg->ysize + mimg->xsize * mimg->zsize + mimg->ysize * mimg->zsize)
                  + 12 * (mimg->xsize + mimg->ysize + mimg->zsize) - 8;
        A = iftSpheric(1.0f);
    } else {
        n_edges = 8 * mimg->xsize * mimg->ysize - 6 * (mimg->xsize + mimg->ysize) + 4;
        A = iftCircular(1.0f);
    }

    iftMaxflowGraph *graph = iftCreateMaxflowGraph(mimg->n, n_edges);

    int p;
    float maxval = -1;
    for(p = 0; p < mimg->n; p++){

        iftAddNodeMaxflow(graph);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for(int i = 1; i < A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if(iftMValidVoxel(mimg, v)){
                int q = iftMGetVoxelIndex(mimg, v);

                if(p < q){
                    float dist = iftMImageDist(mimg, p, q);
                    if (dist > maxval) maxval = dist;
                }
            }
        }
    }

    for(p = 0; p < mimg->n; p++){

        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for(int i = 1; i < A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if(iftMValidVoxel(mimg, v)){
                int q = iftMGetVoxelIndex(mimg, v);

                if(p < q){
                    /* fabs because of numerical error */
                    float edge_weight = powf(fabs(maxval - iftMImageDist(mimg, p, q)) / maxval, beta);
                    iftAddEdgeMaxflow(graph, p, q, edge_weight, edge_weight);
                }
            }
        }
    }

    for(iftLabeledSet *S = seeds; S != NULL; S = S->next){
        p = S->elem;
        int label = S->label;
        if(label != 0)
            iftAddTWeightsMaxflow(graph, p, IFT_INFINITY_FLT, 0);
        else
            iftAddTWeightsMaxflow(graph, p, 0, IFT_INFINITY_FLT);
    }

    iftDestroyAdjRel(&A);

    return graph;
}


iftMaxflowGraph* iftMaxflowGraphFromMImageProb(const iftMImage* mimg, const iftImage *regions, const iftFImage *obj_prob,
                                               const iftFImage *bkg_prob, double beta)
{
    int n_edges = 0;
    const int gamma = 50; /* fixed by grabcut paper */
    iftAdjRel *A = NULL;

    if (mimg->n != obj_prob->n || mimg->n != bkg_prob->n || regions->n != mimg->n) {
        iftError("Multi-band image and probability images have different size", "iftMaxflowGraphFromMImageProb");
    }

    if (iftIs3DMImage(mimg)) {
        n_edges = 26 * mimg->xsize * mimg->ysize * mimg->zsize
                  - 18 * (mimg->xsize * mimg->ysize + mimg->xsize * mimg->zsize + mimg->ysize * mimg->zsize)
                  + 12 * (mimg->xsize + mimg->ysize + mimg->zsize) - 8;
        A = iftSpheric(1.0f);
    } else {
        n_edges = 8 * mimg->xsize * mimg->ysize - 6 * (mimg->xsize + mimg->ysize) + 4;
        A = iftCircular(1.0f);
    }

    iftMaxflowGraph *graph = iftCreateMaxflowGraph(mimg->n, n_edges);

    for (int p = 0; p < mimg->n; p++)
        iftAddNodeMaxflow(graph);

    for (int p = 0; p < mimg->n; p++)
    {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if (iftMValidVoxel(mimg, v)) {
                int q = iftMGetVoxelIndex(mimg, v);

                if (p < q) {
                    /* fabs because of numerical error */
                    float edge_weight = gamma * expf(-beta * iftMImageDist(mimg, p, q));
                    iftAddEdgeMaxflow(graph, p, q, edge_weight, edge_weight);
                }
            }
        }
    }


//    SURE_BKG = 0, SURE_OBJ = 1, MAYBE_BKG = 2, MAYBE_OBJ = 3
    for (int p = 0; p < mimg->n; p++) {
        float p_obj;
        float p_bkg;
        switch (regions->val[p]) {
            case 0:
                iftAddTWeightsMaxflow(graph, p, 0.0, 9 * gamma);
                break;
            case 1:
                iftAddTWeightsMaxflow(graph, p, 9 * gamma, 0.0);
                break;
            default:
                p_obj = -logf(obj_prob->val[p]);
                p_bkg = -logf(bkg_prob->val[p]);
// Are the arc-weights allowed to be negative?
//                if (p_obj < 0)
//                    p_obj = 0;
//                if (p_bkg < 0)
//                    p_bkg = 0;
                iftAddTWeightsMaxflow(graph, p, p_bkg, p_obj);
                break;
        }
    }

    iftDestroyAdjRel(&A);

    return graph;
}



iftMaxflowGraph* iftMaxflowGraphWithObjmap(const iftMImage* mimg, const iftImage *objmap,
                                           iftLabeledSet* seeds, float alpha, int beta)
{
    int n_edges = 0;
    iftAdjRel *A = NULL;

    if (iftIs3DMImage(mimg)) {
        n_edges = 26 * mimg->xsize * mimg->ysize * mimg->zsize
                  - 18 * (mimg->xsize * mimg->ysize + mimg->xsize * mimg->zsize + mimg->ysize * mimg->zsize)
                  + 12 * (mimg->xsize + mimg->ysize + mimg->zsize) - 8;
        A = iftSpheric(1.0f);
    } else {
        n_edges = 8 * mimg->xsize * mimg->ysize - 6 * (mimg->xsize + mimg->ysize) + 4;
        A = iftCircular(1.0f);
    }

    iftMaxflowGraph *graph = iftCreateMaxflowGraph(mimg->n, n_edges);

    int p;
    float max_img_dist = -1;
    int max_obj_dist = -1;
    for(p = 0; p < mimg->n; p++)
    {
        iftAddNodeMaxflow(graph);
        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for(int i = 1; i < A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if(iftMValidVoxel(mimg, v)){
                int q = iftMGetVoxelIndex(mimg, v);

                if(p < q){
                    float img_dist = iftMImageDist(mimg, p, q);
                    if (img_dist > max_img_dist) max_img_dist = img_dist;
                    int obj_dist = abs(objmap->val[p] - objmap->val[q]);
                    if (obj_dist > max_obj_dist) max_obj_dist = obj_dist;
                }
            }
        }
    }

    for(p = 0; p < mimg->n; p++){

        iftVoxel u = iftMGetVoxelCoord(mimg, p);

        for(int i = 1; i < A->n; i++){
            iftVoxel v = iftGetAdjacentVoxel(A,u,i);

            if(iftMValidVoxel(mimg, v)){
                int q = iftMGetVoxelIndex(mimg, v);

                if(p < q) {
                    float img_dist = fabs(max_img_dist - iftMImageDist(mimg, p, q)) / max_img_dist; /* fabs because of numerical error */
                    float obj_dist = (max_obj_dist - abs(objmap->val[q] - objmap->val[p])) / max_obj_dist;
                    float edge_weight = powf( (img_dist) * (1 - alpha) + (obj_dist) * alpha, beta);
                    iftAddEdgeMaxflow(graph, p, q, edge_weight, edge_weight);
                }
            }
        }
    }

    for (iftLabeledSet *S = seeds; S != NULL; S = S->next) {
        p = S->elem;
        int label = S->label;
        if(label != 0)
            iftAddTWeightsMaxflow(graph, p, IFT_INFINITY_FLT, 0);
        else
            iftAddTWeightsMaxflow(graph, p, 0, IFT_INFINITY_FLT);
    }

    iftDestroyAdjRel(&A);

    return graph;
}


iftImage *iftMMaxflowSegment(iftMaxflowGraph *graph, const iftMImage *mimg)
{
    iftImage *seg = iftCreateImage(mimg->xsize, mimg->ysize, mimg->zsize);

    for (int p = 0; p < seg->n; p++) {
        if (iftMaxflowNodeLabel(graph, p) != 0)
            seg->val[p] = 1;
    }
    return seg;
}


iftImage *iftMaxflowSegment(iftMaxflowGraph *graph, const iftImage *img)
{
    iftImage *seg = iftCreateImage(img->xsize, img->ysize, img->zsize);

    for (int p = 0; p < seg->n; p++) {
        if (iftMaxflowNodeLabel(graph, p) != 0)
            seg->val[p] = 1;
    }
    return seg;
}


iftImage* iftGraphCut(const iftImage* gradient, iftLabeledSet* seeds, int beta){

    iftMaxflowGraph *graph = iftMaxflowGraphFromGradient(gradient, seeds, beta);

    iftMaxflow(graph);

    iftImage* seg = iftMaxflowSegment(graph, gradient);

    iftDestroyMaxflowGraph(&graph);

    return seg;
}


iftImage* iftGraphCutFromMImage(const iftMImage* mimg, iftLabeledSet* seeds, int beta){

    iftMaxflowGraph *graph = iftMaxflowGraphFromMImage(mimg, seeds, beta);

    iftMaxflow(graph);

    iftImage* seg = iftMMaxflowSegment(graph, mimg);

    iftDestroyMaxflowGraph(&graph);

    return seg;
}


iftImage* iftGraphCutWithObjmap(const iftMImage* mimg, const iftImage *objmap, iftLabeledSet* seeds, float alpha, int beta){

    iftMaxflowGraph *graph = iftMaxflowGraphWithObjmap(mimg, objmap, seeds, alpha, beta);

    iftMaxflow(graph);

    iftImage* seg = iftMMaxflowSegment(graph, mimg);

    iftDestroyMaxflowGraph(&graph);

    return seg;
}


void iftProcessSinkOrphanOptPathMaxflow(iftMaxflowGraph *graph, iftGQueue *Q, iftMaxflowNode *p) {
    iftMaxflowArc *a_min = NULL;
    int d_min = IFT_INFINITY_INT;

    for (iftMaxflowArc *a0 = p->first; a0 != NULL; a0 = a0->next) {
        if (a0->r_cap != 0) {
            iftMaxflowNode *q = a0->head;
            iftMaxflowArc *arc = q->parent;

            if (q->is_sink == 1 && arc != NULL) {
                int d = 0;

                while (true) {
                    if (q->time_stamp == graph->time) {
                        d += q->dist;
                        break;
                    }
                    arc = q->parent;
                    d++;

                    if (arc == IFT_TERMINAL) {
                        q->time_stamp = graph->time;
                        q->dist = 1;
                        break;
                    }
                    if (arc == IFT_ORPHAN) {
                        d = IFT_INFINITY_INT;
                        break;
                    }
                    q = arc->head;
                }
                if (d < IFT_INFINITY_INT) {
                    /* q originates from the sink - done */

                    if (d < d_min) {
                        a_min = a0;
                        d_min = d;
                    }
                    /* set marks along the path */
                    for (q = a0->head; q->time_stamp != graph->time; q = q->parent->head) {
                        q->time_stamp = graph->time;
                        q->dist = d;
                        d--;
                    }
                }
            }
        }
    }
    p->parent = a_min;
    if (p->parent != NULL) {
        p->time_stamp = graph->time;
        p->dist = d_min + 1;
    } else {
        for (iftMaxflowArc *a0 = p->first; a0 != NULL; a0 = a0->next) {
            iftMaxflowNode *q = a0->head;
            iftMaxflowArc *arc = q->parent;
            if (q->is_sink == 1 && arc != NULL) {
                if(a0->r_cap != 0) {
                    iftTryInsertGQueue(Q, q->index);
                }
                if (arc != IFT_TERMINAL && arc != IFT_ORPHAN && arc->head == p) {
                        iftSetOrphanRearMaxflow(graph, q);
                }
            }
        }
    }
}


void iftProcessSourceOrphanOptPathMaxflow(iftMaxflowGraph *graph, iftGQueue *Q, iftMaxflowNode *p) {
    iftMaxflowArc *a_min = NULL;
    int d_min = IFT_INFINITY_INT;

    for (iftMaxflowArc *a0 = p->first; a0 != NULL; a0 = a0->next) {
        if (a0->sister->r_cap != 0) {
            iftMaxflowNode *q = a0->head;
            iftMaxflowArc *arc = q->parent;

            if (q->is_sink == 0 && arc != NULL) {
                int d = 0;
                while (true) {
                    if (q->time_stamp == graph->time) {
                        d += q->dist;
                        break;
                    }
                    arc = q->parent;
                    d++;
                    if (arc == IFT_TERMINAL) {
                        q->time_stamp = graph->time;
                        q->dist = 1;
                        break;
                    }
                    if (arc == IFT_ORPHAN) {
                        d = IFT_INFINITY_INT;
                        break;
                    }
                    q = arc->head;
                }

                if (d < IFT_INFINITY_INT) {
                    /* q originates from the sink - done */

                    if (d < d_min) {
                        a_min = a0;
                        d_min = d;
                    }
                    /* set marks along the path */
                    for (q = a0->head; q->time_stamp != graph->time; q = q->parent->head) {
                        q->time_stamp = graph->time;
                        q->dist = d;
                        d--;
                    }
                }
            }
        }
    }

    p->parent = a_min;
    if (p->parent != NULL) {
        p->time_stamp = graph->time;
        p->dist = d_min + 1;
    } else {
        for (iftMaxflowArc *a0 = p->first; a0 != NULL; a0 = a0->next) {
            iftMaxflowNode *q = a0->head;
            iftMaxflowArc *arc = q->parent;
            if (q->is_sink == 0 && arc != NULL) {
                if (a0->sister->r_cap != 0) {
                    iftTryInsertGQueue(Q, q->index);
                }
                if (arc != IFT_TERMINAL && arc != IFT_ORPHAN && arc->head == p) {
                    iftSetOrphanRearMaxflow(graph, q);
                }
            }
        }
    }
}


iftMaxflowNode *iftGetActiveNodeMaxflow(iftMaxflowGraph *graph, iftGQueue *Q) {

    while(!iftEmptyGQueue(Q)) {
        int index = iftRemoveGQueue(Q);
        if(graph->nodes[index].parent != NULL) {
            /* if parent different from NULL node is active */
            return &graph->nodes[index];
        }
    }

    return NULL;
}


float iftOptimumPathMaxflow(iftMaxflowGraph *graph, const iftImage* gradient) {

    /* init */
    int *weights_map = iftAlloc(graph->node_max, sizeof *weights_map);
    int max_grad = iftMaximumValue(gradient);
    iftGQueue *Q = iftCreateGQueue(max_grad, graph->node_max, weights_map);

    graph->time = 0;
    graph->orphan_first = NULL;
    graph->orphan_last = NULL;

    for (int i = 0; i < graph->node_last; i++) {
        iftMaxflowNode *p = &graph->nodes[i];

        p->index = i;
        p->next = NULL;
        p->time_stamp = graph->time;
        weights_map[p->index] = gradient->val[p->index];

        if (p->tr_cap > 0) {
            /* p is connected to the source */
            p->is_sink = 0;
            p->parent = IFT_TERMINAL;
            iftTryInsertGQueue(Q, p->index);
            p->dist = 1;
        } else if (p->tr_cap < 0) {
            /* p is connected to the sink */
            p->is_sink = 1;
            p->parent = IFT_TERMINAL;
            iftTryInsertGQueue(Q, p->index);
            p->dist = 1;
        } else {
            p->parent = NULL;
        }
    }
    /* end init */

    /* maxflow */
    iftMaxflowNode *current = NULL;
    iftMaxflowArc *a = NULL;

    int count = 0;

    while (!iftEmptyGQueue(Q)) {
        count++;

        iftMaxflowNode* p = current;

        if (p != NULL && p->parent == NULL) {
            iftRemoveGQueueElem(Q, p->index);
            p = NULL;
        }

        if (p == NULL) {
//            if(iftEmptyGQueue(Q)){
//                printf("Revisar o loop principal e alguma coisa não deve estar certa\n");
//                break;
//            }
            p = iftGetActiveNodeMaxflow(graph, Q);
        }

        /* growth */
        if (p->is_sink == 0) {
            /* grow source tree */
            for (a = p->first; a != NULL; a = a->next) {
                if (a->r_cap != 0) {
                    iftMaxflowNode *q = a->head;
                    if (q->parent == NULL) {
                        q->is_sink = 0;
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                        iftTryInsertGQueue(Q, q->index);
                    } else if (q->is_sink == 1) {
                        break;
                    } else if (q->time_stamp <= p->time_stamp && q->dist > p->dist) {
                        /* VERIFICAR NECESSIDADE DESSA HEURISTICA */
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                    }
                }
            }
        } else {
            /* grow sink tree */
            for (a = p->first; a != NULL; a = a->next) {
                if (a->sister->r_cap != 0) {
                    iftMaxflowNode *q = a->head;
                    if (q->parent == NULL) {
                        q->is_sink = 1;
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                        iftTryInsertGQueue(Q, q->index);
                    } else if (q->is_sink == 0) {
                        a = a->sister;
                        break;
                    } else if (q->time_stamp <= p->time_stamp && q->dist > p->dist) {
                        /* VERIFICAR NECESSIDADE DESSA HEURISTICA */
                        q->parent = a->sister;
                        q->time_stamp = p->time_stamp;
                        q->dist = p->dist + 1;
                    }
                }
            }
        }
        /* end growth */

        graph->time++;

        if (a != NULL) {
            iftTryInsertGQueue(Q, p->index);
            current = p;

            iftAugmentMaxflow(graph, a);

            /* adoption */
            iftMaxflowNodePtr *ptr, *ptr_next;
            while ( (ptr = graph->orphan_first) != NULL) {
                ptr_next = ptr->next;
                ptr->next = NULL;

                while ( (ptr = graph->orphan_first) != NULL) {
                    graph->orphan_first = ptr->next;
                    p = ptr->node;
                    iftFree(ptr);
                    if (graph->orphan_first == NULL) graph->orphan_last = NULL;
                    if (p->is_sink == 1) {
                        iftProcessSinkOrphanOptPathMaxflow(graph, Q, p);
                    } else {
                        iftProcessSourceOrphanOptPathMaxflow(graph, Q, p);
                    }
                }

                graph->orphan_first = ptr_next;
            }
            /* adoption end */
        } else {
            current = NULL;
        }
    }

    printf("count %d\n", count);

    iftFree(weights_map);
    iftDestroyGQueue(&Q);

    graph->iteration++;
    return graph->flow;
}


iftImage* iftOptimumPathGraphCut(const iftImage* gradient, iftLabeledSet* seeds, int beta){

    iftMaxflowGraph *graph = iftMaxflowGraphFromGradient(gradient, seeds, beta);

    iftOptimumPathMaxflow(graph, gradient);

    iftImage* seg = iftCreateImage(gradient->xsize, gradient->ysize, gradient->zsize);

    for(int p = 0; p < seg->n; p++){
        if (iftMaxflowNodeLabel(graph, p) != 0)
            seg->val[p] = 1;
    }

    iftDestroyMaxflowGraph(&graph);

    return seg;
}


double iftGraphCutBeta(iftMImage *mimg)
{
    double beta = 0;

    iftAdjRel *A = iftIs3DMImage(mimg) ? iftSpheric(1.0f) : iftCircular(1.0f);
    int count = 0;
    for (int p = 0; p < mimg->n; p++)
    {
        iftVoxel u = iftMGetVoxelCoord(mimg, p);
        for (int i = 1; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);
            if (iftValidVoxel(mimg, v)) {
                int q = iftGetVoxelIndex(mimg, v);
                beta += iftMImageSqDist(mimg, p, q);
                count++;
            }
        }
    }

    if (beta < IFT_EPSILON) beta = 0.0;
    else {
        beta = 1 / (2 * beta/count);
    }

    iftDestroyAdjRel(&A);
    return beta;
}
