//Hu-Shing method of Matrix Chain Optimization
//Implemented in C by Kathiravan - CS22B2052
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdbool.h>

#define ll long long
#define testsize 5000000  //input array size. mind the size, since it must be malloc-ed
                          // for testsize = 5,000,000, approx 5.6 GB of RAM is used. 
#define DEBUG false


//gcc hu-shing-3.c -pg -o hu-shing-3 && ./hu-shing-3 && gprof hu-shing-3 gmon.out > hsc3_gprof_report.txt && cat hsc3_gprof_report.txt

typedef struct {
    int id, u, v, low;
    ll base, mul, num, den;
} HArc;

typedef struct {
    int size;
    int capacity;
    HArc* data;
} PriorityQueue;

typedef struct {
    int size;
    int capacity;
    int* data;
} IntVector;

typedef struct {
    int size;
    int capacity;
    HArc* data;
} HArcVector;

int n;
ll *w, *CP;
HArc *h;
int *sub, *qid;
int n_arcs;
int n_pqs;
IntVector *child;
PriorityQueue *pq;
HArcVector *con;


void initPriorityQueue(PriorityQueue* pq) {
    pq->size = 0;
    pq->capacity = 10;
    pq->data = (HArc*)malloc(pq->capacity * sizeof(HArc));
}

void pushPriorityQueue(PriorityQueue* pq, HArc arc) {
    if (pq->size == pq->capacity) {
        pq->capacity *= 2;
        pq->data = (HArc*)realloc(pq->data, pq->capacity * sizeof(HArc));
    }
    pq->data[pq->size++] = arc;
    //if(arc.den == 0) arc.den = 1; //added to solve FPE
    // Simple insertion sort for priority queue (can be optimized)
    for (int i = pq->size - 1; i > 0; i--) {
        //if (pq->data[i].den == 0 || pq->data[i-1].den == 0) continue; //also for FPE
        if (pq->data[i].num / pq->data[i].den < pq->data[i - 1].num / pq->data[i - 1].den) {
            HArc temp = pq->data[i];
            pq->data[i] = pq->data[i - 1];
            pq->data[i - 1] = temp;
        } else {
            break;
        }
    }
}

HArc topPriorityQueue(PriorityQueue* pq) {
    return pq->data[0];
}

void popPriorityQueue(PriorityQueue* pq) {
    for (int i = 1; i < pq->size; i++) {
        pq->data[i - 1] = pq->data[i];
    }
    pq->size--;
}

int emptyPriorityQueue(PriorityQueue* pq) {
    return pq->size == 0;
}

void initIntVector(IntVector* vec) {
    vec->size = 0;
    vec->capacity = 10;
    vec->data = (int*)malloc(vec->capacity * sizeof(int));
}

void pushIntVector(IntVector* vec, int value) {
    if (vec->size == vec->capacity) {
        vec->capacity *= 2;
        vec->data = (int*)realloc(vec->data, vec->capacity * sizeof(int));
    }
    vec->data[vec->size++] = value;
}

void popIntVector(IntVector* vec) {
    if (vec->size > 0) {
        vec->size--;
    }
}

void initHArcVector(HArcVector* vec) {
    vec->size = 0;
    vec->capacity = 10;
    vec->data = (HArc*)malloc(vec->capacity * sizeof(HArc));
}

void pushHArcVector(HArcVector* vec, HArc arc) {
    if (vec->size == vec->capacity) {
        vec->capacity *= 2;
        vec->data = (HArc*)realloc(vec->data, vec->capacity * sizeof(HArc));
    }
    vec->data[vec->size++] = arc;
}

void popHArcVector(HArcVector* vec) {
    vec->size--;
}

void new_arc(int u, int v) {
    assert(u <= v);
    n_arcs++;
    h[n_arcs].id = n_arcs;
    h[n_arcs].u = u;
    h[n_arcs].v = v;
    h[n_arcs].low = w[u] < w[v] ? u : v;
    h[n_arcs].mul = (ll)w[u] * w[v];
    h[n_arcs].base = CP[v] - CP[u] - h[n_arcs].mul;
}

void build_tree(int lst[][2], int lst_size) {
    IntVector stk;
    initIntVector(&stk);
    new_arc(1, n + 1); // h[1] is root
    for (int i = 0; i < lst_size; i++) {
        new_arc(lst[i][0], lst[i][1]);
        while (stk.size > 0 && h[n_arcs].u <= h[stk.data[stk.size - 1]].u && h[n_arcs].v >= h[stk.data[stk.size - 1]].v) {
            pushIntVector(&child[n_arcs], stk.data[stk.size - 1]);
            popIntVector(&stk); // Fixed: popIntVector is now defined
        }
        pushIntVector(&stk, n_arcs);
    }
    while (stk.size > 0) {
        pushIntVector(&child[1], stk.data[stk.size - 1]);
        popIntVector(&stk); // Fixed: popIntVector is now defined
    }
    free(stk.data);
}

void one_sweep(int lst[][2], int* lst_size, int size) {
    IntVector stk;
    initIntVector(&stk);
    int (*tmp)[2] = (int (*)[2])malloc(size * sizeof(int[2]));
    int tmp_size = 0;
    for (int i = 1; i <= n; i++) {
        while (stk.size >= 2 && w[stk.data[stk.size - 1]] > w[i]) {
            tmp[tmp_size][0] = stk.data[stk.size - 2];
            tmp[tmp_size][1] = i;
            tmp_size++;
            popIntVector(&stk);
        }
        pushIntVector(&stk, i);
    }
    while (stk.size >= 4) {
        int Vt_1 = stk.data[stk.size - 2];
        tmp[tmp_size][0] = 1;
        tmp[tmp_size][1] = Vt_1;
        tmp_size++;
        popIntVector(&stk);
    }
    *lst_size = 0;
    for (int i = 0; i < tmp_size; i++) {
        if (tmp[i][0] == 1 || tmp[i][1] == 1) continue;
        lst[*lst_size][0] = tmp[i][0];
        lst[*lst_size][1] = tmp[i][1];
        (*lst_size)++;
    }
    free(tmp);
    free(stk.data);
}

void prepare(int size) {
    int V1 = 1;
    for (int i = 2; i <= n; i++) {
        if (w[i] < w[V1]) V1 = i;
    }
    // Rotate w so that w[1] is the smallest
    ll *temp = (ll*)malloc(size * sizeof(ll));
    if(temp == NULL) printf("Malloc Failed for Temp in void prepare\n");
    for (int i = 1; i <= n; i++) {
        temp[i] = w[i];
    }
    for (int i = V1; i <= n; i++) {
        w[i - V1 + 1] = temp[i];
    }
    for (int i = 1; i < V1; i++) {
        w[n - V1 + 1 + i] = temp[i];
    }
    w[n + 1] = w[1];
    for (int i = 1; i <= n + 1; i++) {
        CP[i] = (ll)w[i] * w[i - 1];
        CP[i] += CP[i - 1];
    }
    free(temp);
}

ll get_mn_mul(int node) {
    if (node == 1) return (ll)w[1] * w[2] + (ll)w[1] * w[n];
    HArc* cur = &h[node];
    if (cur->u == cur->low) {
        if (con[cur->u].size == 0 || cur->u > con[cur->u].data[con[cur->u].size - 1].u || cur->v < con[cur->u].data[con[cur->u].size - 1].v) {
            return (ll)w[cur->u] * w[cur->u + 1];
        } else {
            return con[cur->u].data[con[cur->u].size - 1].mul;
        }
    } else {
        if (con[cur->v].size == 0 || cur->u > con[cur->v].data[con[cur->v].size - 1].u || cur->v < con[cur->v].data[con[cur->v].size - 1].v) {
            return (ll)w[cur->v] * w[cur->v - 1];
        } else {
            return con[cur->v].data[con[cur->v].size - 1].mul;
        }
    }
    assert(0); // never happens
    return 0;
}

void add_arc(int cur_node, HArc* arc) {
    if(arc->den == 0) arc->den = 1;  //added for FPE, might tamper with output results
    pushPriorityQueue(&pq[qid[cur_node]], *arc);
    pushHArcVector(&con[arc->u], *arc);
    pushHArcVector(&con[arc->v], *arc);
}

void remove_arc(int cur_node) {
    HArc hm = topPriorityQueue(&pq[qid[cur_node]]);
    popHArcVector(&con[hm.u]);
    popHArcVector(&con[hm.v]);
    popPriorityQueue(&pq[qid[cur_node]]);
}

void merge_pq(int node) {
    int max_child = -1;
    for (int i = 0; i < child[node].size; i++) {
        if (max_child == -1 || sub[max_child] < sub[child[node].data[i]]) {
            max_child = child[node].data[i];
        }
    }
    qid[node] = qid[max_child];
    PriorityQueue* cur_pq = &pq[qid[node]];
    for (int i = 0; i < child[node].size; i++) {
        if (child[node].data[i] == max_child) continue;
        PriorityQueue* child_pq = &pq[qid[child[node].data[i]]];
        while (!emptyPriorityQueue(child_pq)) {
            pushPriorityQueue(cur_pq, topPriorityQueue(child_pq));
            popPriorityQueue(child_pq);
        }
    }
}

void dfs(int node) {
    HArc* cur = &h[node];
    sub[node] = 1;
    if (child[node].size == 0) {
        qid[node] = ++n_pqs;
        cur->den = cur->base;
        cur->num = (ll)w[cur->low] * (cur->den + cur->mul - get_mn_mul(node));
        add_arc(node, cur);
        return;
    }
    cur->den = cur->base;
    for (int i = 0; i < child[node].size; i++) {
        dfs(child[node].data[i]);
        sub[node] += sub[child[node].data[i]];
        cur->den -= h[child[node].data[i]].base;
    }
    cur->num = (ll)w[cur->low] * (cur->den + cur->mul - get_mn_mul(node));
    merge_pq(node);
    PriorityQueue* cur_pq = &pq[qid[node]];
    while (!emptyPriorityQueue(cur_pq) && topPriorityQueue(cur_pq).num / topPriorityQueue(cur_pq).den >= w[cur->low]) {
        HArc hm = topPriorityQueue(cur_pq);
        cur->den += hm.den;
        remove_arc(node);
        cur->num = (ll)w[cur->low] * (cur->den + cur->mul - get_mn_mul(node));
    }
    while (!emptyPriorityQueue(cur_pq) && cur->num / cur->den <= topPriorityQueue(cur_pq).num / topPriorityQueue(cur_pq).den) {
        HArc hm = topPriorityQueue(cur_pq);
        cur->den += hm.den;
        remove_arc(node);
        cur->num += hm.num;
    }
    add_arc(node, cur);
}

ll getans() {
    dfs(1);
    ll ans = 0;
    PriorityQueue* cur_pq = &pq[qid[1]];
    while (!emptyPriorityQueue(cur_pq)) {
        ans += topPriorityQueue(cur_pq).num;
        popPriorityQueue(cur_pq);
    }
    return ans;
}

void init(int size) {
    w = (ll*)malloc((size + 2) * sizeof(ll));  // Extra space for w[n+1]
    CP = (ll*)malloc((size + 2) * sizeof(ll));
    h = (HArc*)malloc((size + 2) * sizeof(HArc));
    sub = (int*)malloc((size + 2) * sizeof(int));
    qid = (int*)malloc((size + 2) * sizeof(int));

    if (!w || !CP || !h || !sub || !qid) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    memset(w, 0, (size + 2) * sizeof(ll));
    memset(CP, 0, (size + 2) * sizeof(ll));
    memset(sub, 0, (size + 2) * sizeof(int));
    memset(qid, 0, (size + 2) * sizeof(int));

    n_arcs = n_pqs = 0;
    
    child = (IntVector*)malloc((size + 2) * sizeof(IntVector));
    pq = (PriorityQueue*)malloc((size + 2) * sizeof(PriorityQueue));
    con = (HArcVector*)malloc((size + 2) * sizeof(HArcVector));

    if (!child || !pq || !con) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < size + 2; i++) {
        initIntVector(&child[i]);
        initPriorityQueue(&pq[i]);
        initHArcVector(&con[i]);
    }

}


void input(ll arr[], int size) {
    n = size;
    for (int i = 1; i <= n; i++) {
        w[i] = arr[i - 1];
    }
}

void cleanup(int size) {
    free(w);
    free(CP);
    free(h);
    free(sub);
    free(qid);

    for (int i = 0; i < size + 2; i++) {
        free(child[i].data);
        free(pq[i].data);
        free(con[i].data);
    }

    free(child);
    free(pq);
    free(con);
}

ll solve(ll arr[], int size) {
    if (size < 2) return 0;
    if (size == 2) return arr[0] * arr[1];

    init(size);
    input(arr, size);
    prepare(size);

    int (*lst)[2] = (int (*)[2])malloc(size * sizeof(int[2]));
    int lst_size;

    one_sweep(lst, &lst_size, size);
    build_tree(lst, lst_size);
    free(lst);

    ll ans = getans();
    cleanup(size);

    return ans;
}


int main() {
    //ll arr[] = {30, 35, 15, 5, 10, 20};
    ll *arr = (ll*) calloc(testsize, sizeof(ll));
    if(!arr) {
      printf("Calloc for Input Array Failed\n");
      return -1;
    }
    for(int i=0;i<testsize;i++) {
      arr[i] = rand()% (100 + 1 - 5) + 5;
    }
    int size = testsize;
    if(DEBUG) for(int i=0;i<testsize;i++) printf("%lld, ", arr[i]);
    long long sol = solve(arr, size);
    if(DEBUG) printf("\nMin Ops: %lld\n", sol);
    printf("Finished Execution\n");
    free(arr);
    return 0;
}
