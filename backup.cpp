#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

//g++ parallel-mcm.cpp -fopenmp -pg -o parallel-mcm && gprof parallel-mcm gmon.out > gprof_report.txt && cat gprof_report.txt
//g++ parallel-mcm.cpp -fopenmp -o parallel-mcm && ./parallel-mcm

#define FILENAME "inp.txt"
#define NUM_TRIALS 1

void read_data(vector<long long>& arr) {
    std::ifstream file(FILENAME);
    long long x;
    while (file >> x) { 
        arr.push_back(x);
    }
}

namespace HuShing {
	using ll=long long;
	using pii=pair<int,int>;
	/* About Polygon */
	int n;
	vector<ll> w, CP; 
  
	/* About H-arcs */
	struct HArc {
		int id; // calc in new_arc
		int u,v; // "
		int low; // "
		ll base; // "
		ll mul; // "
		ll num,den; // calc in dfs
		inline bool contains(HArc& b)const{
			return u<=b.u&&b.v<=v;
		}
		inline ll get_support()const{
			assert(den!=0);
			return num/den;
		}
		bool operator<(const HArc& b)const{
			return get_support()<b.get_support();
		}
		bool operator<=(const HArc& b)const{
			return get_support()<=b.get_support();
		}
		bool operator==(const HArc& b)const{
			return get_support()==b.get_support();
		}
	}; 

        vector<HArc> h;
        
	int n_arcs;
	
	vector<int> sub; // calc in dfs
	vector<vector<int>> child; // calc in build_tree
	
	vector<int> qid; // calc in dfs
	
	int n_pqs;
	
	vector<priority_queue<HArc>> pq;
        vector<vector<HArc>> con;
        
        //

	void new_arc(int u,int v) {
		assert(u<=v);
		++n_arcs;
		h[n_arcs].id=n_arcs;
		h[n_arcs].u=u;
		h[n_arcs].v=v;
		h[n_arcs].low=w[u]<w[v]?u:v;
		h[n_arcs].mul=(ll)w[u]*w[v];
		h[n_arcs].base=CP[v]-CP[u]-h[n_arcs].mul;
	}
	void build_tree(vector<pii>& lst) {
		vector<int> stk;
		new_arc(1,n+1); // h[1] is root
		for(auto& it:lst) {
			new_arc(it.first,it.second);
			while(!stk.empty()&&
				h[n_arcs].contains(h[stk.back()])) {
				child[n_arcs].push_back(stk.back());
				stk.pop_back();
			}
			stk.push_back(n_arcs);
		}
		while(!stk.empty()) {
			child[1].push_back(stk.back());
			stk.pop_back();
		}
	}

	void one_sweep() {
		vector<int> stk;
		vector<pii> tmp,lst;
		for(int i=1;i<=n;i++) {
			while(stk.size()>=2 && w[stk.back()]>w[i]) {
				tmp.push_back({stk[stk.size()-2],i});
				stk.pop_back();
			}
			stk.push_back(i);
		}
		while(stk.size()>=4) {
			int Vt_1=stk[stk.size()-2];
			tmp.push_back({1,Vt_1});
			stk.pop_back();
		}
		for(auto& it:tmp) {
			if(it.first==1||it.second==1) continue;
			lst.push_back(it);
		}
		build_tree(lst);
	}
	void prepare() {
		int V1 = min_element(w.begin() + 1, w.begin() + n + 1) - w.begin();
                rotate(w.begin() + 1, w.begin() + V1, w.begin() + n + 1);

		w[n+1]=w[1];
		for(int i=1;i<=n+1;i++) {
			CP[i]=(ll)w[i]*w[i-1];
			CP[i]+=CP[i-1];
		}
	}
	ll get_mn_mul(int node) {
		if(node==1) return (ll)w[1]*w[2]+(ll)w[1]*w[n];
		HArc& cur=h[node];
		if(cur.u==cur.low) {
			if(con[cur.u].empty()||!cur.contains(con[cur.u].back())) {
				return (ll)w[cur.u]*w[cur.u+1];
			} else return con[cur.u].back().mul;
		} else {
			if(con[cur.v].empty()||!cur.contains(con[cur.v].back())) {
				return (ll)w[cur.v]*w[cur.v-1];
			} else return con[cur.v].back().mul;
		}
		assert(0); // never happens
		return 0;
	}
	inline void add_arc(int cur_node,HArc& arc) {
		pq[qid[cur_node]].push(arc);
		con[arc.u].push_back(arc);
		con[arc.v].push_back(arc);
	}
	inline void remove_arc(int cur_node) {
		const HArc& hm=pq[qid[cur_node]].top();
		con[hm.u].pop_back();
		con[hm.v].pop_back();
		pq[qid[cur_node]].pop();
	}
	void merge_pq(int node) {
		int max_child=-1;
		for(auto& it:child[node])
			if(max_child==-1||sub[max_child]<sub[it])
				max_child=it;
		qid[node]=qid[max_child];
		auto& cur_pq=pq[qid[node]];
		for(auto& it:child[node]) {
			if(it==max_child) continue;
			auto& child_pq=pq[qid[it]];
			while(!child_pq.empty()) {
				cur_pq.push(child_pq.top());
				child_pq.pop();
			}
		}
	}
	void dfs(int node) {
		HArc& cur=h[node];
		sub[node]=1;
		if(child[node].empty()) {
			qid[node]=++n_pqs;
			cur.den=cur.base;
			cur.num=(ll)w[cur.low]*(cur.den+cur.mul-get_mn_mul(node));
			add_arc(node,cur);
			return;
		}
		
		cur.den=cur.base;
		for(auto& it:child[node]) {
			dfs(it);
			sub[node]+=sub[it];
			cur.den-=h[it].base;
		}
		cur.num=(ll)w[cur.low]*(cur.den+cur.mul-get_mn_mul(node));
		merge_pq(node);
		auto& cur_pq=pq[qid[node]];
		while(!cur_pq.empty()&&cur_pq.top().get_support()>=w[cur.low]) {
			auto hm=cur_pq.top();
			cur.den+=hm.den;
			remove_arc(node); // this must be done before calculating cur.num!
			cur.num=(ll)w[cur.low]*(cur.den+cur.mul-get_mn_mul(node));
		}
		while(!cur_pq.empty()&&cur<=cur_pq.top()) {
			auto hm=cur_pq.top();
			cur.den+=hm.den;
			remove_arc(node);
			cur.num+=hm.num;
		}
		add_arc(node,cur);
	}
	ll getans() {
		dfs(1);
		ll ans=0;
		auto& cur_pq=pq[qid[1]];
		while(!cur_pq.empty()) {
			ans+=cur_pq.top().num;
			cur_pq.pop();
		}
		return ans;
	}
	void init(int n, int num_threads) {
		int i;
		w.resize(n+2);
		CP.resize(n+2);
		h.resize(n+2);
		
		sub.resize(n+2);
		child.assign(n+2, vector<int>());
	        
	        qid.resize(n+2);
	        
		pq.assign(n + 2, priority_queue<HArc>());
                con.assign(n + 2, vector<HArc>());
                
		n_arcs=n_pqs=0;
		
		omp_set_num_threads(num_threads);
	}
	void input(const vector<ll>& arr) {
		n=arr.size();
		for(int i=1;i<=n;i++) w[i]=arr[i-1];
	}
	ll solve(const vector<ll>& arr, int num_threads) {
		if(arr.size()<2) return 0;
		if(arr.size()==2) return arr[0]*arr[1];
		init(arr.size(), num_threads);
		input(arr);
		prepare();
		one_sweep();
		ll sol = getans();
		printf("%lld\n",sol);
		return sol;
	}
}

void benchmark(long long int(*func)(const std::vector<long long>&, int), std::vector<long long>& A, std::string method) {
    //std::vector<int> thread_counts = {1, 2, 4, 6, 8, 10, 12, 16, 20, 32, 64};
    std::vector<int> thread_counts = {1, 4, 8, 1};
    int num_tests = thread_counts.size();
    std::vector<double> times(num_tests);
    double base_time;
    
    std::cout<<method<<"\n";
    printf("Threads\tTime (s)\n");
    for (int i = 0; i < num_tests; i++) {
        int num_threads = thread_counts[i];
        double total_time = 0.0;
        
        for (int j = 0; j < NUM_TRIALS; j++) {
            double start = omp_get_wtime();
            
            //Test segment code
            func(A, num_threads);
            
            //End of Test segment
            double end = omp_get_wtime();
            total_time += (end - start);
        }
        times[i] = total_time / NUM_TRIALS;
        if (i == 0) base_time = times[i]; // Time with 1 thread
        printf("%d\t%.6f\n", num_threads, times[i]);
    }

    // Speedup calculation
    std::vector<double> parallel_fracs(num_tests);
    std::cout<<"\nSpeedup vs Processors:\n";
    std::cout<<"Threads\tSpeedup\n";
    for (int i = 0; i < num_tests; i++) {
        // Estimate Parallelization Fraction (Amdahl's Law)
        if(i>0) parallel_fracs[i] = ((times[i] / base_time) - 1.0) / ((1.0 / thread_counts[i]) - 1.0);
        printf("%d\t%.6f\t", thread_counts[i], base_time / times[i]);
        //printf("%.6f, ", base_time / times[i]);
        printf("Parallelization Fraction (Amdahl's Law): %.3f\n", parallel_fracs[i]);
    }
    
    printf("Estimated Avg Parallelization Fraction (Amdahl's Law): %.3f\n", std::accumulate(parallel_fracs.begin(), parallel_fracs.end(), 0.0)/(num_tests - 1));

}

int main() {
        vector<long long> arr;
        read_data(arr);
	
	benchmark(HuShing::solve, arr, "Parallel MCM");
	return 0;
}
