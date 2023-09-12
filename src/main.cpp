// **************************************************************************
//  File       [main.cpp]
//  Author     [YU-AN, CHEN]
//  Synopsis   [The main program of PA2]
// **************************************************************************

#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
#include <climits>

using namespace std;

struct UD_vertex {

	int key = -101;
	int index;
	bool visited = false;
	UD_vertex* parent = nullptr;
	vector<UD_vertex*> adj;

};

struct UD_edge {

	UD_vertex* p;
	UD_vertex* s;

};

struct UD_cmp {

	bool operator() (UD_vertex* u, UD_vertex* v) {
		return u->key < v->key;
	}

};

class UD_Graph {

public:
	UD_Graph(int);   // constructor  v size 
	int W(int, int);
	void AppendEdges(int, int, int);

	vector<UD_vertex*> V;
	vector<UD_edge*> E;

private:
	vector<vector<int> > weight;

};


struct D_vertex {

	int start;
	int finish;
	int index;
	int original_index;
	int new_index;
	int house;

	char color = 'W';

	D_vertex* dfs_parent = nullptr;
	vector<D_vertex*> next;    // adj

};

struct D_edge {

	D_vertex* src;
	D_vertex* dest;

	int w;

};

class D_Graph {

public:
	D_Graph(int);    // constructor v size
	int W(int, int);
	void AppendEdges(int, int, int);

	vector<D_vertex*> V;
	vector<D_edge*> E;

private:
	vector<vector<int> > weight;

};

void DFS_SCC(D_Graph&);
int DFSVisit_SCC(D_Graph&, D_vertex*, int);

vector<D_Graph*> DFS_SCC_T(D_Graph&, D_Graph&);
int DFSVisit_SCC_T(D_Graph&, D_vertex*, int, D_Graph&, vector<D_vertex*>&);

void help_message() {
	cout << "usage: mps <input_file> <output_file>" << endl;
}

int main(int argc, char* argv[])
{

	//////////// read the input file /////////////
	
    fstream fin(argv[1]);
    fstream fout;
    fout.open(argv[2],ios::out);
	char type;
	int amount_vertices;
	int amount_edges;
	fin >> type >> amount_vertices >> amount_edges;

	// Undirected Graph
	if (type == 'u') {

		UD_Graph G(amount_vertices);
		int u_index;
		int v_index;
		int w_uv;

		for (int i = 0; i < amount_edges; i++) {               // record edge information
			fin >> u_index >> v_index >> w_uv;
			G.AppendEdges(u_index, v_index, w_uv);
		}

		UD_vertex* root = G.V[0];
		root->key = 101;

		priority_queue <UD_vertex*, vector<UD_vertex*>, UD_cmp> Q;

		for (int i = 0; i < amount_vertices; i++) {
			Q.push(G.V[i]);
		}

		// Use Prim Algorithm to find maximum spanning tree
		vector<bool> lazy_pq_checker(amount_vertices, false);
		while (!Q.empty()) {

			UD_vertex* u = Q.top();
			Q.pop();
			if (lazy_pq_checker[u->index]) continue;
			else lazy_pq_checker[u->index] = true;

			cout << u->index << " : {";
			for (int i = 0; i < amount_vertices; i++) {
				cout << "v" << i << " : " << G.V[i]->key << ", ";
			}
			cout << "end }" << endl;

			u->visited = true;

			int adj_size = u->adj.size();
			for (int i = 0; i < adj_size; i++) {

				UD_vertex* v = u->adj[i];
				if (!v->visited && G.W(u->index, v->index) > v->key) {
					v->parent = u;
					v->key = G.W(u->index, v->index);
					Q.push(v);    // the adjusted key would be popped out first anyway 

				}

			}

		}

		// Record edges excluded from MSP and the total weight
		vector<UD_edge*> solution;
		int total_weight = 0;
		UD_edge* e;
		for (int i = 0; i < amount_edges; i++) {

			e = G.E[i];
			if (!(e->p->parent == e->s || e->s->parent == e->p)) {

				total_weight += G.W(e->p->index, e->s->index);
				solution.push_back(e);

			}

			else if (G.W(e->p->index, e->s->index) < 0) {

				total_weight += G.W(e->p->index, e->s->index);
				solution.push_back(e);

			}

		}

		fout << total_weight << endl;
		for (int i = 0; i < solution.size(); i++) {

			e = solution[i];
			int src = e->p->index;
			int dest = e->s->index;
			fout << src << " " << dest << " " << G.W(src, dest) << endl;

		}

	}

	// Directed Graph
	else {

		vector<D_edge*> solution;
		int weight = 0;

		D_Graph G(amount_vertices);
		D_Graph GT(amount_vertices);

		int u_index;
		int v_index;
		int w_uv;
		D_edge* neg_e;

		for (int i = 0; i < amount_edges; i++) {               // record edge information

			fin >> u_index >> v_index >> w_uv;
            
			if (w_uv >= 0) {

				G.AppendEdges(u_index, v_index, w_uv);
				GT.AppendEdges(v_index, u_index, w_uv);
			}

			else {

				neg_e = new D_edge;
				neg_e->src = G.V[u_index];
				neg_e->dest = G.V[v_index];
				neg_e->w = w_uv;

				solution.push_back(neg_e);
				weight += w_uv;

			}

		}

		// Arrange SCC vertex
		DFS_SCC(G);
		vector<D_Graph*> SCC_Set = DFS_SCC_T(GT, G);


		// Arrange SCC edges
		D_edge* e;
		for (int i = 0; i < amount_edges; i++) {
		
			e = G.E[i];
			if (e->src->house != e->dest->house) continue;
		    u_index = e->src->new_index;
			v_index = e->dest->new_index;
			(*SCC_Set[e->src->house]).AppendEdges(u_index, v_index, e->w);

		}

		// Solve in each SCC
		D_Graph* SCC;
		int v_size;
		int min_weight;
		int temp_weight;
		int adj_size;
		D_vertex* u;
		D_vertex* v;
		vector<int> perm, min_perm;
		
		for (int i = 0; i < SCC_Set.size(); i++) {
		
			SCC = SCC_Set[i];

			cout << "Set" << i << " = {";
			for (int i = 0; i < SCC->V.size(); i++) {
				cout << SCC->V[i]->original_index << ", ";
			}
			cout << "end}" << endl;

			v_size = SCC->V.size();
			min_weight = INT_MAX;
			perm.clear();
			min_perm.clear();
			
			for (int i = 0; i < v_size; i++) {

				perm.push_back(i);              // perm[index] = order
				min_perm.push_back(i);

			}

			do {
			
				temp_weight = 0;
				for (int i = 0; i < v_size; i++) {
				
					u = SCC->V[i];
					adj_size = u->next.size();
					for (int j = 0; j < adj_size; j++) {
					
						v = u->next[j];
						if (perm[u->index] > perm[v->index]) {
							temp_weight += SCC->W(u->index, v->index);
						}
					
					}

				}

				if (temp_weight < min_weight) {
				
					min_weight = temp_weight;
					min_perm = perm;

				}

			} while (next_permutation(perm.begin(), perm.end()));

			weight += min_weight;

			for (int i = 0; i < v_size; i++) {

				u = SCC->V[i];
				adj_size = u->next.size();


				for (int j = 0; j < adj_size; j++) {

					v = u->next[j];
					if (min_perm[u->index] > min_perm[v->index]) {

						e = new D_edge;
						e->src = G.V[u->original_index];
						e->dest = G.V[v->original_index];
						e->w = (*SCC).W(u->index, v->index);

						solution.push_back(e);

					}

				}

			}

			cout << "minPerm" << i << " = {";
			for (int i = 0; i < min_perm.size(); i++) {
				cout << min_perm[i] << ", ";
			}
			cout << "end}" << endl;

		}

		fout << weight << endl;
		for (int i = 0; i < solution.size(); i++) {

			e = solution[i];
			int src = e->src->index;
			int dest = e->dest->index;
			fout << src << " " << dest << " " << G.W(src, dest) << endl;

		}

	}

	fin.close();
	fout.close();
	return 0;

}

// **************************************************************************

UD_Graph::UD_Graph(int v_size) {

	for (int i = 0; i < v_size; i++) {
		UD_vertex* u = new UD_vertex;
		u->index = i;
		V.push_back(u);
	}

	vector<vector<int> > weight_table(v_size, vector<int>(v_size));
	weight = weight_table;

}

int UD_Graph::W(int u, int v) {
	return weight[u][v];
}

void UD_Graph::AppendEdges(int u, int v, int w) {

	// weight
	weight[u][v] = w;
	weight[v][u] = w;

	// adjacency
	UD_vertex* vertex_u = V[u];
	UD_vertex* vertex_v = V[v];
	vertex_u->adj.push_back(vertex_v);
	vertex_v->adj.push_back(vertex_u);

	// edge
	UD_edge* e = new UD_edge;
	e->p = vertex_u;
	e->s = vertex_v;
	E.push_back(e);

}

//*****************************************************************************

D_Graph::D_Graph(int v_size) {
	for (int i = 0; i < v_size; i++) {
		D_vertex* u = new D_vertex;
		u->index = i;
		V.push_back(u);

	}

	vector<vector<int> > weight_table(v_size, vector<int>(v_size));
	weight = weight_table;

}

void D_Graph::AppendEdges(int u, int v, int w) {

	// weight
	weight[u][v] = w;

	// adjacency
	D_vertex* vertex_u = V[u];
	D_vertex* vertex_v = V[v];
	vertex_u->next.push_back(vertex_v);

	// edge
	D_edge* e = new D_edge;
	e->src = vertex_u;
	e->dest = vertex_v;
	e->w = w;
	E.push_back(e);

}

int D_Graph::W(int u, int v) {
	return weight[u][v];
}
//*****************************************************************************
void DFS_SCC(D_Graph& G) {

	int v_size = G.V.size();
	D_vertex* u;
	for (int i = 0; i < v_size; i++) {

		u = G.V[i];
		u->color = 'W';
		u->dfs_parent = nullptr;

	}

	int time = 0;
	for (int i = 0; i < v_size; i++) {

		u = G.V[i];
		if (u->color == 'W')
			time = DFSVisit_SCC(G, u, time);

	}

}

int DFSVisit_SCC(D_Graph& G, D_vertex* u, int t) {

	int time = t + 1;
	u->start = time;
	u->color = 'G';

	int adj_size = u->next.size();
	D_vertex* v;
	for (int i = 0; i < adj_size; i++) {

		v = u->next[i];
		if (v->color == 'W') {
		
			v->dfs_parent = u;
			time = DFSVisit_SCC(G, v, time);

		}

	}

	u->color = 'B';
	time = ++time;
	u->finish = time;

	return time;
	
}

vector<D_Graph*> DFS_SCC_T(D_Graph& GT, D_Graph& G) {

	int v_size = GT.V.size();
	D_vertex* u;
	D_vertex* original_u;
	vector<D_vertex*> Increasing_Finishing_Order(v_size * 2, nullptr);

	for (int i = 0; i < v_size; i++) {

		u = GT.V[i];
		u->color = 'W';
		u->dfs_parent = nullptr;

		original_u = G.V[u->index];
		Increasing_Finishing_Order[original_u->finish - 1] = u;

	}


	vector<D_vertex*> v_queue;
	vector<D_edge*> e_queue;

	vector<D_Graph*> SCC_Set;
	D_Graph* SCC;

	int time = 0;
	for (int i = v_size * 2 - 1; i >= 0; i--) {

		u = Increasing_Finishing_Order[i];

		if (u == nullptr) continue;

		if (u->color == 'W') {

			v_queue.clear();
			e_queue.clear();
			v_queue.push_back(u);

			time = DFSVisit_SCC_T(GT, u, time, G, v_queue);
			SCC = new D_Graph(v_queue.size());

			for (int i = 0; i < v_queue.size(); i++) {

				v_queue[i]->new_index = i;
				G.V[v_queue[i]->index]->new_index = i;

				SCC->V[i]->original_index = v_queue[i]->index;
				SCC->V[i]->house = SCC_Set.size();

				v_queue[i]->house = SCC_Set.size();
				G.V[v_queue[i]->index]->house = SCC_Set.size();

			}

			SCC_Set.push_back(SCC);

		}

	}

	return SCC_Set;

}

int DFSVisit_SCC_T(D_Graph& GT, D_vertex* u, int t, D_Graph& G, vector<D_vertex*>& v_queue) {

	int time = t + 1;
	u->start = time;
	u->color = 'G';

	int adj_size = u->next.size();
	D_vertex* v;
	for (int i = 0; i < adj_size; i++) {

		v = u->next[i];
		if (v->color == 'W') {
			
			v_queue.push_back(v);
			v->dfs_parent = u;
			time = DFSVisit_SCC_T(GT, v, time, G, v_queue);

		}

	}

	u->color = 'B';
	time = ++time;
	u->finish = time;

	return time;

}
