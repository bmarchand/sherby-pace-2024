#include <vector>
#include <list>
#include <deque>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <iostream>
#include <fstream>

#include "Util.h"

#include "EvalMaxSAT.h"

#define i32 int32_t
#define uset unordered_set
#define umap unordered_map

using namespace std;

struct pair_hasher
    {
        size_t operator()(const pair<int , int> &x) const{
        return x.first ^ x.second;
    }

};


class Digraph
{
public:

    Digraph()
    {
        
    }
    
    ~Digraph()
    {
        
    }

    vector< int > edge_weights;    //eg edge_weights[i] = weight of it-th edge
    vector< pair< int, int > > all_edges;
    umap< pair<int, int>, int, pair_hasher > edge_indexes;  //edge_indexes[ (i, j) ] = the index of that edge in all_edges

    int nb_nodes;
    vector< uset<int> > ins;    
    vector< uset<int> > outs;
    
    int last_cost = 0;
    
    uset< pair<int, int>, pair_hasher > deleted_edges;
    
    long max_weight = 0;
    
    
    //note: deleted edges are not written
    void write_graph(string filename)
    {
        //todo: output streaming would be better
        string strout = "";
        for (int i = 0; i < nb_nodes; ++i)
        {
            for (int j : outs[i])
            {
                int e_index = edge_indexes[ make_pair(i, j) ];
                strout += Util::ToString(i + 1) + " " + Util::ToString(j + 1) + " " + Util::ToString( edge_weights[ e_index ] );
                strout += "\n";
            }
        }
        
        Util::WriteFileContent(filename, strout);
    }
    
    
    
    void write_solution_info(vector<int> &topo_sort, string filename)
    {
        string strout = "";
        
        strout += "cost=" + Util::ToString(last_cost) + "\n";
        
        strout += "topo_sort=";
        for (int v : topo_sort)
        {
            strout += Util::ToString(v) + " ";
        }
        
        Util::WriteFileContent(filename, strout);
        
    }
    
    
    void to_wcnf( string filename, vector<deque<int>> &cycles )
    {
        ofstream f;
        f.open (filename);
        
        //format = p wcnf nbvars nbclauses maxw
        f << "p wcnf "<< all_edges.size() << " " << all_edges.size() + cycles.size() << " " << max_weight << "\n";
        
        //var constraints
        for (int i = 0; i < all_edges.size(); ++i)
        {
            
            if (edge_weights[i] > 10000)    //TODO: cleaner thing for high weight
            {
                f << max_weight << " -" << i + 1 << " 0\n";    
            }
            else 
            {
                f << edge_weights[i] << " -" << i + 1 << " 0\n";
            }
        }
        
        
        //cycle constraints 
        for (deque<int> cyc : cycles)
        {
            f << max_weight << " ";
            for (int trindex = 0; trindex < cyc.size(); ++trindex)
            {
                pair<int, int> arc = make_pair(cyc[trindex], cyc[ (trindex + 1) % cyc.size() ]);
                int e_index = edge_indexes[ arc ];
                
                f << e_index + 1 << " ";
            }
            f << "0\n";
            
        }
        
        f.close();
        
    }
    
    void read_graph_lines(string filelines) {

        vector<string> unfilteredLines = Util::Split(filelines, "\n");

        vector<string> lines;
        for (int i = 0; i < unfilteredLines.size(); i++)
        {
            if (unfilteredLines[i] != "")
                lines.push_back(unfilteredLines[i]);
        }
        
        //first line should have 
        //nbnodes nbarcs
        vector<string> px = Util::Split(lines[0], " ", false);
        
        nb_nodes = Util::ToInt(px[0]);
        
        cout<<"nbnodes="<<nb_nodes<<endl;
        
        //we assume node indices are indexed at 1 (ie first node has index 1, not 0)
        ins.resize(nb_nodes);
        outs.resize(nb_nodes);
        
        for (int l = 1; l < lines.size(); ++l)
        {
            string line = lines[l];
            if (line != "")
            {
                vector<string> pz = Util::Split(line, " ", false);
                
                //NOTE: we assume node ids in the file are indexed at 1, not 0
                pair<int, int> e = make_pair( Util::ToInt(pz[0]) - 1, Util::ToInt(pz[1]) - 1 );
                
                outs[e.first].insert(e.second);
                ins[e.second].insert(e.first);
                    
                all_edges.push_back( e );
                edge_indexes[ e ] = all_edges.size() - 1;
                
                edge_weights.push_back( Util::ToInt(pz[2]) );
            }
        }


    }
    
    
    void read_graph( string filename )
    {
        
        vector<string> lines = Util::GetFileLines(filename);
        
        //first line should have 
        //nbnodes nbarcs
        vector<string> px = Util::Split(lines[0], " ", false);
        
        nb_nodes = Util::ToInt(px[0]);
        
        cout<<"nbnodes="<<nb_nodes<<endl;
        
        //we assume node indices are indexed at 1 (ie first node has index 1, not 0)
        ins.resize(nb_nodes);
        outs.resize(nb_nodes);
        
        for (int l = 1; l < lines.size(); ++l)
        {
            string line = lines[l];
            if (line != "")
            {
                vector<string> pz = Util::Split(line, " ", false);
                
                //NOTE: we assume node ids in the file are indexed at 1, not 0
                pair<int, int> e = make_pair( Util::ToInt(pz[0]) - 1, Util::ToInt(pz[1]) - 1 );
                
                outs[e.first].insert(e.second);
                ins[e.second].insert(e.first);
                    
                all_edges.push_back( e );
                edge_indexes[ e ] = all_edges.size() - 1;
                
                edge_weights.push_back( Util::ToInt(pz[2]) );
            }
        }
        
    }
    
    
    
    void print_cycle( deque<int> cyc )
    {
        for (int c : cyc)
        {
            cout<<" "<<c;
        }
        cout<<endl;
    }
    
    
    
    void solve_fas()
    {
        deque<int> test_cycle;
        bool testret = get_cycle(test_cycle);
        
        if (testret)
        {
            cout<<"Graph has a cycle ";
            for (int v : test_cycle)
                cout<<v + 1<<" ";
            cout<<endl;
        }
        else 
        {
            cout<<"Graph has no cycle!"<<endl;
        }
        
        //first calculate max weight for hard clauses
        int tid = 0;
        for (int i = 0; i < all_edges.size(); ++i)
        {
            max_weight += edge_weights[i];    //TODO: read edge weights
        }
        max_weight *= 2;
        
        
        //get all triangles, which become hard clauses
        vector< deque<int> > current_cycles;    // = get_triangles();
        
        //params are max_per_source, minlen, maxlen, max_per_source, max_ttl
        get_short_cycles( current_cycles, 3, 3, 100 );
        get_short_cycles( current_cycles, 4, 4, 100 );
        get_short_cycles( current_cycles, 5, 5, 100 );
        
        if (current_cycles.size() > 0)
        {
            cout<<"Last cycle added = ";
            print_cycle( current_cycles[ current_cycles.size() - 1 ] );
            cout<<endl;
        }
        
        
        vector< vector<int> > hard_clauses;
        vector< vector<int> > soft_clauses;

        bool done = false;
        
        
        EvalMaxSAT solver;
        umap< int, int > solverid_to_edge_index;
        umap< int, vector<int> > edge_index_to_solverids;
        
        
        /**********************************************/
        //setup solver vars
        for (int i = 0; i < all_edges.size(); ++i)
        {
            int tid = solver.newVar();
            edge_index_to_solverids[i].push_back(tid);
            solverid_to_edge_index[tid] = i;
            
            
            //add one singleton clause per edge.  true means "delete", false means "don't delete"
            if (edge_weights[i] > 10000)    //TODO: cleaner thing for high weight
            {
                hard_clauses.push_back( { -tid } );    
            }
            else 
            {
                soft_clauses.push_back( { -tid } );
                
                for (int w = 1; w < edge_weights[i]; ++w)
                {
                    int tidw = solver.newVar();
                    edge_index_to_solverids[i].push_back(tid);
                    solverid_to_edge_index[tidw] = i;

                    soft_clauses.push_back( { -tidw } );
                    soft_clauses.push_back( {-tid, tidw} );
                }
                
            }
        }

        //done with vars at this point
        
        for (vector<int> hclause : hard_clauses)
        {
            solver.addClause(hclause);
        }
        
        for (vector<int> sclause : soft_clauses)
        {
            /*cout<<"add soft clause ";
            for (int x : sclause) cout<<" "<<x;
            cout<<endl;*/
            solver.addWeightedClause(sclause, 1);
        }
        
        
        vector<int> clause_cycle;
        for (deque<int> cyc : current_cycles)
        {
            clause_cycle.clear();
            //add clause for cycle
            for (int trindex = 0; trindex < cyc.size(); ++trindex)
            {
                pair<int, int> arc = make_pair(cyc[trindex], cyc[ (trindex + 1) % cyc.size() ]);
                int e_index = edge_indexes[ arc ];
                clause_cycle.push_back( edge_index_to_solverids[e_index][0] );
            }
            
            //cout<<"adding cycle clause"<<endl;
            solver.addClause( clause_cycle );  //no weight = hard clause   //, max_weight );
        }
        
        
        cout<<"Number of cycles is "<<current_cycles.size()<<endl;

        
        int nb_cycles = current_cycles.size();

        
        while (!done)
        {    
            
            
            //to_wcnf( "scc1.wcnf", current_cycles );
            
            
            cout<<"ready to solve"<<endl;
            
            //solve 
            if(!solver.solve()) 
            {
                cout << "UNSATISFIABLE" << endl;
                return;
            }
            
            last_cost = solver.getCost();
            cout<<"cost="<<last_cost<< endl;
            
            vector<int> edges_to_delete;
            for (int i = 0; i < all_edges.size(); ++i)
            {
                int tid = edge_index_to_solverids[i][0];
                
                bool v = solver.getValue(tid);
                //cout<<"tid="<<tid<<" val="<<v<<endl;
                
                if (solver.getValue(tid))
                {
                    edges_to_delete.push_back(i);
                }
            }
            
            
            for (int i = 0; i < edges_to_delete.size(); ++i)
            {
                pair<int, int> e = all_edges[ edges_to_delete[i] ];
                delete_edge( e.first, e.second );
                
                //cout<<"Delete ("<<e.first + 1<<","<<e.second + 1<<")"<<endl;
            }
            
            
            deque<int> next_cycle;
            bool ret = get_cycle(next_cycle);
            
            //still a cycle somewhere
            if (ret)
            {
                cout<<"Still has cycle ";
                for (int v : next_cycle)
                    cout<<v<<" ";
                cout<<endl;
                
                vector<deque<int>> newcycles = {next_cycle};
                
                get_short_cycles( newcycles, 3, 5, 5, 100 );
                
                
                
                //current_cycles.push_back(next_cycle);
                
                cout<<"Adding "<<newcycles.size()<<" new cycles"<<endl;
                nb_cycles += newcycles.size();
                
                //add clause cycle
                for (deque<int> cyc: newcycles)
                {
                    vector<int> clause_next_cycle;
                    for (int cindex = 0; cindex < cyc.size(); ++cindex)
                    {
                        pair<int, int> arc = make_pair(cyc[cindex], cyc[ (cindex + 1) % cyc.size() ]);
                        int e_index = edge_indexes[ arc ];
                        clause_next_cycle.push_back( edge_index_to_solverids[e_index][0] );
                    }
                    solver.addClause(clause_next_cycle);
                }
                
                
                for (int i = 0; i < edges_to_delete.size(); ++i)
                {
                    pair<int, int> e = all_edges[ edges_to_delete[i] ];
                    undelete_edge( e.first, e.second );
                    
                    //cout<<"Undelete ("<<e.first + 1<<","<<e.second + 1<<")"<<endl;
                }
                
                
                cout<<"Number of cycles is now "<<nb_cycles<<endl;
                
            }
            else 
            {
                cout<<"no more cycles"<<endl;
                done = true;
            }
            
            
            
            edges_to_delete.clear();
            
        }
        
        
        
    }
    
    
    
    
    
    
    

    //returns all the triangles, each element is a vector with the three node indexes forming the triangle, in order
    vector< deque<int> > get_triangles()
    {
        vector< deque<int> > retvec;
        for (int i = 0; i < nb_nodes; ++i)
        {
            for (int j : outs[i])
            {
                
                for (int k : outs[j])
                {
                    if (outs[k].count(i))    //count acts as contains
                    {
                        deque<int> triangle{i, j, k};
                        retvec.push_back(triangle);
                    }
                }
            }
        }
        
        return retvec;
        
    }
    
    
    
    //does not gurantee it will find them all!
    void get_short_cycles( vector< deque<int> > &cycles, int minlen, int maxlen, int max_per_source = 999999, int max_nb_cycles = 999999 )
    {
        
        for (int i = 0; i < nb_nodes; ++i)
        {
            bfs_cycles( i, cycles, minlen, maxlen - 1, max_per_source );
            if (cycles.size() > max_nb_cycles)
                return;
        }
        
    }
    
    void bfs_cycles( int startnode, vector< deque<int> > &cycles, int minlen, int maxdepth, int max_nb_cycles = 999999 )
    {
        deque<int> queue;
        umap<int,int> visited_depths;
        umap<int,int> visited_pred;
        
        int nb_added = 0;
        
        queue.push_back(startnode);
        visited_depths[startnode] = 0;
        
        //yeah a while true but whatever, I'm not a software engineer anymore
        while (true)
        {
            if (queue.empty())
                break;
            
            int v = queue.front();
            queue.pop_front();
            
            if (visited_depths[v] >= maxdepth)
                continue;
            
            for (int w : outs[v])
            {
                //NOTE: we only visit those w > start_node.  The reason being that we consider the start of a cycle as the min index.
                if (!visited_depths.count(w) && w > startnode)
                {
                    visited_depths[w] = visited_depths[v] + 1;
                    visited_pred[w] = v;
                    queue.push_back(w);

                    //found a short cycle
                    if (visited_depths[w] >= minlen - 1 && outs[w].count(startnode))
                    {
                        deque<int> cyc = {w};
                        int cur = w;
                        while ( cur != startnode )
                        {
                            cur = visited_pred[cur];
                            cyc.push_front(cur);
                        }
                        cycles.push_back(cyc);
                        nb_added += 1;
                        
                        if (nb_added >= max_nb_cycles)
                        {
                            return;
                        }
                    }
                }
            }
        }
    }
    
    
    
    void delete_edge(int i, int j)
    {
        deleted_edges.insert( make_pair(i, j) );
        outs[i].erase(j);
        ins[j].erase(i);
    }
    
    void undelete_edge(int i, int j)
    {
        if (deleted_edges.count( make_pair(i, j) ))
        {
            deleted_edges.erase( make_pair(i, j) );
            outs[i].insert(j);
            ins[j].insert(i);
        }
    }
    
     
    
    bool get_cycle(deque<int> &cycle)
    {
        
        for (int i = 0; i < nb_nodes; ++i)
        {
            uset<int> visited;
            if (dfs_cycle(i, i, cycle, visited))
            {
                return true;
            }
        }
        
        return false;
        
    }
    
    
    //return true if topo_sort succeeds, false otherwise
    bool get_topo_sort( vector< int > &topo_sort )
    {
        vector< int > indegs;
        indegs.resize(nb_nodes);
        
        list< int > zeroin_queue;
        uset< int > unhandled;
        
        for (int i = 0; i < nb_nodes; ++i)
        {
            indegs[i] = ins[i].size();
            
            if (indegs[i] == 0)
                zeroin_queue.push_back(i);
            else 
                unhandled.insert(i);
        }
        
        
        while (zeroin_queue.size() > 0)
        {
            int v = zeroin_queue.front();
            zeroin_queue.pop_front();
            
            topo_sort.push_back(v);
            unhandled.erase(v);
            
            for (int w : outs[v])
            {
                if (unhandled.count(w))
                {
                    indegs[w]--;
                    
                    if (indegs[w] == 0)
                    {
                        zeroin_queue.push_back(w);
                    }
                }
                else 
                {
                    cout<<"ERROR: topo sort found a backwards arc"<<endl;
                }
                
            }
        }
        
        
        return (topo_sort.size() == nb_nodes);
    }
    
    
    
    bool dfs_cycle(int cur_node, int start_node, deque<int> &outvec, uset<int> &visited )
    {
        if (visited.count(cur_node))
        {
            return false;
        }
        
        visited.insert(cur_node);
        
        if (cur_node != start_node && outs[cur_node].count(start_node))
        {
            outvec.push_front(cur_node);
            return true;
        }
        
        
        for (int j : outs[cur_node])
        {
            if (dfs_cycle(j, start_node, outvec, visited))
            {
                outvec.push_front(cur_node);
                return true;
            }
        }
        return false;
        
    }
};

string solve_instance(string infilelines) {

    Digraph g;
    g.read_graph_lines(infilelines);
    g.solve_fas();
        
    vector<int> topo_sort;
    g.get_topo_sort(topo_sort);

    string strout = "";
    for (int v : topo_sort)
    {
        strout += Util::ToString(v) + " ";
    }

    return strout;
}

int main(int argc, char* argv[])
{
    vector<string> raw_args;
    umap<string, string> all_args;

    for (int i = 1; i < argc; ++i)
    {
        string bef = "", aft = "";
        string s = argv[i];
        const size_t pos = s.find('=');
        if (pos != string::npos) {
            bef = s.substr(0, pos);
            aft = s.substr(pos + 1);
            all_args[bef] = aft;
        }
    }
    
    if (!all_args.count("--in") || !all_args.count("--out"))
    {
        cout<<"You need to specify --in=[ingraph] and --out=[outfile]"<<endl;    
        return 0;        
    }
    
    string infile = all_args["--in"];
    string outfile = all_args["--out"];
    string outgraph_file = "";
    if (all_args.count("--outgraph"))
        outgraph_file = all_args["--outgraph"];
    
    Digraph g;
    g.read_graph(infile);
    g.solve_fas();
    
    
    if (outgraph_file != "")
        g.write_graph( outgraph_file );
        
    vector<int> topo_sort;
    g.get_topo_sort(topo_sort);
    
    g.write_solution_info(topo_sort, outfile);
    
    return 0;
}
