#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <fstream>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Point_set_2.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef K::Point_2 Point;
//typedef CGAL::Delaunay_triangulation_2<K> Triangulation;
// As we only consider finite vertices and edges
// we need the following filter
template <typename T>
struct Is_finite {
    const T* t_;
    Is_finite()
    : t_(NULL)
    {}
    Is_finite(const T& t)
    : t_(&t)
    { }
    template <typename VertexOrEdge>
    bool operator()(const VertexOrEdge& voe) const {
        return ! t_->is_infinite(voe);
    }
};


//jądro obliczeń arytmetyczych

//opisy wierzchołków
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;

//reprezentacja triangulacji (ściany i wierzchołki)
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay0;
typedef CGAL::Point_set_2<K, Tds> Delaunay;
typedef Is_finite<Delaunay0> Filter;
typedef boost::filtered_graph<Delaunay0,Filter,Filter> Finite_triangulation;
typedef boost::graph_traits<Finite_triangulation>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Finite_triangulation>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Finite_triangulation>::edge_descriptor edge_descriptor;



using namespace std;



typedef Delaunay::Point Point;                 //punkt triangulacji Delaunay'a
typedef Delaunay::Edge_iterator Edge_iterator; //iterator po krawedziach trian. Delaunaya
typedef Delaunay::Vertex_handle Vertex_handle; //iterator po punktach trian. Delaunaya

typedef boost::graph_traits<Delaunay0>::vertex_descriptor vertex_descriptor1;


// The BGL makes use of indices associated to the vertices
// We use a std::map to store the index
typedef std::map<vertex_descriptor,int> VertexIndexMap;

class FMap {
public:
    typedef vertex_descriptor1                        key_type;
    typedef std::pair<const vertex_descriptor1, int>        value_type;
    int& operator[](const key_type& key) const {
        return key->info();
    }
};
FMap vertex_id_map;

// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<FMap> VertexIdPropertyMap;

VertexIdPropertyMap vertex_index_pmap(vertex_id_map);




struct MPoint {
    Point point;
    int id;
    vector<int> neighbours;
};


struct MultiPoint {
    Point point;
    int id;
    map<int, int> multiNeighbours;
};
struct PointDir {
    Point point;
    int id;
    set<int> neighbours;
};

struct Edge {
    int from;
    int to;
};



class Graph {
public:
    vector<MPoint> points;
    
    Graph(int points_num): points(points_num) {
        points.resize(points_num);
    }
    void add_edge_undir(int from, int to) {
        points[from].neighbours.push_back(to);
        points[to].neighbours.push_back(from);
    }
    void add_edge_dir(int from, int to) {
        points[from].neighbours.push_back(to);
        points[to].neighbours.push_back(from);
    }
};


class DirGraph {
public:
    vector<PointDir> points;
    
    DirGraph(int points_num): points(points_num) {}
    void add_edge_dir(int from, int to) {
        points[from].neighbours.insert(to);
    }
    
};


class MultiGraph {
public:
    vector<MultiPoint> points;
    MultiGraph(int points_num): points(points_num) {
        points.resize(points_num);
    }
    void add_edge_undir(int from, int to) {
        auto neigh = points[from].multiNeighbours;
        auto it = neigh.find(from);
        if (it != neigh.end()) {
            it->second += 1;
            points[to].multiNeighbours[from] += 1;
        }
        else {
            neigh[to] = 1;
            points[to].multiNeighbours[from] = 1;
        }
    }
};


MultiGraph mst(vector<pair<Point, int>>& points) {
    Delaunay t;
    t.insert(points.begin(), points.end());

    Filter is_finite((Delaunay0)t);
    Finite_triangulation ft((Delaunay0)t, is_finite, is_finite);
    
    std::list<edge_descriptor> mstList;
    boost::kruskal_minimum_spanning_tree(ft,
                                         std::back_inserter(mstList),
                                         vertex_index_map(vertex_index_pmap));
    
    MultiGraph result(points.size());
    for(edge_descriptor e: mstList) {
        vertex_descriptor1 svd = source(e, t);
        vertex_descriptor1 tvd = target(e, t);
        result.add_edge_undir(svd->info(),tvd->info());
    }
    return result;
}

void dup_leaf_edges(MultiGraph& graph) {
    for(MultiPoint& mpoint: graph.points) {
        if(mpoint.multiNeighbours.size() == 1) {
            auto it = mpoint.multiNeighbours.begin();
            it->second++;
            graph.points[it->first].multiNeighbours[mpoint.id]++;
        }
    }
}

vector<pair<Point, int>> odd_deg(MultiGraph& graph){
    vector<pair<Point, int>> result;
    for(MultiPoint& mpoint: graph.points) {
        int sum = 0;
        for (auto& kv : mpoint.multiNeighbours) sum += kv.second;
        if(sum%2) result.emplace_back(mpoint.point, mpoint.id);
    }
    return result;
}


vector<Edge> convex_hull(Delaunay& del) {
    vector<Edge> result;
    Delaunay::Vertex_circulator start, done;
    start = del.incident_vertices(del.infinite_vertex());
    done = start;
    if (start != 0)
    {
        auto from = start;
        auto to = ++start;
        do {
            result.push_back({from->info(), to->info()});
            from++;
            to++;
        } while(from != done);
    }
    return result;
}



vector<Edge> sub_tsp(vector<pair<Point, int>> points) {
    //TODO
    Delaunay del;
    //insert
    vector<Edge> path = convex_hull(del);
    return vector<Edge>();
    
}

inline int euclid(int v, int w, const MultiGraph& g) {
    double x = g.points[v].point.x() - g.points[w].point.x();
    double y = g.points[v].point.y() - g.points[w].point.y();
    return floor(sqrt(x*x+y*y) + 0.5);
}


vector<Edge> better_match(vector<Edge> edges, const MultiGraph& g) {
    vector<Edge> result;
    double sums[] = {0,0};
    for (int i = 0; i< edges.size(); i++) {
        Edge edge = edges[i];
        sums[i%2] += euclid(edge.from, edge.to, g);
    }
    bool secondBetter = sums[1] < sums[0];
    for (int i = secondBetter; i< edges.size(); i+=2) {
        result.push_back(edges[i]);
    }
    return result;
}

vector<Edge> euler_cycle(MultiGraph& g) {
    vector<Edge> acc;
    int curr = 0;
    while(true) {
        auto& neighs = g.points[curr].multiNeighbours;
        if (neighs.empty()) break;
        auto a = neighs.begin();
        int neigh = a->first;
        acc.push_back({curr, neigh});
        auto mirrorNeigh = g.points[neigh].multiNeighbours;
        a->second--;
        mirrorNeigh[curr]--;
        if(a->second == 0) {
            neighs.erase(neigh);
            mirrorNeigh.erase(curr);
        }
        curr = neigh;
    }
    return acc;
}

vector<Edge> shortcuts(vector<Edge>& cycle) {
    //TODO
    return vector<Edge>();
}

vector<Edge> tsp(vector<pair<Point, int>>& points) {
    MultiGraph g1 = mst(points);
    dup_leaf_edges(g1);
    vector<pair<Point, int>> odd_points = odd_deg(g1);
    vector<Edge> odd_match = better_match(sub_tsp(odd_points), g1);
    for(Edge e : odd_match) {
        g1.add_edge_undir(e.from, e.to);
    }
    vector<Edge> euler = euler_cycle(g1);
    return shortcuts(euler);
}


int main1(int argc,char* argv[])
{
    
    Delaunay t;
    Filter is_finite((Delaunay0)t);
    Finite_triangulation ft((Delaunay0)t, is_finite, is_finite);
    
    vertex_iterator vit, ve;
    // Associate indices to the vertices
    int index = 0;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vit and ve
    for(boost::tie(vit,ve)=boost::vertices(ft); vit!=ve; ++vit ){
        vertex_descriptor  vd = *vit;
    }
    // We use the default edge weight which is the squared length of the edge
    // This property map is defined in graph_traits_Triangulation_2.h
    // In the function call you can see a named parameter: vertex_index_map
    std::list<edge_descriptor> mst;
    boost::graph_traits<Delaunay0>::vertex_descriptor d;
    boost::graph_traits<CGAL::Delaunay_triangulation_2<K, Tds>>::vertex_descriptor d1;
//    boost::graph_traits<CGAL::Point_set_2<K, Tds>>::vertex_descriptor d2;
    

    
    
    
    std::cout << num_vertices(t);
    edges(t);
    auto vrts = vertices(t);
    auto vrtx = *(vrts.begin());
    auto x = vrtx.operator*();
   
    std::cout <<"YOLO"<< std::endl;
    boost::kruskal_minimum_spanning_tree(ft,
                                         std::back_inserter(mst),
                                         vertex_index_map(vertex_index_pmap));
    std::cout << "The edges of the Euclidean mimimum spanning tree:" << std::endl;
    for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it){
        edge_descriptor ed = *it;
        vertex_descriptor svd = source(ed,t);
        vertex_descriptor tvd = target(ed,t);
//        Triangulation::Vertex_handle sv = svd;
//        Triangulation::Vertex_handle tv = tvd;
    }


    return 0;
}
