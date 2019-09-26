#ifndef DIGRAPH_H
#define DIGRAPH_H

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <stack>

// function std::min()
#include <algorithm>

/**
 * Template class for directed graphs.
 *
 * This template class contains functions that are useful for directed graphs like the computation of strongly
 * connected components.
 */
template<typename T>
class DiGraph {
public:
  /**
  * Add a node v to the digraph.
  * @param v node to be added.
  */
  void addNode(T v);

  /**
  * Add an edge (u,v) to the digraph and the nodes u and v if they do not exist yet.
  * @param u initial node of the edge.
  * @param v terminal node of the edge.
  */
  void addEdge(T u, T v);


  /**
  * Returns the strongly connected components of the digraph.
  * @return strongly connected components.
  */
  std::vector<std::set<T>>
  stronglyConnectedComponents() const;

  /**
  * Returns the first strongly connected components without out-edges to other s.c.c. that is found or an empty
  * list if no such s.c.c. is found.
  * @return strongly connected component without out-edges.
  */
  std::vector<T>
  getSccWithoutOutedges() const;

  /**
  * Returns the adjecency list of the digraph.
  * @return adjecency list.
  */
  const std::unordered_map<T, std::vector<T>>&
  getAdjacency() const { return m_adjList; };
  
  /**
  * Returns the (out-)neighbors of a node v.
  * @param v node
  * @return (out-)neighbors of v.
  */
  const std::vector<T>&
  getNeighbors(T v) const { return m_adjList.at(v); };

private:
  std::unordered_map<T, std::vector<T>> m_adjList; /*!< the adjacency list of the digraph */
};

template<typename T>
void
DiGraph<T>::addNode(T v){
  if(m_adjList.find(v) == m_adjList.end()){
    m_adjList.insert(std::pair<T, std::vector<T>>(v, std::vector<T>()));
  }
}

template<typename T>
void
DiGraph<T>::addEdge(T u, T v){

  // check if u is in the graph and if not insert it first
  if(m_adjList.find(u) == m_adjList.end()){
    m_adjList.insert(std::pair<T, std::vector<T>>(u, std::vector<T>()));
  }

  // check if v is in the graph and if not insert it first
  if(m_adjList.find(v) == m_adjList.end()){
    m_adjList.insert(std::pair<T, std::vector<T>>(v, std::vector<T>()));
  }
  //
  // insert the directed edge (check for existence not implemented!)
  m_adjList[u].push_back(v);
}

template<typename T>
std::vector<std::set<T>>
DiGraph<T>::stronglyConnectedComponents() const {
  // Uses Tarjan's algorithm with Nuutila's modifications.
  // Nonrecursive version of algorithm.
  // (reimplemented from Python package "NetworkX")

  auto sccs = std::vector<std::set<T>>();

  //auto nbrs = std::map<T, std::vector<T>*>();
  auto preorder = std::map<T, int>();
  auto lowLink = std::map<T, int>();
  auto sccFound = std::set<T>();
  auto sccStack = std::stack<T>();
  auto scc = std::stack<T>();
  int i = 0;

  for(auto source = m_adjList.begin(); source != m_adjList.end(); ++source){
    if(sccFound.find(source->first) == sccFound.end()){
      auto stack = std::stack<T>();
      stack.push(source->first);

      while(!stack.empty()){
        T v = stack.top();

        if(preorder.find(v) == preorder.end()){
          ++i;
          preorder[v] = i;
        }

        int done = 1;
        for(auto w : m_adjList.at(v)){
          if(preorder.find(w) == preorder.end()){
            stack.push(w);
            done = 0;
            break;
          }
        }

        if(done == 1){
          lowLink[v] = preorder[v];

          for(auto w : m_adjList.at(v)){
            if(sccFound.find(w) == sccFound.end()){
              if(preorder[w] > preorder[v]){
                lowLink[v] = std::min(lowLink[v], lowLink[w]);
              } else {
                lowLink[v] = std::min(lowLink[v], preorder[w]);
              }
            }
          }

          stack.pop();

          if(lowLink[v] == preorder[v]){
            sccFound.insert(v);
            sccs.push_back(std::set<T>());
            auto& scc = sccs.back();
            scc.insert(v);

            while((!sccStack.empty()) &&
                  (preorder[sccStack.top()] > preorder[v])){
              T k = sccStack.top();
              sccStack.pop();
              sccFound.insert(k);
              scc.insert(k);
            }

          } else {
            sccStack.push(v);
          }
        }
      }
    }
  }

  return sccs;
}

template<typename T>
std::vector<T>
DiGraph<T>::getSccWithoutOutedges() const {

  auto sccs = stronglyConnectedComponents();
  auto result = std::vector<T>();

  for(auto scc : sccs){
    bool noOutedges = true;

    for(auto v : scc){
      for(auto w : m_adjList.at(v)){
        if(scc.find(w) == scc.end()){
          noOutedges = false;
        }
      }
    }

    if(noOutedges){
      for(auto v : scc){
        result.push_back(v);
      }
      break;
    }
  }

  return result;
}

#endif /* DIGRAPH_H */
