#include "gtest/gtest.h"

#include <vector>

#include "DiGraph.h"

namespace {

TEST(DiGraphTest, Neighbors) {
  auto dg = DiGraph<int*>();

  int i1 = 1;
  int i2 = 2;
  int i3 = 3;
  int i4 = 4;

  dg.addEdge(&i1, &i2);
  dg.addEdge(&i2, &i1);
  dg.addEdge(&i1, &i3);

  dg.addEdge(&i2, &i3);

  dg.addNode(&i4);

  auto neighbors = dg.getNeighbors(&i2);
  std::vector<int*> neighbors_expected{ &i1, &i3 };

  EXPECT_EQ(neighbors, neighbors_expected);
}

TEST(DiGraphTest, SCCS) {
  auto dg = DiGraph<int*>();

  int i1 = 1;
  int i2 = 2;
  int i3 = 3;
  int i4 = 4;
  int i5 = 5;

  dg.addEdge(&i1, &i2);
  dg.addEdge(&i2, &i3);
  dg.addEdge(&i3, &i1);
  dg.addEdge(&i1, &i3);

  dg.addEdge(&i4, &i3);
  dg.addEdge(&i4, &i5);
  dg.addEdge(&i5, &i4);

  auto neighbors = dg.getNeighbors(&i1);
  std::vector<int*> neighbors_expected{ &i2, &i3 };

  /*
  auto sccs = dg.stronglyConnectedComponents();
  for(auto scc : sccs){
    for(auto item : scc){
      std::cout << *item << "  ";
    }
    std::cout << std::endl;
  }
  */

  /*
  auto sccWithoutOutedges = dg.getSccWithoutOutedges();
  std::cout << "scc without out-edges: ";
  for(auto item : sccWithoutOutedges){
    std::cout << *item << "  ";
  }
  std::cout << std::endl;
  */

  EXPECT_EQ(neighbors, neighbors_expected);
}

}  // namespace
