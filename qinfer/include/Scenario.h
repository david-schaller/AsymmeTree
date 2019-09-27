#ifndef SCENARIO_H
#define SCENARIO_H

#include <vector>
#include <deque>
#include <unordered_map>
#include <string>

#include "Gene.h"
#include "Matrix.h"
#include "Tree.h"
#include "Parameters.h"
#include "Benchmark.h"

/**
 * Class for the scenario.
 *
 * This class holds the information on the genes, species, etc.
 */
class Scenario {
public:
  /**
   * Class constructor which creates a new instance.
   * @param ptrParam pointer to the Parameters instance.
   * @param bm pointer to the Benchmark instance.
   */
  Scenario(Parameters* ptrParam, Benchmark* bm = nullptr)
  : m_ptrParam(ptrParam)
  , m_genes()
  , m_distanceMatrix()
  , m_geneAssignments()
  , m_speciesGenes()
  , m_STree()
  , m_benchmark(bm) { };

  /**
   * Replace the list of genes by a new list.
   * @param g new list of genes to add.
   */
  void
  addGenes(std::deque<Gene>& g);

  /**
   * Parses all input files.
   */
  void
  parseFiles();

  /**
   * Returns the species of the gene with index i.
   * @param i index of the gene.
   */
  const std::string&
  getGeneSpecies(int i) const { return m_genes[i].getSpecies();};

  /**
   * Returns the list of genes.
   */
  std::deque<Gene>&
  getGenes() { return m_genes; };

  /**
   * Returns the number of genes.
   */
  size_t
  getGeneNumber() const { return m_genes.size(); };

  /**
   * Returns the index of a gene.
   */
  Gene*
  getGenePtr(size_t idx) { return &(m_genes[idx]); }

  /**
   * Returns the list of genes of a specific species.
   * @param species the species' (unique) identifier.
   */
  const std::vector<Gene*>&
  getSpeciesGenes(std::string species) const { return m_speciesGenes.at(species); };

  /**
   * Returns the number of species.
   */
  size_t
  getSpeciesNumber() const { return m_speciesGenes.size(); };

  /**
   * Returns an entry of the distance matrix.
   * @param row row index.
   * @param column column index.
   */
  const double&
  getDistance(size_t row, size_t column) const { return m_distanceMatrix.at(row, column); };

  /**
   * Returns a reference to the species tree.
   */
  Tree&
  getSTree() { return m_STree; };

  /**
   * Returns a list of outgroup genes for a specific gene (w.r.t. the root of the species tree).
   * @param genePtr pointer to the gene.
   */
  const std::vector<Gene*>&
  getOutgroups(Gene* genePtr) const;

  /**
   * Returns a list of species in a subtree of the root of the species tree.
   * @param subtreeIdx index of the subtree under the root.
   */
  const std::vector<std::string>&
  getSpeciesSubtree(size_t subtreeIdx) const;

private:
  Parameters* m_ptrParam;                                             /*!< pointer to the Parameters instance */
  std::deque<Gene> m_genes;                                           /*!< list of genes */
  Matrix<double> m_distanceMatrix;                                    /*!< distance matrix */
  std::unordered_map<std::string, Gene*> m_geneAssignments;           /*!< maps gene identifiers to the genes */
  std::unordered_map<std::string, std::vector<Gene*>> m_speciesGenes; /*!< maps species identifiers to lists of genes */
  Tree m_STree;                                                       /*!< species tree */
  Benchmark* m_benchmark;                                             /*!< pointer to the Parameters instance */

  std::vector<std::vector<std::string>> m_STreeSubtrees;              /*!< list of lists of species identifiers
                                                                      belonging to the same m_STreeSubtree of the
                                                                      species tree root and corresponding outgroup
                                                                      genes */
  std::vector<std::vector<Gene*>> m_outgroups;                        /*!< outgroup lists */

  /**
   * Splits a string at whitespaces.
   * @param s string to split.
   */
  std::vector<std::string> split(std::string s);

  /**
   * Parses the distance matrix.
   * @param filepath path to the distance matrix file.
   */
  void parseDistanceMatrix(std::string filepath);

  /**
   * Parses the species to genes mapping file.
   * @param filepath path to species to genes mapping file.
   */
  void parseSpeciesGenes(std::string filepath);

  /**
   * Parses the species tree file.
   * @param filepath path to species tree file.
   * @param species tree mode (whole tree or subtrees)
   */
  void parseSTreeSubtrees(std::string filepath, bool subtreeFiles);

  /**
   * Rebuild the map mapping gene identifiers to the genes.
   */
  void rebuildMap();

  /**
   * Parses a row of the distance matrix.
   * @param row string to parse.
   * @param rowIdx index of the row in the matrix.
   */
  void parseDistanceMatrixRow(std::string row, size_t rowIdx);

  /**
   * Parses a line of the species to genes mapping file.
   * @param line string to parse.
   */
  void parseSpeciesGenesLine(std::string line);

  /**
   * Parses a Newick tree.
   * @param line string to parse.
   */
  void parseNewickAndCheck(std::string line);

  /**
   * Parses a line with identifiers of the species in a subtree under the root.
   * @param line string to parse.
   */
  void parseSTreeSubtreeLine(std::string line);

  /**
   * Checks if all species could be assigned to a subtree.
   */
  void checkSpeciesAvailability();

  /**
   * Finds the (root) outgroups for the subtrees under the root of the species tree.
   */
  void buildOutgroupInfo();
};

#endif /* SCENARIO_H */
