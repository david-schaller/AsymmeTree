#ifndef GENE_H
#define GENE_H

#include <string>

/**
 * Class for genes.
 *
 * This class holds information of a gene.
 */
class Gene {
public:
  /**
   * Class constructor which creates a new instance.
   * @param identifier identifier of the gene.
   * @param idx (unique) index of the gene.
   */
  Gene(std::string identifier, size_t idx)
    : m_identifier(identifier)
    , m_idx (idx) { };

  /**
   * Returns the identifier of the gene.
   */
  std::string
  getIdentifier() const { return m_identifier; };

  /**
   * Returns the identifier of the genes' species.
   */
  const std::string&
  getSpecies() const { return m_species; };

  /**
   * Returns the index of the gene.
   */
  size_t
  getIndex() const { return m_idx; };

  /**
   * Returns the index of the subtree under the root of the species tree.
   */
  size_t
  getSubtree() const { return m_subtreeIdx; };

  /**
   * Set the species of the gene.
   * @param species (new) species identifier of the gene.
   */
  void
  setSpecies(std::string species) { m_species = species; }

  /**
   * Set the index of the subtree under the root of the species tree.
   * @param subtree index of the subtree.
   */
  void
  setSubtree(size_t subtree) { m_subtreeIdx = subtree; }

  /**
   * Compares if two genes have an equal identifier.
   * @param otherGene gene to compare with.
   */
  bool operator==(const Gene &otherGene) const
    {
        return (m_identifier == otherGene.m_identifier);
    }

private:
  std::string m_identifier;   /*!< the genes' identifier */
  std::string m_species;      /*!< the identifier of the genes' species*/
  size_t m_idx;               /*!< unique index */
  size_t m_subtreeIdx;        /*!< the index of the subtree under the root of the species tree */
};

#endif /* GENE_H */
