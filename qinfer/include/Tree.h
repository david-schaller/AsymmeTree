#ifndef TREE_H
#define TREE_H

#include <memory>
#include <vector>
#include <string>

/**
 * Class for (species) tree nodes.
 *
 * This class contains functions concerning the nodes of a (species) tree.
 */
class TreeNode {
public:
  friend class Tree;

  /**
   * Class constructor which creates a new instance.
   * @param value value/identifier of the node.
   */
  TreeNode(std::string value)
    : m_value(value)
    , m_nodeIdx(-1)
    , m_leafIdx(-1) { };

  /**
   * Returns the value of the node.
   */
  std::string
  getValue() const { return m_value; };

  /**
   * Returns the unique node index.
   */
  size_t
  getNodeIdx() const { return m_nodeIdx; };

  /**
   * Returns the unique leaf index (e.g. for the last common ancestor matrix).
   */
  size_t
  getLeafIdx() const { return m_leafIdx; };

  /**
   * Returns true if the node is not the root.
   */
  bool
  hasParent() const { return !m_parent.expired(); };

  /**
   * Returns true if the node is not a leaf.
   */
  bool
  hasChildren() const { return !m_children.empty(); };

  /**
   * Returns the children of the node.
   */
  const std::vector<std::shared_ptr<TreeNode>>&
  getChildren() const { return m_children; };

  /**
   * Returns the leaves in the subtree rooted at this node.
   */
  const std::vector<std::shared_ptr<TreeNode>>&
  getLeaves() const { return m_leaves; };

  /**
   * Returns the parent of the node.
   */
  const std::weak_ptr<TreeNode>&
  getParent() const { return m_parent; };

  /**
   * Adds a child to the node.
   */
  std::shared_ptr<TreeNode>
  addChild(std::string value, std::shared_ptr<TreeNode> parent);

private:
  std::string m_value;    /*!< value */
  size_t m_nodeIdx;       /*!< unique node index */
  size_t m_leafIdx;       /*!< unique leaf index */
  std::vector<std::shared_ptr<TreeNode>> m_children;  /*!< list of children */
  std::weak_ptr<TreeNode> m_parent;                   /*!< parent */
  std::vector<std::shared_ptr<TreeNode>> m_leaves;    /*!< leaves in the subtree */
};


/**
 * Class for (species) trees.
 *
 * This class contains functions concerning (species) trees.
 */
class Tree {
public:
  /**
   * Class constructor which creates a new instance and a root of the tree.
   */
  Tree()
    : m_root(std::shared_ptr<TreeNode>())
    , m_nodeCounter(0)
    , m_leafCounter(0)
    , m_preorderUpToDate(false) { };

  /**
   * Class constructor which creates a new instance for an existing root.
   * @param root existing root node.
   */
  Tree(std::shared_ptr<TreeNode> root)
    : m_root(root)
    , m_nodeCounter(0)
    , m_leafCounter(0)
    , m_preorderUpToDate(false) { };

  /**
   * Returns the number of nodes in the tree.
   */
  size_t
  getNodeNumber() const { return m_nodeCounter; };

  /**
   * Returns the number of leaves in the tree.
   */
  size_t
  getLeafNumber() const { return m_leafCounter; };

  /**
   * Adds a child node to another node.
   * @param parent pointer to the parent node.
   * @param child pointer to the child node.
   */
  void
  addChild(std::shared_ptr<TreeNode> parent, std::shared_ptr<TreeNode> child);

  /**
   * Computes the leaves in the subtrees for all internal nodes.
   */
  void
  supplyLeaves(std::shared_ptr<TreeNode> node = nullptr);

  /**
   * Returns a preorder traversal of the tree.
   */
  const std::vector<std::shared_ptr<TreeNode>>&
  getPreorder();

  /**
   * Returns the tree in standard Newick format.
   */
  std::string
  toNewick() const;

  /**
   * Parses a tree from standard Newick format.
   */
  static Tree
  parseNewick(std::string newick);

  /**
   * Returns a list of lists with the species in the subtrees under the root.
   */
  std::vector<std::vector<std::string>>
  speciesInSubtrees() const;

private:
  std::shared_ptr<TreeNode> m_root;                     /*!< root of the tree */
  std::vector<std::shared_ptr<TreeNode>> m_preorder;    /*!< preorder list */
  size_t m_nodeCounter;                                 /*!< number of nodes */
  size_t m_leafCounter;                                 /*!< number of leaves */
  bool m_preorderUpToDate;                              /*!< true if preorder list is up-to-date */

  /**
   * Builds a new preorder list.
   */
  void
  buildPreorder(std::shared_ptr<TreeNode> node = nullptr);

  /**
   * Initializes preorder list and indeces.
   */
  void
  initPreorderAndIndex(std::shared_ptr<TreeNode> node = nullptr);

  /**
   * Recursively constructs the tree in Newick format.
   */
  std::string
  constructNewick(std::shared_ptr<TreeNode> node) const;

  /**
   * Parses a subtree in Newick format.
   */
  static void
  parseSubtree(Tree& tree, std::shared_ptr<TreeNode> subRoot, std::string subtreeString);

  /**
   * Splits the Newick string w.r.t. the children of a node.
   */
  static std::vector<std::string>
  splitChildren(std::string childString);

  /**
   * Necessary for constructing a list of lists with the species in the subtrees under the root.
   */
  void
  fillSubtreeVector(std::shared_ptr<TreeNode> node, std::vector<std::string>& subtree) const;
};

#endif /* TREE_H */
