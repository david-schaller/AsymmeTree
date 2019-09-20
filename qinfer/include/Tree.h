#ifndef TREE_H
#define TREE_H

#include <memory>
#include <vector>
#include <string>

class TreeNode {
public:
  friend class Tree;
  TreeNode(std::string value)
    : m_value(value)
    , m_nodeIdx(-1)
    , m_leafIdx(-1) { };

  std::string
  getValue() const { return m_value; };

  size_t
  getNodeIdx() const { return m_nodeIdx; };

  size_t
  getLeafIdx() const { return m_leafIdx; };

  bool
  hasParent() const { return !m_parent.expired(); };

  bool
  hasChildren() const { return !m_children.empty(); };

  const std::vector<std::shared_ptr<TreeNode>>&
  getChildren() const { return m_children; };

  const std::vector<std::shared_ptr<TreeNode>>&
  getLeaves() const { return m_leaves; };

  const std::weak_ptr<TreeNode>&
  getParent() const { return m_parent; };

  std::shared_ptr<TreeNode>
  addChild(std::string value, std::shared_ptr<TreeNode> parent);

private:
  std::string m_value;
  size_t m_nodeIdx;
  size_t m_leafIdx;
  std::vector<std::shared_ptr<TreeNode>> m_children;
  std::weak_ptr<TreeNode> m_parent;
  std::vector<std::shared_ptr<TreeNode>> m_leaves;
};


class Tree {
public:
  Tree()
    : m_root(std::shared_ptr<TreeNode>())
    , m_nodeCounter(0)
    , m_leafCounter(0)
    , m_preorderUpToDate(false) { };

  Tree(std::shared_ptr<TreeNode> root)
    : m_root(root)
    , m_nodeCounter(0)
    , m_leafCounter(0)
    , m_preorderUpToDate(false) { };

  size_t
  getNodeNumber() const { return m_nodeCounter; };

  size_t
  getLeafNumber() const { return m_leafCounter; };

  void
  addChild(std::shared_ptr<TreeNode> parent, std::shared_ptr<TreeNode> child);

  void
  supplyLeaves(std::shared_ptr<TreeNode> node = nullptr);

  const std::vector<std::shared_ptr<TreeNode>>&
  getPreorder();

  std::string
  toNewick() const;

  static Tree
  parseNewick(std::string newick);

  std::vector<std::vector<std::string>>
  speciesInSubtrees() const;

private:
  std::shared_ptr<TreeNode> m_root;
  std::vector<std::shared_ptr<TreeNode>> m_preorder;
  size_t m_nodeCounter;
  size_t m_leafCounter;
  bool m_preorderUpToDate;

  void
  buildPreorder(std::shared_ptr<TreeNode> node = nullptr);

  void
  initPreorderAndIndex(std::shared_ptr<TreeNode> node = nullptr);

  std::string
  constructNewick(std::shared_ptr<TreeNode> node) const;

  static void
  parseSubtree(Tree& tree, std::shared_ptr<TreeNode> subRoot, std::string subtreeString);

  static std::vector<std::string>
  splitChildren(std::string childString);

  void
  fillSubtreeVector(std::shared_ptr<TreeNode> node, std::vector<std::string>& subtree) const;
};

#endif /* TREE_H */
