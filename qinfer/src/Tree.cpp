#include <iostream>
#include <stdexcept>
#include <vector>

#include "Tree.h"

void
Tree::addChild(std::shared_ptr<TreeNode> parent, std::shared_ptr<TreeNode> child){
  ++m_nodeCounter;
  if(parent->hasChildren()){
    ++m_leafCounter;
  }
  parent->m_children.push_back(std::shared_ptr<TreeNode>(child));
  child->m_parent = parent;

  m_preorderUpToDate = false;
}

void
Tree::buildPreorder(std::shared_ptr<TreeNode> node){

  if(node == nullptr){
    m_preorder.clear();
    buildPreorder(m_root);
    m_preorderUpToDate = true;

  } else {
    m_preorder.push_back(node);
    for(auto child : node->getChildren()){
      buildPreorder(child);
    }
  }
}

void
Tree::initPreorderAndIndex(std::shared_ptr<TreeNode> node){

  if(node == nullptr){
    m_preorder.clear();
    m_nodeCounter = 0;
    m_leafCounter = 0;
    initPreorderAndIndex(m_root);
    m_preorderUpToDate = true;

  } else {
    m_preorder.push_back(node);
    node->m_nodeIdx = m_nodeCounter;
    ++m_nodeCounter;
    if(!node->hasChildren()){
      node->m_leafIdx = m_leafCounter;
      ++m_leafCounter;
    }

    for(auto child : node->getChildren()){
      initPreorderAndIndex(child);
    }
  }
}

const std::vector<std::shared_ptr<TreeNode>>&
Tree::getPreorder(){
  if(!m_preorderUpToDate){
    buildPreorder();
  }
  return m_preorder;
}

void
Tree::supplyLeaves(std::shared_ptr<TreeNode> node){
  if(node == nullptr){
    supplyLeaves(m_root);
  } else if(!node->hasChildren()){
    node->m_leaves = {node};
  } else {
    node->m_leaves.clear();
    for(auto child : node->getChildren()){
      supplyLeaves(child);
      node->m_leaves.insert(node->m_leaves.end(),
                            child->m_leaves.begin(),
                            child->m_leaves.end());
    }
  }
}

std::string
Tree::toNewick() const {
  return constructNewick(m_root) + ";";
}

std::string
Tree::constructNewick(std::shared_ptr<TreeNode> node) const {
  if(!node->hasChildren()){
    return node->getValue();
  } else {
    std::string s = "(";
    for(auto child : node->getChildren()){
      s += constructNewick(child) + ",";
    }
    return s.substr(0, s.size()-1) + ")" + node->getValue();
  }
}

Tree
Tree::parseNewick(std::string newick){

  if(newick.find(';') != std::string::npos){
    newick.erase(newick.find(';'), newick.size());
  }

  auto tempRoot = std::make_shared<TreeNode>(std::basic_string("-1"));
  auto tree = Tree(std::shared_ptr<TreeNode>());

  parseSubtree(tree, tempRoot, newick);

  if(tempRoot->hasChildren()){
    // remove the temporary root
    tree.m_root = tempRoot->m_children[0];
    tempRoot->m_children.clear();
    tree.m_root->m_parent.reset();
  } else {
    throw std::invalid_argument( "invalid Newick string! (empty)" );
  }

  tree.initPreorderAndIndex();

  return tree;
}

void
Tree::parseSubtree(Tree& tree,
                   std::shared_ptr<TreeNode> subRoot,
                   std::string subtreeString){
  auto children = splitChildren(subtreeString);

  for(std::string& child : children){

    auto node = std::make_shared<TreeNode>(std::string(""));
    tree.addChild(subRoot, node);

    size_t subtreeEnd = std::string::npos;
    if(child.size() > 1 && child.at(0) == '('){
      subtreeEnd = child.rfind(')');
      if(subtreeEnd == std::string::npos){
        throw std::invalid_argument( "invalid Newick string! (brackets 1)" );
      } else {
        parseSubtree(tree, node, child.substr(1, subtreeEnd-1));
        child.erase(0, subtreeEnd+1);
      }
    }

    // assign value to node (trim distance information)
    if(child.find(":") != std::string::npos){
      node->m_value = child.substr(0, child.find(":"));
    } else {
      node->m_value = child;
    }
  }
}

std::vector<std::string>
Tree::splitChildren(std::string childString){
  auto children = std::vector<std::string>();
  int stack = 0;
  std::string current = "";

  for(char c : childString){
    if((stack == 0) && (c == ',')){
      children.push_back(current);
      current.clear();
    } else if(c == '('){
      ++stack;
      current.append(1, c);
    } else if(c == ')'){
      if(stack <= 0){
        throw std::invalid_argument( "invalid Newick string! (brackets 2)" );
      }
      --stack;
      current.append(1, c);
    } else {
      current.append(1, c);
    }
  }
  children.push_back(current);

  return children;
}

std::vector<std::vector<std::string>>
Tree::speciesInSubtrees() const {
  auto subtrees = std::vector<std::vector<std::string>>();

  // walk to the true root i.e. the first with multiple children (rho_S)
  std::shared_ptr<TreeNode> pos(m_root);
  while (pos->m_children.size() == 1) {
    pos = pos->m_children[0];
  }

  // fill the lists
  for(auto& child : pos->m_children){
    subtrees.push_back(std::vector<std::string>());
    fillSubtreeVector(child, subtrees.back());
  }

  return subtrees;
}

void
Tree::fillSubtreeVector(std::shared_ptr<TreeNode> node,
                        std::vector<std::string>& subtree) const {
  if(!node->hasChildren()){
    subtree.push_back(node->getValue());
  } else {
    for(auto& child : node->m_children){
      fillSubtreeVector(child, subtree);
    }
  }
}
