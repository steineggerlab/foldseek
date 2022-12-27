#include <iostream>
#include <regex>
#include "Newick.h"


void parseNewick(std::string newick) {
    std::cout << "Tree:\n\t" << newick << std::endl; 

    std::vector<TreeNode> nodes;
    std::vector<TreeNode *> ancestors;
    
    nodes.emplace_back();
    TreeNode *node = &nodes.back();
    TreeNode *subtree = nullptr;

    std::string pattern("\\s*(;|\\(|\\)|,|:)\\s*");
    std::regex regex(pattern);
    std::sregex_token_iterator iter(newick.begin(), newick.end(), regex, { -1, 0 } );
    std::sregex_token_iterator end;
    std::string prevToken;
    
    while (iter != end) {
        std::string token = *iter++;
        
        if (token == "")
            continue;
       
        if (token == "(") {
            nodes.emplace_back();
            subtree = &nodes.back();
            node->children.push_back(subtree);
            ancestors.push_back(node);
            node = subtree;
            
        } else if (token == ",") {
            nodes.emplace_back();
            subtree = &nodes.back();
            ancestors.back()->children.push_back(subtree);
            node = subtree;

        } else if (token == ")") {
            node = ancestors.back();
            ancestors.pop_back();

        } else if (token == ":") {

        } else {
            if (prevToken == "")
                continue;
            if (prevToken == ")" || prevToken == "(" || prevToken == ",") {
                node->name = token;
            } else {
                node->branchLength = std::stof(token);
            }
        }
        
        prevToken = token;
    }
    
    std::cout << "Ancestors: " << ancestors.size() << std::endl;
    for (size_t i = 0; i < nodes.size(); i++) {
        std::cout << "Node: " << nodes[i].name << ", children: ";
        for (size_t j = 0; j < nodes[i].children.size(); j++) {
            std::cout << nodes[i].children[j]->name << ", ";
        }
        std::cout << std::endl;
    }

}