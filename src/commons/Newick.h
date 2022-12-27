#include <string>
#include <vector>
#include <map>

struct TreeNode {
    std::string name;
    float branchLength;
    std::vector<TreeNode *> children;
    TreeNode *parent = nullptr;
};

class NewickParser {
    private:
        std::vector<TreeNode> nodes;
        TreeNode tree;
        int treeSize;
        std::map<std::string, int> label;
    public:
};


void parseNewick(std::string newick);