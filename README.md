# UMApHCP
This project implements metaheuristics to solve the Uncapacitated Multiple Allocation p-hub Center Problem (UMApHCP). The problem involves the efficient distribution of products among nodes, minimizing transportation costs.

# UMApHCP Problem
UMApHCP is a common problem in distribution systems planning, such as telecommunications and transportation. It consists of determining which nodes will be hubs and calculating the transportation cost of all node pairs, passing through at least one hub, to minimize the highest transportation cost among all pairs.

# Key Features
Input Processing: Reads node coordinates from an input file (inst200.txt).  
Distance Calculation: Computes Euclidean distances between all nodes.  
Hub Selection: Selects hubs based on total distance metrics using a sorting approach.  
Route Optimization: Determines optimal routes between nodes via hubs, minimizing transportation costs.  
Solution Construction: Builds an initial solution with selected hubs and assigned routes.  
Solution Management: Supports saving and loading solutions for further analysis.  
Performance Monitoring: Measures execution time for key operations (hub selection, cost calculation).  
Memory Management: Handles dynamic allocation and deallocation of resources to prevent leaks.

# Important Note
Due to the fact that this is an academic project, I am unable to provide some of the necessary files for you to run the program locally. However, the results have been made available.
