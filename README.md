# Travelling Salesperson
R scripts to solve* the Travelling salesman problem.

A Travelling Salesman problem asks which is the shortest possible route that visits each city and returns to the origin city from a given list of cities. It is an NP-Hard problem, which means it is computationally difficult. It is a problem easy to state, difficult to solve. TSP finds its application in vehicle routing, PCB/IC Design etc. 

Christofides Algorithm:
	Christofides Algorithm is an approximation algorithm that guarantees that its solution will be within a factor of 3/2 of the optimal solution length. Hence, it is also known as the 3/2 approximation / 1.5 approximation algorithm. It was discovered in 1976 by Nicos Christofides and Anatoliy I. Serdyukov. Firstly, we create a minimum spanning tree graph by using one of the many algorithms, such as Prim’s algorithm, Kruskal’s algorithm, etc. Then, we find all the odd-degree vertices(vertices with odd number of edges connected to them) and match these vertices, by the application of minimum weight perfect matching techniques. Then we combine these two graphs to form an Eulerian cycle (A cycle that visits every given node exactly once). To this day, this algorithm has provided the best approximation ratio that has been proven for the traveling salesman problem on general metric spaces. 
	
*	The input coordinates of the nodes can be provided in the form of a file or rather a random coordinate-generating function employing a certain distribution technique (normalised / uniform random distribution).
*	Read the coordinates of the cities and plot them on a 2D graph.
*	Calculate the Distance Matrix - an upper-diagonal 2D matrix that contains the euclidean distance calculated between all the cities taken in pairs. 
*	Construct the Minimum Spanning Tree (MST) - which is essentially a graph with no cycles and minimum possible weight (here euclidean distance). To implement Kruskal’s algorithm: We consider each vertex to be a seperate tree. We traverse all the node pairs in the graph in the ascending order of distance. If two vertices are from different trees, we join these nodes and combine them into a single tree. On repeating this process for all the nodes, we obtain the MST of the undirected graph of nodes.
*	Find all the odd-degree vertices in the MST.
*	Match all the odd-degree vertices using minimum weight perfect matching techniques. There cannot be an odd number of odd-degree vertices, since each vertex must have 2 edges. Therefore, no vertices will be left out in this process.
*	Construct a circuit by combining the MST and the perfect-matched edges. The resultant graph will have vertices with more than two edges (ie. a single city might be visited more than once). We remedy this by removing repeating vertices.
*	As an additional optimisation step, we remove all the intersections by performing a simple Subtour reversal of nodes. We traverse all the edges of the graph and check for intersections with other edges. If an intersection is found between node x and y, we reverse all the list of nodes between x and y. We repeat this process until all the nodes are visited without correction.
*	Finally, we calculate the total distance covered by the tour and the time taken.
* 	Plot the solution path and save it in the output file.
 
  ![Performance_Christofides_Time](https://github.com/RahulRavi1997/travelling-salesperson/blob/main/images/Performance_Christophides_Time.png)
It is seen that the time computation power of the Christofides algorithm is quite faster and linear.
Whereas on the negative side, the algorithm differs widely with respect to the circuit distance wherein it ranges from 13000 to 30000 which clearly proves that the algorithm only proves the lower bound remains suboptimal.
  ![Performance_Christofides_Distance](https://github.com/RahulRavi1997/travelling-salesperson/blob/main/images/Performance_Christphides_Distance.png)
The choice of algorithm clearly depends on the requirement. At the same time, the computation power of the system used to execute the code also plays a vital role, for example, when I executed the code in a computer with lower processing capabilities, it took me way longer than when I executed it on a system with higher processing capabilities.



Cheapest Insertion Heuristic
	Cheapest Insertion Heuristic is an intuitive approach to the Travelling Salesman Problem. We start with a tour on small subsets of nodes and then extend this tour by inserting the remaining nodes one after the other until all nodes have been inserted. In case of Cheapest Insertion Heuristics, among all nodes inserted so far, we choose a node whose insertion causes the lowest increase in the tour length. This algorithm provides mostly sub-optimal solutions. But this solution can serve as a valid upper bound to further optimizations
One of the biggest advantage of this method is that we do not tend to have the issue of crossover in between nodes and this provides a comparably optimal solution, whereas on the other hand, due to the number of iterations, this algorithm takes reasonably more time to provide the solution for the problem.


*	The input coordinates of the nodes can be provided in the form of a file or rather a random coordinate-generating function employing a certain distribution technique (normalised / uniform random distribution).
*	Read the coordinates of cities and plot them on a 2D graph.
*	Calculate the Distance Matrix - an upper-diagonal 2D matrix that contains the euclidean distance calculated between all the cities taken in pairs. 
*	Construct a Convex Hull - the smallest convex polygon, that encloses all the given points (here cities), with the fewest number of points as its vertices.  Mark the cities in the convex hull as visited.
*	From the remaining cities, using the distance matrix, poll the city which is the closest to the convex hull. This ensures that the increase in weight (here euclidean distance) is small when added to the current set of cities.
*	Now, we find the index of the tour path, at which the selected city must be inserted, so that we obtain the cheapest insertion cost. Insert this city into the set of visited cities and mark it as visited.
*	Repeat Steps 5 & 6 until all the given input cities are visited.
*	Calculate the total distance covered by the tour and the time taken.
* 	Plot the solution path and save it in the output file.


Performance and optimality are very high for lower number of nodes whereas this stores only the data matrix on the memory hence the process is quite simple
 
 ![Performance_Cheapest_Time](https://github.com/RahulRavi1997/travelling-salesperson/blob/main/images/Performance_Cheapest_Time.png)
 
It is seen that the time increment for the cheapest insertion algorithm is quite linear as the number of nodes increases, whereas, due to the complexity of iterations; the algorithm needs higher run time to complete the computation.

 ![Performance_Cheapest_Distance](https://github.com/RahulRavi1997/travelling-salesperson/blob/main/images/Performance_Cheapest_Distance.png)

On the other hand, the distance taken for the circuit to complete remains within a range in between 11000 – 12500. 
Hence it can be clearly inferred that this algorithm has better efficiency in terms of distance but on the other hand has lesser efficiency in terms of computation time.
 
   ![Performance_Comparison_Distance](https://github.com/RahulRavi1997/travelling-salesperson/blob/main/images/Performance_Comparison_Distance.png)
   ![Performance_Comparison_Distance](https://github.com/RahulRavi1997/travelling-salesperson/blob/main/images/Performance_Comparison_Distance.png)


A few input and output testcases are also commited in this repository under the testcases folder, along with the circuit images.
The data provided are subject to the specifications of the computer in which this script is run.

*These do not solve the problem, but provide valid upper bounds and approximation. But the word "solve" does look good doesn't it. 😀️
*All the distances are calculated with Euclidean path 
** Time taken may vary with the speed of your processor.


