library(tseries)

# TSP - Cheapest Insertion Heursistic
# This program tries to find solutions to the famous np-hard 'Travelling Salesman Problem' using Cheapest Insertion Heuristic.
# This is a greedy approach to solving this problem, in which the cheapest(least distance) node is visited first.
# Firstly, we create a convex hull. Then we select the nearest point to the convex hull, and add the point to the visited nodes.
# We repeat this process until all nodes are visited.
# This provides sub-optimal solutions in most cases. But this solution can serve as a valid upper bound to further optimsations.

# <-- R Syntax --->
# set.seed(x) # used to randomise data based on X
# c -> combine values into vectors
# matrix(data, nrow=, ncol=) # matrix initialisation
# array(data, c(nrow, ncol?, ndepth?)) # array initialisation 
# runif/rnorm -> random function generates numbers with uniform/normal distribution

# working directory
wd = getwd()
initTSP = function() {
  
  startTime = Sys.time()
  print("TSP cheapest Insertion heuristic started")
  print(startTime)

  # Pre-written test cases
  testcase = 10
  inputfile = paste(wd, "/testcases/input/testcase_" ,testcase, sep="")
  outputfile = paste(wd, "/testcases/output/cheapest/result_" ,testcase, sep="")
  nodes = matrix(,,ncol=2)
  nodes <- read.matrix(inputfile)
  rownames(nodes) <- NULL
  # Uncomment below to create nodes in a random fashion
  # numberOfNodes = 50
  # nodes = matrix(c(runif(numberOfNodes*2)*100), ,ncol = 2)
  
  numberOfNodes = nrow(nodes)
  
  # plot(data, ...graphical_parmas)
  # text(x-axis, y-axis, names) but here x-y-coords is used so y defaults to NULL
  plot(nodes, asp = 1)
  colnames(nodes) = c("X", "Y")
  #rownames(nodes) = c(1:numberOfNodes)
  #text(nodes-2,,rownames(nodes))
  
  # array(data, dimensions(x, y, z)) => makes 3d arrays 
  # nrow -> number of rows in the vector
  # sample -> sample(x, size)
  ## array indices start from 1 in R 
  optimumOrder = solveTsp(nodes)
  pathdistances = array(, numberOfNodes-1)
  pathdistances = calculatePathDistances(nodes, optimumOrder)
  distance = sum(pathdistances)
  plot(nodes, asp=1)
  text(nodes-2,,rownames(nodes))
  lines(getNodesFromNodeOrder(nodes, c(optimumOrder, optimumOrder[1])), col= "red")
  
  endTime = Sys.time()
  print("TSP cheapest Insertion heuristic completed")
  print(endTime)
  print("Time taken")
  print(endTime - startTime)
  print("Best Tour")
  print(optimumOrder)
  print("Tour Length")
  print(distance)
  write(optimumOrder ,file = outputfile, sep=' -> ',append = FALSE)
}

solveTsp = function(nodes) {
  distanceMatrix = constructDistanceMatrix(nodes)
  optimumOrder = chull(nodes)
  optimumOrder = constructOptimumOrder(nodes, distanceMatrix, optimumOrder)
  return(optimumOrder)
}
constructDistanceMatrix = function (nodes) {
  numberOfNodes = nrow(nodes)
  dist = matrix(,numberOfNodes, numberOfNodes)
  row = 1
  while (row <= numberOfNodes) {
    column = row
    while (column <= numberOfNodes) {
      if (row == column) {
        dist[row,column] = 0
      } else {
        dis = getDistance(nodes, row, column);
        dist[row,column] = dis
        dist[column,row] = dis
      }
      column = column + 1
    }
    row = row + 1
  }
  return(dist)
}

getNodesFromNodeOrder = function(nodes, nodeOrder) {
  numberOfNodes = nrow(nodes)
  plotNodes = array(, c(numberOfNodes+1, 2));
  for (i in 1:(numberOfNodes+1)) {
    plotNodes[i,] = nodes[nodeOrder[i],]
  }
  return(plotNodes)
}

#Function to calculate distances between nodes => √((x2-x1)²+(y2-y1)²) // leaving sqrt function for optimsation
getDistance = function(nodes, firstPoint, secondPoint) {
  dist = (nodes[secondPoint, 1]-nodes[firstPoint,1])^2 + (nodes[secondPoint, 2]-nodes[firstPoint, 2])^2
  return(dist)
}

#Function to calculate path distance
calculatePathDistances = function(nodes, points) {
  distance = 0
  numberOfNodes = nrow(nodes)
  for(i in 1:(numberOfNodes-1)) {
    distance[i] = getDistance(nodes, points[i], points[i+1])
  }
  return(distance)
}
#Function to calculate path distance
calculateFullPathDistance = function(nodes, points) {
  distance = 0
  i = 0
  while (i <(length(points))) {
    distance[i] = getDistance(nodes, points[i], points[i+1])
    i = i + 1
  }
  return(sum(distance))
}

constructOptimumOrder = function(nodes, distanceMatrix, startNodes) {
  nodeOrder = 0
  nodeOrder = startNodes
  shortestOrder = nodeOrder
  while (length(nodeOrder) < nrow(nodes)) {
    shortestDistance = Inf
    currentNode = nodeOrder[length(nodeOrder)]
    closestNode = getOptimalIndex(nodes, distanceMatrix, nodeOrder)
    arr = c(distanceMatrix[nodeOrder,closestNode])
    limit = min(c(length(nodeOrder), 11))
    mins = sort(arr, partial=limit)[limit]
    minNodes = which(arr <= (mins))[-1]
    # fetch 10 closest nodes in nodeOrder to find optimal insertion point
    for (i in 1: (length(minNodes))) {
      if (minNodes[i] != 1) {
        start = nodeOrder[1:(minNodes[i]-1)]
        order = c(nodeOrder[minNodes[i]:length(nodeOrder)], start, closestNode);
        currentDistance = calculateFullPathDistance(nodes, c(order, order[1]))
        res = isShorterPath(nodes, currentDistance, order, shortestDistance, nodeOrder)
        if (res != 0) {
          shortestDistance = res
          shortestOrder = order
        }
      }
      if (minNodes[i] != length(nodeOrder)) {
        end = nodeOrder[(minNodes[i]+1):length(nodeOrder)]
        order = c(end, nodeOrder[1:minNodes[i]], closestNode);
        currentDistance = calculateFullPathDistance(nodes, c(order, order[1]))
        res = isShorterPath(nodes, currentDistance, order, shortestDistance, nodeOrder)
        if (res != 0) {
          shortestDistance = res
          shortestOrder = order
        }
      }
    }
    nodeOrder = shortestOrder
    # plot the current optimal solution
    plot(nodes, asp=1)
    text(nodes-2,,rownames(nodes))
    lines(getNodesFromNodeOrder(nodes, c(nodeOrder, nodeOrder[1])), col="red")
    Sys.sleep(0.0001)
  }
  return(nodeOrder)
}


getOptimalIndex = function(nodes, distanceMatrix, currentNodeOrder) {
  allNodes = c(1:nrow(nodes))
  remaining = allNodes[!(allNodes %in% currentNodeOrder)]
  currentNodeOrderLen = length(currentNodeOrder)
  if ((currentNodeOrderLen %% 4) == 1) {
    # Return min X-axis node => Closest node to left of convex hull
    remainingNodes = nodes[remaining, 1]
    return(remaining[which.min(remainingNodes)])
  } else if ((currentNodeOrderLen %% 4) == 2) {
    # Return max X-axis node => Closest node to right of convex hull
    remainingNodes = nodes[remaining, 1]
    return(remaining[which.max(remainingNodes)])
  } else if ((currentNodeOrderLen %% 4) == 3) {
    # Return min Y-axis node => Closest node to top of convex hull
    remainingNodes = nodes[remaining,2]
    return(remaining[which.min(remainingNodes)])
  } else {
    # Return max Y-axis node => Closest node to bottom of convex hull
    remainingNodes = nodes[remaining, 2]
    return(remaining[which.max(remainingNodes)])
  }
}

isShorterPath = function(nodes, currentDistance, order, shortestDistance, optimumOrder) {
  # check if current distance provides better solutuion
  if (currentDistance < shortestDistance) {
    intersectCheckNodes = array(, c(2,2))
    intersectCheckNodes[,1] = c(order[length(order)-1], order[length(order)])
    intersectCheckNodes[,2] = c(order[length(order)], order[1])
    # check for intersections
    # By Triangle inequality theorem, having intersections will always reduce optimality
    if (length(optimumOrder) < 3 || !areIntersectingNodes(nodes, order, intersectCheckNodes)) {
      shortestDistance = currentDistance
      shortestOrder = order
      return(shortestDistance)
    }
  }
  return(0)
}

areIntersectingNodes = function(nodes, nodeOrder, checkNodes) {
  for (j in 1:2) {
    for (i in 1:(min(5:(length(nodeOrder)-1)))) {
      if (checkIntersecting(nodes, checkNodes[j,1], checkNodes[j,2], nodeOrder[i], nodeOrder[i + 1])) {
        return(TRUE);
      }
    }
    for (i in (length(nodeOrder)-1):(max(1:(length(nodeOrder)-6)))) {
      if (checkIntersecting(nodes, checkNodes[j,1], checkNodes[j,2], nodeOrder[i], nodeOrder[i + 1])) {
        return(TRUE);
      }
    }
  }
  return(FALSE)
}

checkIntersecting = function(nodes, startNode, endNode, startCompareNode, endCompareNode) {
  if (startCompareNode != startNode && startCompareNode != endNode && endCompareNode != startNode && endCompareNode!= endNode) {
    return(isIntersecting(nodes[startNode,], nodes[endNode,], nodes[startCompareNode,], nodes[endCompareNode,]))
  }
  return(FALSE)
}

isIntersecting = function(p1, q1, p2, q2) {
  o1 = getOrientation(p1, q1, p2); 
  o2 = getOrientation(p1, q1, q2); 
  o3 = getOrientation(p2, q2, p1); 
  o4 = getOrientation(p2, q2, q1); 
  if (o1 != o2 && o3 != o4) 
    return(TRUE); 
  if (o1 == 0 && isOnLineSegment(p1, p2, q1)) return(TRUE); 
  if (o2 == 0 && isOnLineSegment(p1, q2, q1)) return(TRUE); 
  if (o3 == 0 && isOnLineSegment(p2, p1, q2)) return(TRUE); 
  if (o4 == 0 && isOnLineSegment(p2, q1, q2)) return(TRUE); 
  return(FALSE);
}

getOrientation = function(p, q, r) {
  px = p[1]
  qx = q[1]
  rx = r[1]
  py = p[2]
  qy = q[2]
  ry = r[2]
  val = (qy - py) * (rx - qx) - (qx - px) * (ry - qy);
  val = val[1] 
  if (val == 0) return(0); 
  if ((val > 0)) {
    return(1)
  }
  return(2); 
}
isOnLineSegment = function(p, q, r){
  px = p[1]
  qx = q[1]
  rx = r[2]
  py = p[2]
  qy = q[2]
  ry = r[2]
  if (qx <= max(px, rx) && qx >= min(px, rx) &&  qy <= max(py, ry) && qy >= min(py, ry)) 
    return(TRUE);
  return(FALSE); 
}

# Utility Function to swap 2 points
swap = function(a, i, j) {
  temp = a[i]
  a[i] = a[j]
  a[j] = temp
  return(a)
}

# Starts the TSP algorithm
initTSP()

# References
# https://en.wikipedia.org/wiki/Travelling_salesman_problem - Wikipedia article - Travelling salesman problem
