library(tseries)

# TSP - Christofides Algorithm - 3/2 Approximation Algorithm
# This program tries to find solutions to the famous np-hard 'Travelling Salesman Problem' using Christofides Algorithm.
# This is an approximation algorithm discovered by Nicos Christofides and Anatoliy I. Serdyukov.
# It is the best known approximation algorithm to this day. We can implement this algorithm using the following steps:
# Firstly, we construct the minimum spanning tree.
# Then we select the vertices containing odd number of edges.
# We construct a minimum-weight perfect matching using the odd-degree vertices (There can't be odd number of odd-degree vertices, since each vertex has 2 edges) 
# Now we merge the MST and the perfect matching graphs to obtain a solution.

# working directory
wd = getwd()
initTSP = function() {
  
  startTime = Sys.time()
  print("TSP Christofides algorithm started")
  print(startTime)
  
  # Pre-written test cases
  testcase = 5
  inputfile = paste(wd, "/testcases/input/testcase_" ,testcase, sep="")
  outputfile = paste(wd, "/testcases/output/christofides/result_" ,testcase, sep="")
  nodes = matrix(,,ncol=2)
  nodes <- read.matrix(inputfile)
  rownames(nodes) <- NULL
  # Uncomment below to create nodes in a random fashion
  # numberOfNodes = 10
  # nodes = matrix(c(runif(numberOfNodes*2)*100), ,ncol = 2)
  
  numberOfNodes = nrow(nodes)
  
  # plot(data, ...graphical_parmas)
  # text(x-axis, y-axis, names) but here x-y-coords is used so y defaults to NULL
  plot(nodes, asp = 1)
  colnames(nodes) = c("X", "Y")
  # rownames(nodes) = c(1:numberOfNodes)
  # text(nodes-2,,rownames(nodes))
  
  # array(data, dimensions(x, y, z)) => makes 3d arrays 
  # nrow -> number of rows in the vector
  # sample -> sample(x, size)
  ## array indices start from 1 in R 
  res = solveTsp(nodes)
  optimumOrder = res$path
  distance = res$length
  # plot(nodes, asp=1)
  # # text(nodes-2,,rownames(nodes))
  lines(nodes[optimumOrder,], col= "red")
  
  endTime = Sys.time()
  print("TSP Christofides algorithm completed")
  print(endTime)
  print("Time taken")
  print(endTime - startTime)
  print("Best Tour")
  print(optimumOrder)
  print("Tour Length")
  print(distance)
  write(optimumOrder ,file = outputfile, sep=' -> ',append = FALSE)
}

solveTsp = function(data){
  # build a graph
  G = build_graph(data)
  print("Graph: ")
  print(G)
  
  # build a minimum spanning tree
  MSTree = minimum_spanning_tree(G, data)
  print("MSTree: ")
  print(MSTree)
  
  # find odd vertexes
  odd_vertexes = find_odd_vertexes(MSTree, data)
  print("Odd vertexes in MSTree: ")
  print(odd_vertexes)
  
  # add minimum weight matching edges to MST
  MSTree = minimum_weight_matching(MSTree, G, odd_vertexes)
  print("Minimum weight perfect matching: ")
  print (MSTree)
  
  # find an eulerian tour
  eulerian_tour = find_eulerian_tour(MSTree, G, data)
  
  print("Eulerian tour: ")
  print(eulerian_tour)
  
  current = eulerian_tour[1]
  path = array(c(current),)
  visited = array(rep(FALSE, length(eulerian_tour)),) 
  visited[eulerian_tour[1]] = TRUE
  
  length = 0
  
  for(v in 2:length(eulerian_tour)) {
    if(!visited[eulerian_tour[v]]) {
      path[length(path)+1] = eulerian_tour[v]
      visited[eulerian_tour[v]] = TRUE
      length =  length + G[current, eulerian_tour[v]]
      current = eulerian_tour[v]
    }
  }
  plot(data, asp=1)
  lines(data[c(path, path[1]),])
  # rownames(data) = c(1:nrow(data))
  # text(data-2,,rownames(data))

  # remove intersections using 'Subtour reversal approach'
  path = subtour_reversal(data, path)

  print("Subtour reversal: ")
  print(path)
  length = length + G[path[length(path)], path[1]]
  path[length(path) + 1] = path[1]
  return(list("path"=path, "length"=length))
}

get_length = function(x1, y1, x2, y2){
  return ((((x1 - x2) ^ 2) +((y1 - y2) ^ 2)))
}


build_graph = function(data) {
  graph = array(, c(nrow(data), nrow(data)))
  for(this in 1: nrow(data)) {
    for(another_point in 1: nrow(data)) {
      if(this != another_point){
        graph[this, another_point] = get_length(data[this, 1], data[this, 2], data[another_point, 1], data[another_point, 2])
      }
    }
  }
  return(graph)
}

parents = rep(NA, 1)
weights = 0
getitem = function(object) {
  if(is.na(parents[object]) || (parents[object] == 0)) {
    parents[object] <<- object
    weights[object] <<- 1
    return(object)
  }
  # find path of objects leading to the root
  path = array(c(object))
  root = parents[object]
  while(root != path[length(path)]) {
    path[length(path)+1] = root
    root = parents[root]
  }
  for(ancestor in path) {
    parents[ancestor] <<- root
  }
  return(root)
}

union = function(u, v){
  roots = c(getitem(u), getitem(v))
  arr = c(weights[roots[1]], weights[roots[2]])
  heaviest = roots[which.max(arr)]
  r = roots[2]
  if(roots[1] != heaviest) {
    r = roots[1]
  }
  weights[heaviest] <<- weights[heaviest]  + weights[r]
  parents[r] <<- heaviest
}

minimum_spanning_tree = function(G, nodes){
  n = nrow(nodes)
  tree = numeric()
  distances = array(, c(n*n, 3))
  dist = 0
  index = 1
  for (u in 1: n) {
    for (v in 1: n) {
      distances[index,] = c(G[u, v], u, v)
      dist[index] = G[u, v]
      index = index + 1
    }
  }
  order = order(dist, decreasing=FALSE)
  sortedDistances = distances[order, ]
  sortedDistances = matrix(sortedDistances[!is.na(sortedDistances[,1])] ,, ncol=3)
  for(el in 1:nrow(sortedDistances)) {
    W = sortedDistances[el,1]
    u = sortedDistances[el,2]
    v = sortedDistances[el,3]
    if (is.na(W)) {
      next
    }
    if(getitem(u) != getitem(v)) {
      tree = rbind(tree, c(u, v, W))
      union(u, v)
    }
  }
  return(tree)
}


find_odd_vertexes = function(MST, nodes){
  n = nrow(nodes)
  tmp_g = array(rep(0,n),c(n))
  vertexes = array()
  for(i in 1:nrow(MST)) {
    edge = MST[i,]
    tmp_g[edge[1]] = tmp_g[edge[1]] +1
    tmp_g[edge[2]] = tmp_g[edge[2]] +1
  }
  
  k = 1;
  for(vertex in 1:length(tmp_g)) {
    if(tmp_g[vertex] %% 2 == 1) {
      vertexes[k] = (vertex)  
      k = k + 1
    }
  }
  return(vertexes)
}


minimum_weight_matching = function(MST, G, odd_vert){
  # odd_vert = sample(odd_vert, length(odd_vert))
  while(length(odd_vert) != 0) {
    v = odd_vert[length(odd_vert)]
    odd_vert = odd_vert[-length(odd_vert)]
    length = Inf
    closest = 0
    for(u in 1:length(odd_vert)) {
      vert = odd_vert[u]
      if ((v != vert) && (G[v, vert] < length)){
        length = G[v, vert]
        closest = vert
      }
    }
    MST = rbind(MST, c(v, closest, length))
    odd_vert = odd_vert[!(odd_vert %in% closest)]
  }
  return(MST)
}

find_eulerian_tour = function(MatchedMSTree, G, nodes) {
  # find neighbours
  n = nrow(nodes)
  neighbours = as.list(rep(0, n))
  for(i in 1:nrow(MatchedMSTree)) {
    edge = MatchedMSTree[i,]
    neighbours[[edge[1]]] <- c(neighbours[[edge[1]]], c(edge[2]))
    neighbours[[edge[2]]] <- c(neighbours[[edge[2]]], c(edge[1]))
  }
  
  # finds the hamiltonian circuit
  start_vertex = MatchedMSTree[1, 1]
  EP = neighbours[[start_vertex]][2]
  
  while(length(MatchedMSTree) > 0) {
    i = 1
    for (s in 1:length(EP)) {
      v = EP[s];
      if(length(neighbours[[v]]) > 1){        
        break
      }
      i = i + 1
    }
    while (length(neighbours[[v]]) > 1) {
      w = neighbours[[v]][2]
      MatchedMSTree = matrix(remove_edge_from_matchedMST(MatchedMSTree, v, w),,3)
      
      neighbours[[v]] <- c(na.omit(neighbours[[v]][-match(w, neighbours[[v]])]))
      neighbours[[w]] <- c(na.omit(neighbours[[w]][-match(v, neighbours[[w]])]))
      EP = append(EP, w, after=i)
      i = i + 1
      v = w
    }
  }
  return(EP)
}


remove_edge_from_matchedMST = function(MatchedMST, v1, v2){
  toRemoveEdges = vector()
  it = 1
  while (it <= nrow(MatchedMST)) {
    edge1 = MatchedMST[it, 1];
    edge2 = MatchedMST[it, 2];
    if (((edge1 == v2) && (edge2 == v1)) || ((edge1 == v1) && (edge2 == v2))) {
      toRemoveEdges = c(toRemoveEdges, it)
    }
    it = it + 1
  }
  if (length(toRemoveEdges) > 0) {
    MatchedMST = MatchedMST[-toRemoveEdges,]
  }
  return(MatchedMST)
}

subtour_reversal = function(nodes, path) {
  n = nrow(nodes)
  k = 0
  iter = 1
  while(iter <= n) {
    k = (k %% (n-1)) + 1
    crossNode = areIntersectingNodes(nodes, path, c(path[k], path[k+1]))
    if (crossNode!= -1) {
      newPath = path
      start = min((k), crossNode)
      end = max((k), crossNode)
      if (start != n) {
        newPath[(start+1):end] = path[end:(start+1)]
      } else {
        newPath[1:end] = path[end:1]
      }
      path = newPath
      plot(nodes, asp=1)
      lines(nodes[c(path, path[1]),])
      # rownames(nodes) = c(1:nrow(nodes))
      # text(nodes-2,,rownames(nodes))
      iter = 1
    }
    iter = iter + 1
  }
  return(path);
}

areIntersectingNodes = function(nodes, nodeOrder, checkNodes) {
  nodeOrder = c(nodeOrder, nodeOrder[1])
    for (i in (1:(length(nodeOrder)-1))) {
      if (checkIntersecting(nodes, checkNodes[1], checkNodes[2], nodeOrder[i], nodeOrder[i + 1])) {
        return(i);
      }
    }
  return(-1)
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
