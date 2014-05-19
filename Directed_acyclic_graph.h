// Author: Christopher Mitcheltree
// Date: Fall 2013

#ifndef DIRECTED_ACYCLIC_GRAPH_H
#define DIRECTED_ACYCLIC_GRAPH_H

#include "Double_hash_table.h"
#include <limits>
#include <iostream>

class Directed_acyclic_graph {
	private:
		bool resetPriority; // Keeps a flag to whether any priority value has been changed
		int edgeCount; // Keeps count of how many edge connections there are
		int vertices; // Contains the number of vertices
		int vIn; // Keeps the starting position of in-degrees for vertecies
		int vOut; // Keeps the starting position of out-degrees for vertecies
		int vToDo; // Starting index of the todo queue for connected function
		int vDone; // Starting index of the done queue for 
		int vQueue; // Starting index of the sort queue
		int vIn2; // Starting index for changing in-degree in the sort
		int length; // Length of the dag array 
		int *dag; // Array which stores all the above information
		double *priority; // Array which stores vertice's priority
		int hashSize; // The power of 2 of the hash size array
		Double_hash_table<double> *hashTable; // hash table type

		void insertQueue(int, int, int) const;

	public:
		Directed_acyclic_graph( int = 10 );
		Directed_acyclic_graph( Directed_acyclic_graph const & );
		~Directed_acyclic_graph();

		void swap( Directed_acyclic_graph & );
		Directed_acyclic_graph &operator = ( Directed_acyclic_graph );

		int in_degree( int ) const;
		int out_degree( int ) const;
		int edge_count() const;
		bool adjacent( int, int ) const;
		bool connected( int, int ) const;
		void topological_sort() const;

		bool set_priority( int, double );
		bool insert_edge( int, int );
		void clear_edges();
		void reset_priorities();

	friend std::ostream &operator << ( std::ostream &, Directed_acyclic_graph const & );
};

// Constructor: creates a new instance of the DAG data structure.
Directed_acyclic_graph::Directed_acyclic_graph(int n):
	resetPriority(false),
	edgeCount(0),
	vertices(n),
	vIn(n*n),
	vOut(vIn + n),
	vToDo(vOut + n), 
	vDone(vToDo + n), 
	vQueue(vDone + n), 
	vIn2(vQueue + n), 
	length(vIn2 + n), 
	dag(new int[length]), 
	priority(new double[vertices]),
	hashSize((int)ceil((log10(double(n)))/(log10(2.0)))),
	hashTable(new Double_hash_table<double>(hashSize)) {
	// Set all vertice's priority to its value
	for(int k = 0; k < vertices; ++k) {
		priority[k] = k;
		hashTable->insert(k);
	}
	// Clear all adjacency, in-deg, out-deg values
	for(int k = 0; k < vToDo; ++k) {
		dag[k] = 0;
	}
}

// Copy Constructor: creates copy of the DAG data structure.
Directed_acyclic_graph::Directed_acyclic_graph(Directed_acyclic_graph const &obj): 
	resetPriority(obj.resetPriority),
	edgeCount(obj.edgeCount),
	vertices(obj.vertices),
	vIn(obj.vIn),
	vOut(obj.vOut),
	vToDo(obj.vToDo),
	vDone(obj.vDone),
	vQueue(obj.vQueue),
	vIn2(obj.vIn2),
	length(obj.length),
	dag(new int[obj.length]), 
	priority(new double[obj.vertices]),
	hashSize(obj.hashSize),
	hashTable(new Double_hash_table<double>(*obj.hashTable)) {
		// Copies over the adjecency relation, in-deg, and out-deg values
		for(int k = 0; k < length; ++k) {
			dag[k] = obj.dag[k];
		}
		// Copies over the priority values
		for(int k = 0; k < vertices; ++k) {
			priority[k] = obj.priority[k];
		}
}

// Deconstructor: Deallocates any memory that was allocated by the DAG constructor.
Directed_acyclic_graph::~Directed_acyclic_graph() {
	delete hashTable;
	delete[] priority;
	delete[] dag;
}

// In Degree Acessor: Returns the in-degree of vertex i (how many edges are entering it)
int Directed_acyclic_graph::in_degree(int i) const {
	// Checks if the vertex value is within bounds
	if(i > vertices - 1) {
		throw illegal_argument();
	}
	if(i < 0) {
		throw illegal_argument();
	}
	return dag[vIn + i];
}

// Out Degree Acessor: Returns the out-degree of vertex i (how many edges are leaving it)
int Directed_acyclic_graph::out_degree(int i) const {
	// Checks if the vertex value is within bounds
	if(i > vertices - 1) {
		throw illegal_argument();
	}
	if(i < 0) {
		throw illegal_argument();
	}
	return dag[vOut + i];
}

// Edge Count Acessor: Returns the number of edges stored
// in the current acyclic graph
int Directed_acyclic_graph::edge_count() const {
	return edgeCount;
}

// Adjacent Acessor: Returns a boolean value whether or not vertex
// i and j are adjecent to each other or not. If i and j are the same 
// vertex, returns false. 
bool Directed_acyclic_graph::adjacent(int i, int j) const {
	if(i > (vertices - 1) || j > (vertices - 1)) { // Checks bounds
		throw illegal_argument();
	}
	if(i < 0 || j < 0) { // Checks bounds
		throw illegal_argument();
	}
	if(dag[vertices*i + j] == 1) { 
		return true;
	}
	return false;
}

// Connected Acessor: Returns true if there exists a path from vertex i to vertex j, otherwise returns false.
bool Directed_acyclic_graph::connected(int i, int j) const {
	if(adjacent(i, j)) { // Check for trivial case of adjacency.
		return true;
	}
	if(i == j) { // A vertex is connected to itself, therefore return true.
		return true;
	}
	if(dag[vOut + i] == 0) { // Check for trivial case if out degree of i is 0.  
		return false;
	}
	if(dag[vIn + j] == 0) { // Check trivial case if in degree of j is 0.
		return false;
	}

	for(int k = 0; k < vertices; ++k) { // Reset done vertices array to 0.
		dag[vDone + k] = 0;
	}
	
	int toDo = (i == 0) ? 1 : i; // Sum of vertices left in queue (vertex 0 must have a weight of 1).
	int indexStart = 0; // Start index of queue.
	dag[vToDo] = i;
	int indexEnd = 1; // End index of queue.
	int element;
	int elementOut; // Number of out degrees of element.
	int destination; // Possible destination vertex.

	while(toDo > 0) { // While queue is not empty.
		element = dag[vToDo + indexStart]; // Pop first element off queue.
		elementOut = dag[vOut + element];
		toDo = (element == 0) ? toDo - 1 : toDo - element; // Decrease queue sum accordingly (vertex 0 must have a weight of 1).
		dag[vDone + element] = 1; // Mark current vertex as done.
		destination = 0;

		for(int k = 0; k < elementOut;) { // Iterate through possible destination vertices of current vertex.
			if(dag[vertices * element + destination] == 0) { // Vertex is not a destination.
				++destination; 
			} else { // Vertex is a destination.
				if(destination == j) { // If path has been found, return true.
					return true;
				}
				if(dag[vDone + destination] == 0) {
					dag[vToDo + indexEnd] = destination; // Add vertex to queue if it has not been checked already.
					++indexEnd;
					toDo = (destination == 0) ? toDo + 1 : toDo + destination; // Increase queue sum accordingly (vertex 0 must have a weight of 1).
				}
				++destination;
				++k;
			}
		}
		++indexStart; // Go to next vertex in queue.
	}
	return false; // Return false if no path was found.
}

// Topological Sort Acessor: Prints a topological sort based on vertex priority of the vertices in the DAG (lowest priority goes first).
void Directed_acyclic_graph::topological_sort() const {
	if((edge_count() == 0) && !resetPriority) { // If no edges have been added, print trivial case.
		std::cout << "0";
		for(int k = 1; k < vertices; ++k) {
			std::cout << "-" << k;
		}
	} else {
		for(int k = 0; k < vertices; ++k) { // Set modifiable indegree array to indegree of each vertex.
			dag[vIn2 + k] = dag[vIn + k];
		}
		
		int start = 0; // Start of sort queue.
		int end = 0; // End of sort queue.
		int count = 0; // Count of number of vertices printed.
		
		for(int k = 0; k < vertices; ++k) { // Push vertices with indegree 0 onto the queue.
			if(dag[vIn2 + k] == 0) {
				insertQueue(k, start, end);
				++end;
			}
		}

		while((end - start) > 0) { // While queue is not empty.
			int element = dag[vQueue + start]; // Pop off first vertex from queue.
			int edgesOut = dag[vOut + element]; // Get outdegree of current vertex.
			int indexOut = 0; // Index of possible destination vertices of current vertex.

			std::cout << element;
			++start; // Pop off current vertex.

			for(int k = 0; k < edgesOut;) { // Iterate through possible destination vertices of current vertex.
				if(dag[vertices * element + indexOut] == 0) { // Vertex is not a destination.
					++indexOut;
				} else {
					dag[vIn2 + indexOut] = dag[vIn2 + indexOut] - 1; // Decrease indegree of destination vertices.
					
					if(dag[vIn2 + indexOut] == 0) { // If indegree reaches 0, pop the vertex onto the queue.
						insertQueue(indexOut, start, end);
						++end;
					}					
					++indexOut;
					++k;
				}
			}
			++count;
			if(count != vertices) { // Prevent dash from being printed after last vertex.
				std::cout << "-";
			}
		}
	}
}

// Insert Queue Mutator: Inserts a vertex into the sort queue in its correct position based on priority.
void Directed_acyclic_graph::insertQueue(int i, int start, int end) const{
	if(start == end) { // If sort queue is empty, simply insert at the front.
		dag[vQueue + start] = i;
		return;
	}
	
	for(int k = start; k < end; ++k) { // Iterate through queue.
		if(priority[i] < priority[dag[vQueue + k]]) { // Found position vertex should be inserted.
			for(int m = 0; m < (end - k); ++m) { // Shift all remaining vertices over on queue.
				dag[vQueue + end - m] = dag[vQueue + end - 1 - m]; 
			}
			dag[vQueue + k] = i; // Insert vertex into queue in correct position.

			return;
		}
	}	
	dag[vQueue + end] = i; // Else insert vertex at end of queue.
}

// Set Priority Mutator: Checks if priority value exists, if not, update specified 
// vertex with the priority value
bool Directed_acyclic_graph::set_priority(int i, double p) {
	// Checks if the vertex value is within bounds
	if(i > vertices - 1) {
		throw illegal_argument();
	}
	if(i < 0) {
		throw illegal_argument();
	}
	
	if(hashTable->member(p)) { // Check if priority is already being used.
		return false;
	}
	hashTable->remove(priority[i]);
	hashTable->insert(p);
	priority[i] = p;
	resetPriority = true; // Reset priorities is now allowed since at least one has been changed.
	return true;
}

// Insert Edge Mutator: Check if an edge already exists between vertex i to j, 
// if not, make sure by inserting the edge a cycle is not created
bool Directed_acyclic_graph::insert_edge(int i, int j) {
	if(adjacent(i, j)) { // First check for already existing adjacency	
		return false;
	}
	if(connected(j, i)) { // Return false if it creates a cycle.
		return false;
	}
	dag[vertices*i + j] = 1;
	dag[vIn + j] = dag[vIn + j] + 1;
	dag[vOut + i] = dag[vOut + i] + 1;
	++edgeCount;
	return true;
}

// Clear Edges Mutator: Clears all the current edges between all the vertices
// Resets the number of edge count back to zero
void Directed_acyclic_graph::clear_edges() {
	// Check if edge count is zero -> no new connections were made
	// Therefore, we have nothing to clear
	if(edge_count() == 0) {
		return;
	}	
	// Check in-degree value for each vertex, if it is zero, skip over
	// the section (as it is already cleared), otherwise clear each entry
	for(int k = 0; k < vertices; ++k) {
		if(dag[vOut + k] > 0) {
			for(int j = 0; j < vertices; ++j){
				dag[vertices*k + j] = 0;
			}
			dag[vOut + k] = 0;
		}
		dag[vIn + k] = 0;
	}
	edgeCount = 0; // Reset edge count
}

// Reset Priorities Mutator: Resets the vertices' priorities back to their original values
void Directed_acyclic_graph::reset_priorities() {
	// Checks if a priority value has been changed, if not, return since unnecessary
	if(!resetPriority) {
		return;
	}
	// Reset all values to their original priority
	hashTable->clear();
	for(int k = 0; k < vertices; ++k) {
		priority[k] = k;
		hashTable->insert(k);
	}
	resetPriority = false; // Reset resetPriority flag
}

// Swap Mutator: Swaps all the member variables
void Directed_acyclic_graph::swap(Directed_acyclic_graph &dag2) {
	std::swap(resetPriority, dag2.resetPriority);
	std::swap(edgeCount, dag2.edgeCount);
	std::swap(vertices, dag2.vertices);	
	std::swap(vIn, dag2.vIn);
	std::swap(vOut, dag2.vOut);
	std::swap(vToDo, dag2.vToDo);
	std::swap(vDone, dag2.vDone);
	std::swap(vQueue, dag2.vQueue);
	std::swap(vIn2, dag2.vIn2);
	std::swap(length, dag2.length);
	std::swap(dag, dag2.dag);
	std::swap(priority, dag2.priority);
	std::swap(hashSize, dag2.hashSize);
	std::swap(hashTable, dag2.hashTable);
}

Directed_acyclic_graph &Directed_acyclic_graph::operator = (Directed_acyclic_graph rhs) {
	swap( rhs );
	return *this;
}

// Cout
std::ostream &operator << ( std::ostream &out, Directed_acyclic_graph const &list ) {
	out << "#";
	return out;
}

#endif