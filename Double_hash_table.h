// Author: Christopher Mitcheltree
// Date: Fall 2013

#ifndef DOUBLE_HASH_TABLE_H
#define DOUBLE_HASH_TABLE_H

enum state { EMPTY, OCCUPIED, DELETED };

template<typename Type>
class Double_hash_table {
	private:
		int count;
		int power; // Used for constructor.
		int array_size;
		int mask; // Used for constructor.
		Type *array;
		state *occupied;

		int h1( Type const & ) const;
		int h2( Type const & ) const;

	public:
		Double_hash_table( int = 5 );
		Double_hash_table( Double_hash_table const & );
		~Double_hash_table();
		int size() const;
		int capacity() const;
		double load_factor() const;
		bool empty() const;
		bool member( Type const & ) const;
		Type bin( int ) const;

		void print() const;

		void insert( Type const & );
		bool remove( Type const & );
		void clear();

	// Friends

	template <typename T>
	friend std::ostream &operator << ( std::ostream &, Double_hash_table<T> const & );
};

// Constructor: creates a new instance of the double hash table data structure.
template<typename Type>
Double_hash_table<Type>::Double_hash_table( int m ):
count( 0 ), power( m ),
array_size( 1 << power ),
mask( array_size - 1 ),
array( new Type[array_size] ),
occupied( new state[array_size] ) {
	for ( int i = 0; i < array_size; ++i ) {
		occupied[i] = EMPTY; // Initialize all bins as being empty.
	}
}

template<typename Type>
Double_hash_table<Type>::Double_hash_table(Double_hash_table<Type> const &obj):
	count( obj.count ), power( obj.power ),
	array_size( obj.array_size ),
	mask( obj.mask ),
	array( new Type[array_size] ),
	occupied( new state[array_size] ) {
		for(int i = 0; i < array_size; ++i) {
			occupied[i] = obj.occupied[i];
			array[i] = obj.array[i];
		}
}

// Deconstructor: deallocates any memory that was allocated by the double hash table constructor.
template<typename Type>
Double_hash_table<Type>::~Double_hash_table() {
	delete[] array;
	delete[] occupied;
}

// Hash Function 1: primary hash function: object is cast as an int and then this value
// is taken modulo the capacity of the hash table (which is then made positive if necessary).
// This function determines the bin the object is placed in.
template<typename Type>
int Double_hash_table<Type>::h1(Type const &obj) const {
	int hashValue1 = (static_cast<int>(obj*1000.0)) % (capacity());
	hashValue1 = (hashValue1 < 0) ? (hashValue1 + capacity()) : (hashValue1); // Make primary hash value positive by adding capacity if necessary.
	return hashValue1;
}

// Hash Function 2: secondary hash function: object is cast as an int, then divided by
// capacity and then this value is taken modulo the capacity of the hash table 
// (which is then made positive if necessary). Lastly this value is then made odd to
// ensure that the jump size is coprime with the capacity of the hash table.
// This function determines the jump size of the object.
template<typename Type>
int Double_hash_table<Type>::h2(Type const &obj) const {
	int hashValue2 = (static_cast<int>(obj*1000.0) / capacity()) % (capacity());
	hashValue2 = (hashValue2 < 0) ? (hashValue2 + capacity()) : (hashValue2); // Make secondary hash value positive by adding capacity if necessary.
	hashValue2 |= 1; // Bit wise OR with 1 in order to make the secondary hash value odd.
	return hashValue2;
}

// Size Accessor: returns how many elements are currently stored in the hash table.
template<typename Type>
int Double_hash_table<Type>::size() const {
	return count;
}

// Capacity Accessor: returns the total capacity of the hash table (how many bins there are).
template<typename Type>
int Double_hash_table<Type>::capacity() const {
	return array_size;
}

// Load Factor Accessor: returns how full the hash table is, that is the size divided by the capacity.
template<typename Type>
double Double_hash_table<Type>::load_factor() const {
	return (static_cast<double>(size()) / static_cast<double>(capacity()));
}

// Empty Accessor: returns true if there are no elements stored in the hash table.
template<typename Type>
bool Double_hash_table<Type>::empty() const {
	return(size() == 0);
}

// Member Accessor: returns true if the specified element is currently stored in the hash table.
template<typename Type>
bool Double_hash_table<Type>::member(Type const &obj) const {
	int index = h1(obj);

	for(int i = 0; i < capacity(); ++i) { // Worst case scenario: search entire hash table.
		if(occupied[index] == EMPTY) { // Found element.
			return false;
		}
		else if(array[index] == obj && occupied[index] == OCCUPIED) { // Found different element.
			return true;
		}
		else if(array[index] == obj && occupied[index] == DELETED) { // Found element but is flagged as erased.
			return false;
		}
		else {
			index = (index + h2(obj)) % (capacity()); // Search next bin determined by jump size (secondary hash value).
		}
	}
	return false; // Did not find element.
}

// Bin Accessor: returns the element being stored in a specific bin.
template<typename Type>
Type Double_hash_table<Type>::bin(int n) const {
	return array[n];
}

// Print Accessor: prints the contents of the entire hash table.
template<typename Type>
void Double_hash_table<Type>::print() const {
	for(int i = 0; i < capacity(); ++i) {
		std::cout << array[i] << " ";
	}
}

// Insert Mutator: inserts an element into the hash table.
template<typename Type>
void Double_hash_table<Type>::insert(Type const &obj) {
	if(load_factor() == 1.0) { 
		throw overflow(); // Prevents insertion when hash table is full.
	}
	if(member(obj)) {
		return; // Prevents duplicate elements from being inserted.
	}
	
	int index = h1(obj); // Primary hash value determines which bin to store element in.

	for(int i = 0; i < capacity(); ++i) { // Worst case scenario: go through entire hash table.
		if(occupied[index] == EMPTY || occupied[index] == DELETED) { // Insert in an empty or deleted bin.
			array[index] = obj;
			occupied[index] = OCCUPIED;
			break; // Insertion is completed once element has been successfully inserted.
		}
		else {
			index = (index + h2(obj)) % (capacity()); // Go to next bin determined by jump size (secondary hash value).
		}
	}
	++count;
}

// Remove Mutator: removes a specified element from the hash table.
template<typename Type>
bool Double_hash_table<Type>::remove(Type const &obj) {
	if(!member(obj)) { // Returns false if element is not in hash table.
		return false; 
	}

	int index = h1(obj); // Start looking for element in the bin determined by primary hash value.

	for(int i = 0; i < capacity(); ++i) { // Worst case scenario: go through entire hash table.
		if(array[index] == obj) { // Found element.
			occupied[index] = DELETED;
			break;
		}
		else {
			index = (index + h2(obj)) % (capacity()); // Go to next bin determined by jump size (secondary hash value).
		}
	}
	--count;
	return true;
}

// Clear Mutator: removes all elements from the hash table.
template<typename Type>
void Double_hash_table<Type>::clear() {
	for(int i = 0; i < capacity(); ++i) {
		occupied[i] = EMPTY; // Sets all bins of the hash table to empty.
	}
	count = 0;
}

// Cout
template <typename T>
std::ostream &operator << ( std::ostream &out, Double_hash_table<T> const &hash ) {
	for ( int i = 0; i < hash.capacity(); ++i ) {
		if ( hash.occupied[i] == EMPTY ) {
			out << "- ";
		} else if ( hash.occupied[i] == DELETED ) {
			out << "x ";
		} else {
			out << hash.array[i] << ' ';
		}
	}

	return out;
}

#endif
