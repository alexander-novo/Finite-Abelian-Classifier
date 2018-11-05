#include "FiniteAbelianClassifier.h"

int main(int argc, char** argv) {
	// Arguments are expected to be given as a b c d ... mod z y x w ...
	if(argc % 2 == 1) {
		std::cerr << "Malformed arguments: Expected odd number of arguments" << std::endl;
		return 1;
	}
	if(argv[argc / 2] != std::string("mod")) {
		std::cerr << "Malformed arguments: Expected argument "<< (argc / 2)
		          << " to be \"mod\", instead found \"" << argv[argc/2]
		          << '"' << std::endl;
		return 1;
	}

	// Group given is G/H, where H = <g> is cyclic
	Group G(argc / 2 - 1);
	Tuple generator(argc / 2 - 1);
	// Keeps track of possibilities we want to consider
	std::vector<Group> possibleIsoGroups;
	// Keeps track of how many elements of each order there are in G
	std::map<unsigned,unsigned> GStats;
	// Keeps track of how many elements of each order in considered groups
	// Only one is needed, since after comparison with GStats,
	//     we have either found the answer or we know the group is wrong
	std::map<unsigned,unsigned> possibleStats;

	// Construct G from arguments and calculate order
	G.order = 1;
	for(int i = 1; i < argc / 2; i++) {
		G.products[i - 1] = atoi(argv[i]);
		G.order *= G.products[i - 1];
	}
	// This will be used to eliminate possibilities later
	G.findLargestOrderElement();

	// Construct generator of H, g
	for(int i = argc / 2 + 1; i <= argc; i++) {
		generator.x[i - argc / 2 - 1] = atoi(argv[i]);
	}

	std::cout << "Given group:\n"
	          << '\t' << G
	          << " / <" << generator << '>' << std::endl;

	std::cout << std::endl;

	// Find a list of possible isomorphic groups from FTFGAG
	// Some narrowing has already been done - specifically
	//     comparison of largest order elements, calculated earlier
	findPossibleIsoGroup(possibleIsoGroups, G, generator);

	std::cout << "Possible Isomorphic Groups (" << possibleIsoGroups.size() << "):\n";
	for(const auto& g : possibleIsoGroups) {
		std::cout << g << std::endl;
	}
	std::cout << std::endl;

	// Calculate the orders of the elements in G
	calcElementOrders(G, generator, GStats);

	// Then use this info to narrow down our possible suspects
	unsigned maxOrder = GStats.rbegin()->first;
	// Rule out any groups which don't have any elements on this order
	for(unsigned i = 0; i < possibleIsoGroups.size(); i++) {
		if(possibleIsoGroups[i].largestOrderElement < maxOrder) {
			possibleIsoGroups.erase(possibleIsoGroups.begin() + i);
			i--;
		}
	}

	std::cout << "For the given group:\n"
	          << printStats(GStats) << std::endl;

	std::cout << "Narrowed to (" << possibleIsoGroups.size() << "):\n";
	for(const auto& g : possibleIsoGroups) {
		std::cout << g << std::endl;
	}

	std::cout << std::endl;

	// Find the group by comparing orders of elements
	for(const auto& g : possibleIsoGroups) {
		possibleStats.clear();
		calcElementOrders(g, g.identity<unsigned>(), possibleStats);

		if(possibleStats == GStats) {
			std::cout << G << "/<" << generator << "> is isomorphic to " << g << std::endl;
			break;
		}
	}
	
	return 0;
}

// GENERAL GROUP FUNCTIONS ////////////////////////////////////////////////////////////////////////

void orderOfCoset(Coset& gh, const Coset& H, const Group& G) {
	Tuple conductor = gh.elements[0];
	Tuple& rep = gh.elements[0];
	gh.order = 1;

	// Simple algorithm - just keep adding the representative to itself until it is an element in H
	while(std::find(H.elements.begin(), H.elements.end(), conductor) == H.elements.end()) {
		gh.order++;
		conductor = G.mod(conductor + rep);
	}
}

void primeFactorize(unsigned n, std::map<unsigned,unsigned>& primes) {
	// Fairly simple algorithm
	// Starting at 2, divide out all factors of every number
	// Then advance by 2
	// Since we divide out primes as we come along them, we will never divide out any composite numbers

	// Divide out powers of two
	// Makes the next for loop slightly more efficient
	if(n % 2 == 0) {
		primes[2] = 0;
		do {
			primes[2]++;
			n /= 2;
		} while(n % 2 == 0);
	}

	// Only check odd numbers
	for(unsigned p = 3; p <= sqrt(n); p += 2) {
		if(n % p == 0) {
			primes[p] = 0;
			do {
				primes[p]++;
				n /= p;
			} while(n % p == 0);
		}
	}

	// The last number might be prime
	if(n != 1) {
		auto found = primes.find(n);
		if(found == primes.end()) {
			primes[n] = 1;
		} else {
			found->second++;
		}
	}
}

unsigned orderOfGenerator(const Group& G, const Tuple& generator) {
	if(G.products.size() != generator.n) throw std::logic_error("n-tuple size mismatch");

	unsigned currentOrder = 1;

	for(unsigned i = 0; i < G.products.size(); i++) {
		currentOrder = std::lcm(currentOrder, G.products[i] / std::gcd(G.products[i], generator.x[i]));
	}

	return currentOrder;
}

void findPossibleIsoGroup(std::vector<Group>& groups, const Group& G, const Tuple& generator) {
	// Partition Node is just a list of partitions, as well as the current partition we are looking at
	// Made a later algorithm easier
	struct partitionNode {
		std::vector<std::vector<unsigned>> partitions;
		unsigned iterator = 0;
	};

	// Map which stores partitions of powers of primes. Maps primes -> partitions of powers
	std::map<unsigned,partitionNode> partitions;
	// The current partition we are constructing, to be inserted into partitions
	std::vector<unsigned> currentPartition;
	// A map which maps primes to their powers. The primes are factors of the order of G/H
	std::map<unsigned,unsigned> factors;

	// Factor the order of G/H, to then partition
	unsigned orderOfGen = orderOfGenerator(G, generator);
	unsigned orderOfQgroup = G.order / orderOfGen;
	primeFactorize(orderOfQgroup, factors);

	// Multiplicatively partition the powers of primes by additively partitioning the powers
	unsigned power;
	int partitioner;
	unsigned rem_val;
	for(const auto& prime : factors) {
		power = prime.second;
		currentPartition.clear();
		currentPartition.push_back(power);
		
		partitioner = 0;

		while(true) {
			// Insert current partition
			partitions[prime.first].partitions.push_back(currentPartition);

			// Find rightmost non-1 value
			rem_val = 0;
			while(partitioner >= 0 && currentPartition[partitioner] == 1) {
				rem_val += 1;
				partitioner--;
			}

			if(partitioner < 0) break;

			currentPartition[partitioner]--;
			rem_val++;

			while(rem_val > currentPartition[partitioner]) {
				if(partitioner + 1 == currentPartition.size())
					currentPartition.push_back(currentPartition[partitioner]);
				else
					currentPartition[partitioner + 1] = currentPartition[partitioner];
				rem_val -= currentPartition[partitioner];
				partitioner++;
			}

			if(partitioner + 1 == currentPartition.size())
				currentPartition.push_back(rem_val);
			else {
				currentPartition.resize(partitioner + 2);
				currentPartition[partitioner + 1] = rem_val;
			}
			partitioner++;
		}
	}

	// Now mix every potential combination prime powers and we have our groups
	// Makes heavy use of iterator member of partitionNode above
	bool allDone = false;
	while(!allDone) {
		allDone = true;
		// Add a new potential group
		groups.push_back(Group(0));

		// The new group will be a sum of partitions of prime factors
		for(auto& part : partitions) {
			std::vector<unsigned>& partition = part.second.partitions[part.second.iterator];
			for(unsigned p : partition) {
				groups[groups.size() - 1].products.push_back(pow((long) part.first, (long) p));
			}

			if(allDone) part.second.iterator++;

			if(part.second.iterator == part.second.partitions.size()) {
				part.second.iterator = 0;
			} else {
				allDone = false;
			}
		}

		// Calculate order of largest element of new group
		// If it is too large - don't consider it at all
		groups[groups.size() - 1].findLargestOrderElement();
		if(G.largestOrderElement < groups[groups.size() - 1].largestOrderElement) {
			groups.erase(groups.end() - 1);
		}
	}
}

// Calculates orders of elements in G and removes any potential isomorphic candidates from the pool
void calcElementOrders(const Group& G,
	                   const Tuple& generator,
	                   std::map<unsigned,unsigned>& orderTracker) {
	Tuple identity = G.identity<unsigned>();
	Coset H;
	Coset gH;
	// unsigned ticket = 2;
	std::vector<Coset> qGroup;
	std::vector<Tuple> elements;
	Tuple element = generator;
	Tuple temp(G.products.size());

	G.generateElements<unsigned>(elements);

	// Generate identity coset
	H.elements.push_back(element);
	elements.erase(std::find(elements.begin(), elements.end(), element));

	while(element != identity) {
		element = G.mod(element + generator);
		H.elements.push_back(element);
		elements.erase(std::find(elements.begin(), elements.end(), element));
	}
	H.order = 1;

	// Add our identity coset to the group
	qGroup.push_back(H);
	orderTracker[1] = 1;

	// Find all cosets
	while(elements.size() > 0) {
		// Start by picking an element remaining in the group
		element = elements[0];

		// Add it to all of the elements in H to find gH
		gH.elements.clear();
		for(const Tuple& e : H.elements) {
			temp = G.mod(e + element);
			gH.elements.push_back(temp);
			elements.erase(std::find(elements.begin(), elements.end(), temp));
		}

		// Calculate the order of gH
		orderOfCoset(gH, H, G);
		qGroup.push_back(gH);

		// ticket++;
		if(orderTracker.find(gH.order) == orderTracker.end()) {
			orderTracker[gH.order] = 1;
		} else {
			orderTracker[gH.order]++;
		}
	}
}

std::string printStats(std::map<unsigned,unsigned>& stats) {
	std::stringstream ss;
	for(const auto& p : stats) {
		ss << "Elements of order " << p.first << ": " << p.second << '\n';
	}

	return ss.str();
}

// OUTPUT OPERATORS ///////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& out, const Tuple& e) {
	out << '(';
	for(unsigned i = 0; i < e.n; i++) {
		out << e.x[i] << ((i < e.n - 1) ? "," : "");
	}
	out << ')';

	return out;
}

std::ostream& operator<<(std::ostream& out, const Coset& e) {
	std::cout << "{";
	for(unsigned index = 0; index < e.elements.size(); index++) {
		std::cout << e.elements[index] << ((index < e.elements.size() - 1) ? "," : "");
	}
	std::cout << "}, order " << e.order;

	return out;
}

std::ostream& operator<<(std::ostream& out, const Group& g) {
	out << "(";
	for(unsigned i = 0; i < g.products.size(); i++) {
		out << "Z_" << g.products[i] << ((i < g.products.size() - 1) ? " x " : "");
	}
	out << ")";

	return out;
}