#include "trie.hpp"

#include <iostream>

int main() {
  using wcoj::BinaryRelation;
  using wcoj::LeapfrogTrieJoin;

  BinaryRelation r{{1, 2}, {1, 3}, {2, 3}};
  BinaryRelation s{{2, 1}, {3, 1}, {3, 2}};
  BinaryRelation t{{1, 1}, {1, 2}, {2, 1}};

  const auto result = LeapfrogTrieJoin::triangle_join(r, s, t);

  std::cout << "Results: " << "\n";
  for (const auto &[a, b, c] : result) {
    std::cout << "(" << a << ", " << b << ", " << c << ")\n";
  }

  return 0;
}
