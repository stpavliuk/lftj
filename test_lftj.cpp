#include "trie.hpp"

#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

using wcoj::BinaryRelation;
using wcoj::LeapfrogTrieJoin;
using wcoj::NaryJoinResult;
using wcoj::NaryRelation;
using wcoj::TriangleResult;

namespace {

void expect_equal(const std::vector<int> &actual, const std::vector<int> &expected) {
  assert(actual == expected);
}

void expect_equal(const TriangleResult &actual, const TriangleResult &expected) {
  assert(actual == expected);
}

void expect_equal(const NaryJoinResult &actual, const NaryJoinResult &expected) {
  assert(actual == expected);
}

void test_normalize_relation() {
  const BinaryRelation relation{{2, 1}, {1, 3}, {2, 1}, {1, 2}};
  const auto normalized = LeapfrogTrieJoin::normalize_relation(relation);
  const BinaryRelation expected{{1, 2}, {1, 3}, {2, 1}};

  assert(normalized == expected);
}

void test_intersect_sorted_basic() {
  const std::vector left{1, 2, 3, 4};
  const std::vector right{2, 4, 5};

  expect_equal(LeapfrogTrieJoin::intersect_sorted(left, right), {2, 4});
}

void test_intersect_sorted_with_duplicates() {
  const std::vector left{1, 2, 2, 3, 3, 4};
  const std::vector right{2, 2, 3, 5};

  expect_equal(LeapfrogTrieJoin::intersect_sorted(left, right), {2, 3});
}

void test_triangle_join_reference_example() {
  const BinaryRelation r{{1, 2}, {1, 3}, {2, 3}};
  const BinaryRelation s{{2, 1}, {3, 1}, {3, 2}};
  const BinaryRelation t{{1, 1}, {1, 2}, {2, 1}};

  const TriangleResult expected{
      std::make_tuple(1, 2, 1),
      std::make_tuple(1, 3, 1),
      std::make_tuple(1, 3, 2),
      std::make_tuple(2, 3, 1),
  };

  const auto actual = LeapfrogTrieJoin::triangle_join(r, s, t);

  expect_equal(actual, expected);
}

void test_triangle_join_deduplication() {
  const BinaryRelation r{{1, 2}, {1, 2}, {1, 3}};
  const BinaryRelation s{{2, 4}, {2, 4}, {3, 4}};
  const BinaryRelation t{{1, 4}, {1, 4}};

  const TriangleResult expected{
      std::make_tuple(1, 2, 4),
      std::make_tuple(1, 3, 4),
  };

  const auto actual = LeapfrogTrieJoin::triangle_join(r, s, t);

  expect_equal(actual, expected);
}

void test_triangle_join_empty_result() {
  const BinaryRelation r{{1, 2}};
  const BinaryRelation s{{3, 4}};
  const BinaryRelation t{{1, 4}};

  const auto actual = LeapfrogTrieJoin::triangle_join(r, s, t);

  assert(actual.empty());
}

void test_nary_join_triangle_equivalence() {
  const std::vector relations{
      NaryRelation{{0, 1}, {{1, 2}, {1, 3}, {2, 3}}},
      NaryRelation{{1, 2}, {{2, 1}, {3, 1}, {3, 2}}},
      NaryRelation{{0, 2}, {{1, 1}, {1, 2}, {2, 1}}},
  };

  const NaryJoinResult expected{
      {1, 2, 1},
      {1, 3, 1},
      {1, 3, 2},
      {2, 3, 1},
  };

  const auto actual = LeapfrogTrieJoin::nary_join(relations, {0, 1, 2});
  expect_equal(actual, expected);
}

void test_nary_join_four_variable_cycle() {
  const std::vector relations{
      NaryRelation{{0, 1}, {{1, 2}, {1, 3}, {2, 3}}},
      NaryRelation{{1, 2}, {{2, 4}, {3, 4}, {3, 5}}},
      NaryRelation{{2, 3}, {{4, 9}, {5, 8}}},
      NaryRelation{{0, 3}, {{1, 9}, {1, 8}, {2, 9}}},
  };

  const NaryJoinResult expected{
      {1, 2, 4, 9},
      {1, 3, 4, 9},
      {1, 3, 5, 8},
      {2, 3, 4, 9},
  };

  const auto actual = LeapfrogTrieJoin::nary_join(relations, {0, 1, 2, 3});
  expect_equal(actual, expected);
}

void test_nary_join_mixed_arity_relations() {
  const std::vector relations{
      NaryRelation{{0, 1, 2}, {{1, 2, 10}, {1, 3, 11}, {2, 3, 10}}},
      NaryRelation{{2, 3}, {{10, 7}, {11, 8}}},
      NaryRelation{{0, 3}, {{1, 7}, {1, 8}, {2, 7}}},
  };

  const NaryJoinResult expected{
      {1, 2, 10, 7},
      {1, 3, 11, 8},
      {2, 3, 10, 7},
  };

  const auto actual = LeapfrogTrieJoin::nary_join(relations, {0, 1, 2, 3});
  expect_equal(actual, expected);
}

} // namespace

int main() {
  test_normalize_relation();
  test_intersect_sorted_basic();
  test_intersect_sorted_with_duplicates();
  test_triangle_join_reference_example();
  test_triangle_join_deduplication();
  test_triangle_join_empty_result();
  test_nary_join_triangle_equivalence();
  test_nary_join_four_variable_cycle();
  test_nary_join_mixed_arity_relations();

  std::cout << "Tests passed" << std::endl;

  return 0;
}
