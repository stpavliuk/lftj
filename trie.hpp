#pragma once

#include <tuple>
#include <utility>
#include <vector>

namespace wcoj {

using BinaryRelation = std::vector<std::pair<int, int>>;
using TriangleTuple = std::tuple<int, int, int>;
using TriangleResult = std::vector<TriangleTuple>;
using NaryTuple = std::vector<int>;
using NaryJoinResult = std::vector<NaryTuple>;

struct NaryRelation {
  std::vector<int> variables;      // Variable ids used by this relation.
  std::vector<NaryTuple> tuples;   // One row per tuple, same width as variables.
};

class LeapfrogTrieJoin {
public:
  // Sort by (first, second) and remove duplicates.
  static BinaryRelation normalize_relation(BinaryRelation relation);

  // Distinct first-attribute keys, sorted.
  static std::vector<int> first_keys(const BinaryRelation &relation);

  // Group second-attribute values by first-attribute key.
  // Each value list is sorted and deduplicated.
  static std::vector<std::pair<int, std::vector<int>>>
  group_by_first(const BinaryRelation &relation);

  // Intersect two sorted vectors and return sorted unique values.
  static std::vector<int> intersect_sorted(const std::vector<int> &left,
                                           const std::vector<int> &right);

  // Cursor-based Leapfrog Triejoin for triangle query:
  // R(a,b), S(b,c), T(a,c) -> (a,b,c)
  static TriangleResult triangle_join(BinaryRelation r,
                                      BinaryRelation s,
                                      BinaryRelation t);

  // Generic cursor-based Leapfrog Triejoin over arbitrary variable order.
  // - `relations`: atoms with variable ids and tuples.
  // - `variable_order`: global variable order for LFTJ.
  // Returns tuples in `variable_order` column order.
  static NaryJoinResult nary_join(const std::vector<NaryRelation> &relations,
                                  const std::vector<int> &variable_order);
};

} // namespace wcoj
