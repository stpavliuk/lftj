#include "trie.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

std::vector<wcoj::NaryTuple> normalize_tuples(std::vector<wcoj::NaryTuple> tuples) {
    std::sort(tuples.begin(), tuples.end());
    tuples.erase(std::unique(tuples.begin(), tuples.end()), tuples.end());
    return tuples;
}

class GenericTrie {
public:
    struct Node {
        std::vector<int> keys;
        std::vector<std::size_t> child_ids;
    };

    explicit GenericTrie(std::vector<wcoj::NaryTuple> tuples, const std::size_t arity)
        : arity_(arity), tuples_(normalize_tuples(std::move(tuples))) {
        if (arity_ == 0 || tuples_.empty()) {
            return;
        }

        for (const auto &tuple: tuples_) {
            if (tuple.size() != arity_) {
                throw std::invalid_argument("Tuple width does not match relation arity");
            }
        }

        build_node(0, tuples_.size(), 0);
    }

    [[nodiscard]] bool empty() const {
        return nodes_.empty();
    }

    [[nodiscard]] const Node &node(const std::size_t node_id) const {
        return nodes_[node_id];
    }

private:
    std::size_t build_node(const std::size_t begin,
                           const std::size_t end,
                           const std::size_t depth) {
        const std::size_t node_id = nodes_.size();
        nodes_.push_back(Node{});

        std::size_t i = begin;
        while (i < end) {
            const int key = tuples_[i][depth];
            std::size_t j = i + 1;
            while (j < end && tuples_[j][depth] == key) {
                j += 1;
            }

            nodes_[node_id].keys.push_back(key);
            if (depth + 1 < arity_) {
                const auto child_id = build_node(i, j, depth + 1);
                nodes_[node_id].child_ids.push_back(child_id);
            }

            i = j;
        }

        return node_id;
    }

    std::size_t arity_;
    std::vector<wcoj::NaryTuple> tuples_;
    std::vector<Node> nodes_;
};

class GenericTrieCursor {
public:
    explicit GenericTrieCursor(const GenericTrie &trie) : trie_(&trie) {
        if (!trie_->empty()) {
            stack_.push_back(State{0, 0});
        }
    }

    [[nodiscard]] bool at_end() const {
        if (stack_.empty()) {
            return true;
        }
        const auto &state = stack_.back();
        const auto &keys = current_node().keys;
        return state.pos >= keys.size();
    }

    [[nodiscard]] int key() const {
        const auto &state = stack_.back();
        return current_node().keys[state.pos];
    }

    void next() {
        if (!at_end()) {
            stack_.back().pos += 1;
        }
    }

    void seek(const int target) {
        if (stack_.empty()) {
            return;
        }

        auto &[node_id, pos] = stack_.back();
        const auto &keys = current_node().keys;
        const auto begin = keys.begin() + static_cast<std::ptrdiff_t>(pos);
        const auto it = std::lower_bound(begin, keys.end(), target);
        pos = static_cast<std::size_t>(it - keys.begin());
    }

    void open() {
        if (stack_.empty() || at_end()) {
            return;
        }

        const auto &state = stack_.back();
        const auto &node = current_node();
        if (node.child_ids.empty()) {
            return;
        }

        const auto child_id = node.child_ids[state.pos];
        stack_.push_back(State{child_id, 0});
    }

    void up() {
        if (stack_.size() > 1) {
            stack_.pop_back();
        }
    }

private:
    struct State {
        std::size_t node_id;
        std::size_t pos;
    };

    [[nodiscard]] const GenericTrie::Node &current_node() const {
        return trie_->node(stack_.back().node_id);
    }

    const GenericTrie *trie_;
    std::vector<State> stack_;
};

template <typename OnMatch>
void for_each_intersection(std::vector<GenericTrieCursor *> &cursors, OnMatch on_match) {
    if (cursors.empty()) {
        return;
    }

    std::size_t pivot = 0;
    while (true) {
        int max_key = std::numeric_limits<int>::min();
        for (const auto *cursor: cursors) {
            if (cursor->at_end()) {
                return;
            }
            max_key = std::max(max_key, cursor->key());
        }

        bool moved = false;
        for (auto *cursor: cursors) {
            if (cursor->key() < max_key) {
                cursor->seek(max_key);
                if (cursor->at_end()) {
                    return;
                }
                moved = true;
            }
        }

        if (moved) {
            continue;
        }

        on_match(max_key);

        cursors[pivot]->next();
        if (cursors[pivot]->at_end()) {
            return;
        }
        pivot = (pivot + 1) % cursors.size();
    }
}

std::vector<wcoj::NaryTuple> binary_to_nary_tuples(const wcoj::BinaryRelation &relation) {
    std::vector<wcoj::NaryTuple> tuples;
    tuples.reserve(relation.size());
    for (const auto &[first, second]: relation) {
        tuples.push_back({first, second});
    }
    return tuples;
}

} // namespace

namespace wcoj {

BinaryRelation LeapfrogTrieJoin::normalize_relation(BinaryRelation relation) {
    std::sort(relation.begin(), relation.end());
    relation.erase(std::unique(relation.begin(), relation.end()), relation.end());
    return relation;
}

std::vector<int> LeapfrogTrieJoin::first_keys(const BinaryRelation &relation) {
    std::vector<int> keys;
    keys.reserve(relation.size());
    for (const auto &[first, second]: relation) {
        (void) second;
        keys.push_back(first);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

std::vector<std::pair<int, std::vector<int>>> LeapfrogTrieJoin::group_by_first(
    const BinaryRelation &relation) {
    std::vector<std::pair<int, std::vector<int>>> groups;

    for (const auto &[first, second]: relation) {
        if (groups.empty() || groups.back().first != first) {
            groups.push_back({first, {}});
        }
        groups.back().second.push_back(second);
    }

    for (auto &[key, values]: groups) {
        (void) key;
        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());
    }

    return groups;
}

std::vector<int> LeapfrogTrieJoin::intersect_sorted(const std::vector<int> &left,
                                                    const std::vector<int> &right) {
    if (left.empty() || right.empty()) {
        return {};
    }

    std::size_t left_idx = 0;
    std::size_t right_idx = 0;

    std::vector<int> result;
    while (left_idx < left.size() && right_idx < right.size()) {
        const int lv = left[left_idx];
        const int rv = right[right_idx];

        if (lv == rv) {
            if (result.empty() || result.back() != lv) {
                result.push_back(lv);
            }
            left_idx += 1;
            right_idx += 1;
        } else if (rv < lv) {
            right_idx = static_cast<std::size_t>(
                std::lower_bound(right.begin() + static_cast<std::ptrdiff_t>(right_idx),
                                 right.end(),
                                 lv) -
                right.begin());
        } else {
            left_idx = static_cast<std::size_t>(
                std::lower_bound(left.begin() + static_cast<std::ptrdiff_t>(left_idx),
                                 left.end(),
                                 rv) -
                left.begin());
        }
    }

    return result;
}

NaryJoinResult LeapfrogTrieJoin::nary_join(const std::vector<NaryRelation> &relations,
                                           const std::vector<int> &variable_order) {
    if (relations.empty() || variable_order.empty()) {
        return {};
    }

    std::unordered_set<int> seen_variables;
    seen_variables.reserve(variable_order.size());
    for (const int var: variable_order) {
        if (!seen_variables.insert(var).second) {
            throw std::invalid_argument("variable_order contains duplicate variable ids");
        }
    }

    for (const auto &relation: relations) {
        if (relation.tuples.empty()) {
            return {};
        }
    }

    std::vector<std::vector<int>> local_variables_by_relation;
    std::vector<std::vector<int>> local_depth_by_global_depth;
    local_variables_by_relation.reserve(relations.size());
    local_depth_by_global_depth.reserve(relations.size());

    std::vector<GenericTrie> tries;
    tries.reserve(relations.size());

    for (const auto &relation: relations) {
        if (relation.variables.empty()) {
            throw std::invalid_argument("relation.variables cannot be empty");
        }

        std::unordered_map<int, std::size_t> input_pos;
        input_pos.reserve(relation.variables.size());
        for (std::size_t i = 0; i < relation.variables.size(); ++i) {
            const int var = relation.variables[i];
            if (!input_pos.emplace(var, i).second) {
                throw std::invalid_argument("relation.variables contains duplicate variable ids");
            }
        }

        std::vector<int> local_variables;
        local_variables.reserve(relation.variables.size());

        std::vector<int> depth_map(variable_order.size(), -1);
        for (std::size_t global_depth = 0; global_depth < variable_order.size(); ++global_depth) {
            const int var = variable_order[global_depth];
            const auto it = input_pos.find(var);
            if (it == input_pos.end()) {
                continue;
            }

            depth_map[global_depth] = static_cast<int>(local_variables.size());
            local_variables.push_back(var);
        }

        if (local_variables.size() != relation.variables.size()) {
            throw std::invalid_argument("relation uses variables not present in variable_order");
        }

        std::vector<NaryTuple> local_tuples;
        local_tuples.reserve(relation.tuples.size());
        for (const auto &row: relation.tuples) {
            if (row.size() != relation.variables.size()) {
                throw std::invalid_argument("tuple width does not match relation.variables size");
            }

            NaryTuple reordered_row;
            reordered_row.reserve(local_variables.size());
            for (const int var: local_variables) {
                reordered_row.push_back(row[input_pos[var]]);
            }
            local_tuples.push_back(std::move(reordered_row));
        }

        tries.emplace_back(std::move(local_tuples), local_variables.size());
        local_variables_by_relation.push_back(std::move(local_variables));
        local_depth_by_global_depth.push_back(std::move(depth_map));
    }

    std::vector<GenericTrieCursor> cursors;
    cursors.reserve(tries.size());
    for (const auto &trie: tries) {
        cursors.emplace_back(trie);
    }

    std::vector<std::vector<std::size_t>> active_relations_by_depth(variable_order.size());
    std::vector<std::vector<std::size_t>> open_relations_by_depth(variable_order.size());

    for (std::size_t relation_idx = 0; relation_idx < relations.size(); ++relation_idx) {
        const auto &depth_map = local_depth_by_global_depth[relation_idx];
        const auto &local_vars = local_variables_by_relation[relation_idx];

        for (std::size_t global_depth = 0; global_depth < variable_order.size(); ++global_depth) {
            const int local_depth = depth_map[global_depth];
            if (local_depth < 0) {
                continue;
            }

            active_relations_by_depth[global_depth].push_back(relation_idx);
            if (static_cast<std::size_t>(local_depth + 1) < local_vars.size()) {
                open_relations_by_depth[global_depth].push_back(relation_idx);
            }
        }
    }

    for (std::size_t depth = 0; depth < variable_order.size(); ++depth) {
        if (active_relations_by_depth[depth].empty()) {
            throw std::invalid_argument("every variable in variable_order must appear in at least one relation");
        }
    }

    NaryJoinResult result;
    std::vector<int> binding(variable_order.size(), 0);

    std::function<void(std::size_t)> dfs = [&](const std::size_t depth) {
        if (depth == variable_order.size()) {
            result.push_back(binding);
            return;
        }

        const auto &active_relations = active_relations_by_depth[depth];

        std::vector<GenericTrieCursor *> active_cursors;
        active_cursors.reserve(active_relations.size());
        for (const auto relation_idx: active_relations) {
            active_cursors.push_back(&cursors[relation_idx]);
        }

        for_each_intersection(active_cursors, [&](const int matched_value) {
            const auto cursor_snapshot = cursors;
            binding[depth] = matched_value;

            const auto &relations_to_open = open_relations_by_depth[depth];
            for (const auto relation_idx: relations_to_open) {
                cursors[relation_idx].open();
            }

            dfs(depth + 1);

            for (std::size_t i = 0; i < cursors.size(); ++i) {
                cursors[i] = cursor_snapshot[i];
            }
        });
    };

    dfs(0);

    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

TriangleResult LeapfrogTrieJoin::triangle_join(BinaryRelation r,
                                               BinaryRelation s,
                                               BinaryRelation t) {
    const std::vector<NaryRelation> relations{
        NaryRelation{{0, 1}, binary_to_nary_tuples(normalize_relation(std::move(r)))} ,
        NaryRelation{{1, 2}, binary_to_nary_tuples(normalize_relation(std::move(s)))} ,
        NaryRelation{{0, 2}, binary_to_nary_tuples(normalize_relation(std::move(t)))} ,
    };

    const auto joined = nary_join(relations, {0, 1, 2});

    TriangleResult output;
    output.reserve(joined.size());
    for (const auto &row: joined) {
        if (row.size() != 3) {
            throw std::runtime_error("triangle_join expected 3-column rows from nary_join");
        }
        output.emplace_back(row[0], row[1], row[2]);
    }

    return output;
}

} // namespace wcoj
