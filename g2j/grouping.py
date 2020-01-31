"""Filter features."""


from collections import defaultdict
from itertools import groupby, combinations


def split_genic_features(features):
    genic, other = [], []
    for feature in features:
        if feature.is_genic():
            genic.append(feature)
        else:
            other.append(feature)
    return genic, other


def attribute_match(one, two):
    """Find equal attribute between two Features."""
    tags = {
        "proteinId",
        "protein_id",
        "transcriptId",
        "transcript_id",
        "name",
        "Name",
        "ID",
        "locus_tag",
        "Parent",
        "protein_source_id",
    }
    common = tags & set(one.qualifiers).intersection(two.qualifiers)
    if not common:
        return False
    return any(one.qualifiers[key] == two.qualifiers[key] for key in common)


def link_features(features):
    """Find feature groups by single linkage."""
    tags = {
        "proteinId",
        "protein_id",
        "transcriptId",
        "transcript_id",
        "name",
        "Name",
        "ID",
        "locus_tag",
        "Parent",
    }

    def get_values(f):
        return set(f.qualifiers[t] for t in tags.intersection(f.qualifiers))

    def adjacency(features):
        """Generate adjacency list of Features.
        i.e. Pre-compute graph edges, so we can identify connected components
        using DFS/BFS
        """
        graph = defaultdict(set)
        total = len(features)
        values = [get_values(f) for f in features]
        for i, j in combinations(range(total), 2):
            if i != j and not values[i].isdisjoint(values[j]):
                graph[i].add(j)
                graph[j].add(i)
        return graph

    def dfs(graph, start, visited=None):
        """Depth-first search an adjacency list."""
        if visited is None:
            visited = set()
        visited.add(start)
        for next in graph[start] - visited:
            dfs(graph, next, visited)
        return visited

    # Build adjacency list of features
    graph, groups = adjacency(features), []

    # Depth first search; related nodes are marked as visited.
    # Save group, remove indices from the graph, repeat until
    # there is nothing left in the graph.
    while graph:
        group = dfs(graph, list(graph)[0])
        for index in group:
            graph.pop(index)
        groups.append([features[i] for i in group])

    return groups


def iter_overlapping_features(features, check_attributes=False):
    """Group overlapping sequence features."""
    if not features:
        return
    sorted_features = sorted(features, key=lambda f: f.location.min())
    first = sorted_features.pop(0)
    group, border = [first], first.location.max()

    for feature in sorted_features:
        within_border = feature.location.min() <= border
        shared_attrib = check_attributes and attribute_match(group[-1], feature)

        if within_border or shared_attrib:
            group.append(feature)
            border = max(border, feature.location.max())
        else:
            yield group
            group, border = [feature], feature.location.max()
    yield group


def group_overlapping_features(features, check_attributes=False):
    """Find grouped sequence features on forward and reverse strand."""
    forward, reverse = [], []
    for feature in features:
        if feature.location.strand == "+":
            forward.append(feature)
        else:
            reverse.append(feature)
    forward = list(iter_overlapping_features(forward, check_attributes))
    reverse = list(iter_overlapping_features(reverse, check_attributes))
    return forward + reverse


def iter_collapsed_features(features):
    """Find collapsible feature groups.
    Able to dig out overlapping same-type features this way.
    """
    result, previous = [], None
    for feature in features:
        if previous and attribute_match(previous, feature):
            result[-1].append(feature)
        else:
            result.append([feature])
        previous = feature
    return result


def collapse_same_type_features(features):
    """Combine features of the same type within an overlap group."""

    final = []
    feats = sorted(features, key=lambda f: f.type)

    for _, type_group in groupby(feats, key=lambda f: f.type):
        type_group = list(type_group)

        if len(type_group) == 1:
            final.extend(type_group)
            continue

        # To collapse, require at least one matching attribute
        # e.g. multiple CDS should share the same ID | Parent | protein_id
        for group in iter_collapsed_features(type_group):
            first, *rest = group
            for feature in rest:
                first += feature
            final.append(first)

    # Now, group different type Feature objects by attribute_link
    return final


def get_minimum(feature):
    """Get minimum of either a single Feature or Feature group."""
    if hasattr(feature, "location"):
        return feature.location.min()
    return min(f.location.min() for f in feature)


def get_maximum(feature):
    """Get maximum of either a single Feature or Feature group."""
    if hasattr(feature, "location"):
        return feature.location.max()
    return min(f.location.max() for f in feature)


def get_length(feature):
    if hasattr(feature, "location"):
        return len(feature.location)
    return get_maximum(feature) - get_minimum(feature)


def group_features(features, check_attributes=True):
    genic, other = split_genic_features(features)
    results = []
    for group in group_overlapping_features(genic, check_attributes):
        if not group:
            continue
        subgroups = link_features(group)
        for subgroup in subgroups:
            results.append(collapse_same_type_features(subgroup))
    return sorted(results + other, key=lambda f: get_minimum(f))
