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
    """Find a matching attribute in two Features."""
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
        values = set()
        for tag in tags.intersection(f.qualifiers):
            value = f.qualifiers[tag]
            if isinstance(value, list):
                values.update(value)
            else:
                values.add(value)
        return values

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
        if (
            feature.location.min() <= border
            or check_attributes and attribute_match(group[-1], feature)
        ):
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


def collapse_features(features):
    """Collapses features that share the same ID.

    For example, CDS features will typically share the same ID= in a GFF file, whereas
    exons/introns/misc features will not. In this case, this function will collapse
    the CDS features (e.g. single multi-interval Feature) and ignore the others.
    """
    final = []
    features.sort(key=lambda f: f.qualifiers["ID"])
    for _, group in groupby(features, key=lambda f: f.qualifiers["ID"]):
        group = list(group)
        if len(group) == 1:
            final.extend(group)
            continue
        feature, *rest = group
        for other in rest:
            feature += other
        final.append(feature)
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
    """Get length of a feature."""
    if hasattr(feature, "location"):
        return len(feature.location)
    return get_maximum(feature) - get_minimum(feature)


def group_features(features, check_attributes=True):
    """Perform grouping on a collection of Feature objects."""
    genic, other = split_genic_features(features)
    results = []
    for group in group_overlapping_features(genic, check_attributes):
        if not group:
            continue
        subgroups = link_features(group)
        for subgroup in subgroups:
            results.append(subgroup)
    return sorted(results + other, key=lambda f: get_minimum(f))
