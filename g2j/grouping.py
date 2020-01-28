"""Filter features."""


from itertools import groupby


def split_genic_features(features):
    genic, other = [], []
    for feature in features:
        if feature.is_genic():
            genic.append(feature)
        else:
            other.append(feature)
    return genic, other


def attributes_match(one, two):
    # Check Parent= is same as ID= of previous feature in proper GFF file
    try:
        return one.qualifiers["ID"] == two.qualifiers["Parent"]
    except KeyError:
        pass

    # Otherwise, look for common tags
    tags = {
        "proteinId",
        "protein_id",
        "transcriptId",
        "transcript_id",
        "name",
        "ID",
        "locus_tag",
    }
    common = set(one.qualifiers).intersection(two.qualifiers)

    if not common & tags:
        return False

    return any(one.qualifiers[key] == two.qualifiers[key] for key in common)


def group_overlapping_features(features, check_attributes=False):
    """Group overlapping sequence features.

    Given a sorted list of features, and format-specific overrides
    for self.features_overlap(), this generic iterator function will
    group sequence features via location overlap.

    This removes the dependence on shared identifier attributes
    (i.e. locus tags or protein IDs), which are unreliable.
    """
    if not features:
        return

    sorted_features = sorted(features, key=lambda f: f.location.min())
    first = sorted_features.pop(0)
    group, border = [first], first.location.max()

    for feature in sorted_features:
        if feature.location.min() <= border or (
            check_attributes and attributes_match(group[-1], feature)
        ):
            group.append(feature)
            border = max(border, feature.location.max())
        else:
            yield group
            group, border = [feature], feature.location.max()

    yield group


def collapse_same_type_features(features):
    """Combine features of the same type within an overlap group."""

    final = []
    feats = sorted(features, key=lambda f: f.type)

    for _, group in groupby(feats, key=lambda f: f.type):
        first, *group = list(group)

        for feature in group:
            first.location.start.extend(feature.location.start)
            first.location.end.extend(feature.location.end)
            first.qualifiers.update(feature.qualifiers)

        final.append(first)

    return sorted(final, key=lambda f: f.location.min())


def get_minimum(feature):
    """Get minimum of either a single Feature or Feature group.

    Facilitates sorting of features by location, regardless of whether they are grouped
    or not.
    """
    if hasattr(feature, "location"):
        return feature.location.min()
    return min(f.location.min() for f in feature)


def group_features(features):
    genic, other = split_genic_features(features)
    genic = [
        collapse_same_type_features(group)
        for group in group_overlapping_features(genic, check_attributes=True)
        if group
    ]
    return sorted(genic + other, key=lambda f: get_minimum(f))
