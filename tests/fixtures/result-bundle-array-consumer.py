"""Independent standard-library consumer of Shennong array metadata.

Paths contain native R selectors, not JSON Pointer components. JSON objects
preserve R list order; row-serialized data frames require column projection.
No Shennong code or numeric Python package is imported here.
"""

import json
import math
import sys


def resolve(root, path):
    current = root
    for selector in path:
        if isinstance(selector, bool) or not isinstance(selector, (str, int)):
            raise AssertionError(f"Invalid native R selector: {selector!r}")
        if isinstance(current, dict):
            if isinstance(selector, str):
                current = current[selector]
            else:
                assert 1 <= selector <= len(current), (path, selector)
                current = list(current.values())[selector - 1]
        elif isinstance(current, list):
            if isinstance(selector, int):
                assert 1 <= selector <= len(current), (path, selector)
                current = current[selector - 1]
            else:
                # A named column of a dataframe serialized as an array of rows.
                assert all(isinstance(row, dict) and selector in row for row in current), path
                current = [row[selector] for row in current]
        else:
            raise AssertionError(f"Cannot traverse {path!r} at {selector!r}")
    return current


def assert_shape(values, dimensions):
    assert isinstance(values, list), (dimensions, values)
    assert len(values) == dimensions[0], (dimensions, len(values))
    if len(dimensions) > 1:
        for value in values:
            assert_shape(value, dimensions[1:])
    else:
        assert all(not isinstance(value, (dict, list)) for value in values), dimensions


def assert_scalar(actual, expected, storage_type):
    if expected is None:
        assert actual is None, (actual, expected)
    elif storage_type in ("double", "integer"):
        assert not isinstance(actual, bool) and isinstance(actual, (int, float)), actual
        assert math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12), (actual, expected)
    else:
        assert type(actual) is type(expected) and actual == expected, (actual, expected)


def consume(bundle, expected):
    assert bundle["schema"] == "shennong.dev/analysis-result-bundle/v1"
    metadata = bundle["serialization"]
    assert metadata["schema"] == "shennong.dev/r-array-metadata/v1"
    records = metadata["arrays"]
    assert isinstance(records, list)
    keyed = {}
    for record in records:
        key = json.dumps(record["path"], ensure_ascii=False)
        assert key not in keyed, f"Duplicate metadata path: {key}"
        keyed[key] = record
    expected_keys = {json.dumps(record["path"], ensure_ascii=False) for record in expected["arrays"]}
    assert set(keyed) == expected_keys, (set(keyed) - expected_keys, expected_keys - set(keyed))
    checked = 0
    for reference in expected["arrays"]:
        record = keyed[json.dumps(reference["path"], ensure_ascii=False)]
        for field in ("storage_type", "dim", "dimnames", "axis_names"):
            assert record[field] == reference[field], (reference["path"], field, record[field], reference[field])
        values = resolve(bundle, record["path"])
        dimensions = record["dim"]
        assert isinstance(dimensions, list) and dimensions
        assert all(isinstance(dimension, int) and dimension >= 0 for dimension in dimensions)
        assert_shape(values, dimensions)
        assert len(reference["coordinates"]) == math.prod(dimensions)
        for coordinate in reference["coordinates"]:
            value = values
            labels = []
            for axis, position in enumerate(coordinate["index"]):
                assert 1 <= position <= dimensions[axis]
                value = value[position - 1]
                names = None if record["dimnames"] is None else record["dimnames"][axis]
                labels.append(None if names is None else names[position - 1])
            assert labels == coordinate["labels"], (record["path"], coordinate, labels)
            assert_scalar(value, coordinate["value"], record["storage_type"])
            checked += 1
    return len(records), checked


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as handle:
        bundle = json.load(handle)
    with open(sys.argv[2], encoding="utf-8") as handle:
        expected = json.load(handle)
    arrays, coordinates = consume(bundle, expected)
    print(f"Validated {arrays} arrays and {coordinates} labeled scalar coordinates")
