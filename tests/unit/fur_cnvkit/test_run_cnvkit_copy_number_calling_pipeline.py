import pytest

from fur_cnvkit.run_cnvkit_copy_number_calling_pipeline import parse_call_thresholds


def test_parse_call_thresholds_valid_string():
    thresholds = parse_call_thresholds("-1.1,-0.4,0.3,0.7")
    assert thresholds == [-1.1, -0.4, 0.3, 0.7]


def test_parse_call_thresholds_none_returns_none():
    assert parse_call_thresholds(None) is None


def test_parse_call_thresholds_requires_at_least_two_values():
    with pytest.raises(ValueError):
        parse_call_thresholds("0.5")


def test_parse_call_thresholds_validates_numeric_entries():
    with pytest.raises(ValueError):
        parse_call_thresholds("-1.0,foo,0.5")


def test_parse_call_thresholds_validates_strictly_increasing():
    with pytest.raises(ValueError):
        parse_call_thresholds("-1.0,-0.4,-0.4,0.5")
