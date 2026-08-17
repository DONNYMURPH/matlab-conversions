"""example pytest module"""

from python_template.example import example


def test_example_passthrough():
    """Test that example returns the given value unchanged."""
    assert example(True) is True


def test_example_with_different_types():
    """Test that example works with various types."""
    assert example(42) == 42
    assert example("hello") == "hello"
    assert example(None) is None
    assert example([1, 2, 3]) == [1, 2, 3]
