"""This module exposes unit-tests of validation.py functionality."""

# standard libraries
import unittest

# project libraries
import VAPr.validation as ns_test

__author__ = 'Birmingham'


class TestFunctions(unittest.TestCase):
    """Test top-level functions of validation.py module."""

    # region convert_to_nullable tests
    def test_convert_to_nullable_null(self):
        """Test that a null input ('.') produces output of None."""
        real_output = ns_test.convert_to_nullable(".", float)
        self.assertIsNone(real_output)

    def test_convert_to_nullable_not_null(self):
        """Test that a castable non-null input produces a cast output."""
        real_output = ns_test.convert_to_nullable("20.1", float)
        self.assertEqual(20.1, real_output)

    def test_convert_to_nullable_not_null_error(self):
        """Test that a non-castable non-null error raises an error."""
        with self.assertRaises(ValueError):
            ns_test.convert_to_nullable("20.B", float)

    # endregion

    # region convert_to_nonneg_int
    def test_convert_to_nonneg_int_nullable_pass_null(self):
        """Test that a null input produces output of utilities.database.NULL when nullable=True."""
        real_output = ns_test.convert_to_nonneg_int(".", nullable=True)
        self.assertIsNone(real_output)

    def test_convert_to_nonneg_int_nullable_fail(self):
        """Test that a null input raises an error when nullable=False."""
        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int(".")

    def test_convert_to_nonneg_int_not_nullable_pass(self):
        """Test that a non-null, non-negative-integer-castable input produces cast output for either nullable value."""
        real_output = ns_test.convert_to_nonneg_int("3.0")
        self.assertEqual(3, real_output)

        real_output = ns_test.convert_to_nonneg_int("3.0", nullable=True)
        self.assertEqual(3, real_output)

    def test_convert_to_nonneg_int_not_integer(self):
        """Test that a non-null, non-integer input raises an error for either nullable value."""
        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int("3.1")

        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int("3.1", nullable=True)

    def test_convert_to_nonneg_int_fail_negative(self):
        """Test that a non-null, negative integer input raises an error for either nullable value."""
        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int("-3.0")

        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int("-3.0", nullable=True)

    def test_convert_to_nonneg_int_fail_negative_non_integer(self):
        """Test that a non-null, negative non-integer input raises an error for either nullable value."""
        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int("-3.1")

        with self.assertRaises(ValueError):
            ns_test.convert_to_nonneg_int("-3.1", nullable=True)

    # endregion
