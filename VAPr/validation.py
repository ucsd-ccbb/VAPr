"""This module exposes utility functions to validate user inputs

By convention, validation functions in this module raise an appropriate Error if validation is unsuccessful.  If it is
successful, they return either nothing or the appropriately converted input value.
"""

# personal libraries

__author__ = 'Birmingham'


def convert_to_nullable(input_val, cast_function):
    """For non-null input_val, apply cast_function and return result if successful; otherwise, return database null.

    Args:
        input_val (Any): The value to attempt to convert to either a null or the type specified by cast_function.
            The recognized null value is '.'
        cast_function (Callable[[Any], Any]): A function to cast the input_val to some specified type; should raise an
            error if this cast fails.

    Returns:
        utilities.database.NULL value if input is the null value.
        An appropriately cast value if input is not null and the cast is successful.

    Raises:
        Error: whatever error is provided by cast_function if the cast fails.
    """
    if input_val in ['.', None, '', 'NULL']:
        return -9999
    else:
        try:
            result = cast_function(input_val)
        except:
            result = -9999
    return result


def convert_to_nonneg_int(input_val, nullable=False):
    """Cast the input to a non-negative integer, if possible, or (optionally) a database null if matches null value.

    Args:
        input_val (Any): The value to attempt to convert to a non-negative integer (or a null, if nullable=True).
            The recognized null value is '.'
        nullable (Optional[bool]): True if the input value may be null, false otherwise.  Defaults to False.

    Returns:
        utilities.database.NULL value if nullable=True and the input is the null value.
        The appropriately cast non-negative integer if input is not null and the cast is successful.

    Raises:
        ValueError: if the input cannot be successfully converted to a non-negative integer (or null, if nullable=True)
    """
    err_count = 0
    if input_val == 'NULL':
        err_count += 1
        return -9999

    if nullable:
        try:
            result = convert_to_nullable(input_val, float)
        except:
            result = -9999
        if result == 'NULL':
            err_count += 1
            return -9999, err_count

    else:
        try:
            result = int(input_val)
        except TypeError:
            err_count += 1
        except ValueError:
            err_count += 1

    if not result.is_integer():
        err_count += 1
    if result < 0:
        err_count += 1
    return result, err_count
