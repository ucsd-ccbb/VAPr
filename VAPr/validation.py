"""This module exposes utility functions to validate user inputs

By convention, validation functions in this module raise an appropriate Error if validation is unsuccessful.  If it is
successful, they return either nothing or the appropriately converted input value.
"""


__author__ = 'Birmingham'


def convert_to_nullable(input_val, cast_function):
    """For non-null input_val, apply cast_function and return result if successful; for null input_val, return None.

    Args:
        input_val (Any): The value to attempt to convert to either a None or the type specified by cast_function.
            The recognized null values are '.', None, '', and 'NULL'
        cast_function (Callable[[Any], Any]): A function to cast the input_val to some specified type; should raise an
            error if this cast fails.

    Returns:
        None if input is the null value.
        An appropriately cast value if input is not null and the cast is successful.

    Raises:
        Error: whatever error is provided by cast_function if the cast fails.
    """
    if input_val in ['.', None, '', 'NULL']:
        result = None
    else:
        result = cast_function(input_val)
    return result


def convert_to_nonneg_int(input_val, nullable=False):
    """For non-null input_val, cast to a non-negative integer and return result; for null input_val, return None.

    Args:
        input_val (Any): The value to attempt to convert to either a non-negative integer or a None (if nullable).
            The recognized null values are '.', None, '', and 'NULL'
        nullable (Optional[bool]): True if the input value may be null, false otherwise.  Defaults to False.

    Returns:
        None if nullable=True and the input is a null value.
        The appropriately cast non-negative integer if input is not null and the cast is successful.

    Raises:
        ValueError: if the input cannot be successfully converted to a non-negative integer or, if allowed, None
    """
    try:
        if nullable:
            result = convert_to_nullable(input_val, float)
            if result is None:
                return result
        else:
            result = float(input_val)

        if not result.is_integer(): raise ValueError()
        if result < 0: raise ValueError()
        return int(result)
    except ValueError:
        raise ValueError("Input ({0}) must be a non-negative integer".format(input_val))
