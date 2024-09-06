import numpy as np

from typing import Tuple


Point = Tuple[float, float]
Vector = Tuple[Point, Point]


def distance(point1: Point, point2: Point) -> float:
    """Compute the Euclidian distance between two points.

    Parameters
    ----------
    point1: Tuple[float, float]
    Point2: Tuple[float, float]

    Returns
    -------
    float
    """
    x1, y1 = point1
    x2, y2 = point2

    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)


def is_close(
    a: float,
    b: float,
    rel_tol: float = 1e-09,
    abs_tol: float = 0.0,
) -> bool:
    """Determine if the difference between two values falls within the given
    relative and absolute tolerances.

    Parameters
    ----------
    a: float
    b: float
    rel_tol: float, optional
        Relative tolerance. Default is 1e-09.
    abs_tol: float, 0.0
        Absolute tolerance. Default is 0.0.

    Returns
    -------
    bool
    """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def next_index(idx: int, n: int) -> int:
    """Get the next index.

    Index loops back to 0 if the next index exceeds the length of the
    list represented by n.

    Parameters
    ----------
    idx: int
        Current index value.
    n: int
        Length of the list or array.

    Returns
    -------
    int
    """
    return (idx + 1) % n
